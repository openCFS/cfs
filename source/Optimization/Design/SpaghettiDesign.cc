#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "Optimization/Design/SpaghettiDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Utils/PythonKernel.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Optimization/Function.hh"
#include "Optimization/TransferFunction.hh"

using std::string;
using std::pair;
using std::make_pair;
using std::to_string;

namespace CoupledField {

DEFINE_LOG(pasta, "pasta")

// The static Enum SpaghettiDesign::tip is in Optimization.cc to be available for DensityFile.cc when not having USE_EMBEDDED_PYTHON
// Once Spaghetti work also without Python, sort thy python stuff to SpaghettiDesignPython.cc and bring back tip.

SpaghettiDesign::SpaghettiDesign(StdVector<RegionIdType>& regionIds, PtrParamNode empn, ErsatzMaterial::Method method)
: FeaturedDesign(regionIds, empn, method)
{
  setup_timer_->Start();

  tip.SetName("SpaghettiDesign::Tip");
  tip.Add(START, "start");
  tip.Add(END, "end");
  tip.Add(INNER, "inner");

  PtrParamNode pn = empn->Get(ErsatzMaterial::method.ToString(method)); // spaghetti or spaghettiParamMat

  combine_ = combine.Parse(pn->Get("combine")->As<string>());
  boundary_ = boundary.Parse(pn->Get("boundary")->As<string>());
  if (method == ErsatzMaterial::SPAGHETTI_PARAM_MAT){
    orientation_ = orientation.Parse(pn->Get("orientation")->As<string>());
  }

  order = pn->Get("order")->As<unsigned int>();
  transition = pn->Get("transition")->As<double>(); // make optional
  radius = pn->Get("radius")->As<double>();
  if (pn->Has("q")){
    q = pn->Get("q")->As<double>();
  }
  if (pn->Has("p")){
    pen = pn->Get("p")->As<double>();
  }

  if(pn->Get("gradplot")->As<bool>()) // todo: move to FeaturedDesign
   gradplot_.open((progOpts->GetSimName() + ".grad.dat").c_str()); // the auto destructor does the job.

  // map_
  SetupMapping();

  // setup shaphetti with noodles and shape_param_ and opt_shape_param_
  SetupDesign(pn);

  // get rhomin from first density element and check transfer function
  assert(!data.IsEmpty());
  const DesignElement& de = data.First();
  rhomin = de.GetLowerBound(); // at least in python use for all
  rhomax = de.GetUpperBound();

  py_timer = info_->Get("spaghetti/python/timer")->AsTimer();
  py_timer->SetSub(); // we are in design_setup and in the eval timers

  if(pn->Has("python"))
    PythonInit(pn->Get("python"));

  assert(mapped_design_ != design_id);

  MapFeatureToDensity(); // only so late because of python -> calls PythonUpdateSpaghetti()

  setup_timer_->Stop(); // python and map are subtimers, so count them here
}

SpaghettiDesign::~SpaghettiDesign()
{
  Py_XDECREF(module_);
}

void SpaghettiDesign::ToInfo(ErsatzMaterial* em)
{
  AuxDesign::ToInfo(em);

  sp_info_ = info_->Get("spaghetti");

  sp_info_->Get("combine")->SetValue(combine.ToString(combine_));
  if(radius > 0.0)
    sp_info_->Get("radius")->SetValue(radius);
  sp_info_->Get("boundary/type")->SetValue(boundary.ToString(boundary_));
  sp_info_->Get("boundary/transition")->SetValue(transition);
  sp_info_->Get("orientation")->SetValue(orientation.ToString(orientation_));

  // python stuff written with PythonInit()

  PtrParamNode msh = sp_info_->Get("mesh");
  msh->Get("n")->SetValue(n_.ToString());

  // get keys of python glob.info_field -> to be repeated e.g. in MapFeatureGradient() as we have initially not much info
  StdVector<string> keys = PythonGetInfoFieldKeys();
  sp_info_->Get("info_field/genericElem")->SetValue(keys.ToString()); // TODO use nice Vector::ToString()

  // the stuff below will not be written to spaghetti
  PtrParamNode base = info_->Get("designVariables");
  for(Noodle& s : spaghetti)
    s.ToInfo(base->Get("noodle", ParamNode::APPEND));
}


StdVector<string> SpaghettiDesign::PythonGetInfoFieldKeys()
{
  // get keys of python glob.info_field
  // no need to time the stuff
  PyObject* afunc = PyObject_GetAttrString(module_, "cfs_info_field_keys");
  PythonKernel::CheckPythonFunction(afunc, "cfs_info_field_keys");

  PyObject* py_list = PyObject_CallObject(afunc, NULL);
  PythonKernel::CheckPythonReturn(py_list);

  if(!PyList_Check(py_list))
    EXCEPTION("python function cfs_info_field_keys did not return a list")

  StdVector<string> res(PyList_Size(py_list));

  for(unsigned int i = 0; i < res.GetSize(); i++)
  {
    PyObject* item = PyList_GetItem(py_list,i);
    const char* c_str = PyUnicode_AsUTF8(item);
    assert(c_str != NULL);
    res[i] = string(c_str);
    // we must not Py_XDECREF(item), it segfaults in the second run!
  }
  Py_XDECREF(py_list);
  Py_XDECREF(afunc);

  return res;
}


void SpaghettiDesign::PythonInit(PtrParamNode pn)
{
  // as spaghetti.py has no callback, we don't need to register
  py_timer->Start();

  pyopts = PythonKernel::ParseOptions(pn->GetList("option"));

  PythonKernel::LoadStatus stat = python->LoadPythonModule(pn);
  module_ = stat.module;

  sp_info_ = info_->Get("spaghetti");
  PtrParamNode py = sp_info_->Get("python");
  py->Get("file")->SetValue(stat.full_file);
  if(progOpts->DoDetailedInfo())
    py->Get("syspath")->SetValue(stat.sys_path.ToString(TS_PLAIN, ":"));
  py->Get("options")->SetValue(pyopts);

  // def cfs_init(rhomin, radius, boundary, transition, combine, orientation, nx, ny, nz, dx, design, dict):
  PyObject* func = PyObject_GetAttrString(module_, "cfs_init");
  PythonKernel::CheckPythonFunction(func, "cfs_init");

  PyObject* arg = PyTuple_New(4);

  Matrix<double> mm;
  domain->GetGrid()->CalcBoundingBoxOfRegion(GetRegionIds().First(), mm);
  string box_lower = "[";
  for(unsigned int i = 0; i < dim_; i++)
    box_lower += to_string(mm[i][0]) + ",";
  box_lower += "]";
  string box_upper = "[";
  for(unsigned int i = 0; i < dim_; i++)
    box_upper += to_string(mm[i][1]) + ",";
  box_upper += "]";

  StdVector<pair<string, string> > settings;
  settings.Push_back(make_pair("rhomin", to_string(rhomin)));
  settings.Push_back(make_pair("rhomax", to_string(rhomax)));
  if(radius > 0)
    settings.Push_back(make_pair("radius", to_string(radius)));
  settings.Push_back(make_pair("boundary", boundary.ToString(boundary_)));
  if(pen > 0)
    settings.Push_back(make_pair("p", to_string(pen)));
  settings.Push_back(make_pair("q", to_string(q)));
  settings.Push_back(make_pair("transition", to_string(transition)));
  settings.Push_back(make_pair("combine", combine.ToString(combine_)));
  settings.Push_back(make_pair("order", to_string(order)));
  settings.Push_back(make_pair("n", "[" + to_string(nx_) + "," + to_string(ny_) + "," + to_string(nz_) + "]"));
  settings.Push_back(make_pair("dx", to_string(dx_)));
  settings.Push_back(make_pair("box_lower", box_lower));
  settings.Push_back(make_pair("box_upper", box_upper));
  settings.Push_back(make_pair("orientation", orientation.ToString(orientation_)));
  PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonDict(settings));

  PyObject* des = PyTuple_New(design.GetSize());
  for(unsigned int i = 0; i < design.GetSize(); i++)
    PyTuple_SetItem(des, i,  PyUnicode_FromString(DesignElement::type.ToString(design[i].design).c_str()));
  PyTuple_SetItem(arg, 1, des); // steals the reference, so no need to decref

  PyObject* opt_ind = PyTuple_New(opt_indices.GetSize());
  for(unsigned int i = 0; i < opt_indices.GetSize(); i++)
    PyTuple_SetItem(opt_ind, i,  PyLong_FromLong(opt_indices[i]));
  PyTuple_SetItem(arg, 2, opt_ind); // steals the reference, so no need to decref

  PyTuple_SetItem(arg, 3, PythonKernel::CreatePythonDict(pyopts));

  PyObject* ret = PyObject_CallObject(func, arg);
  PythonKernel::CheckPythonReturn(ret);

  Py_XDECREF(ret);
  Py_XDECREF(arg);
  Py_XDECREF(func);

  py_timer->Stop();
}

PyObject* SpaghettiDesign::GetPythonModule()
{
  return module_;
}

void SpaghettiDesign::PythonUpdateSpaghetti()
{
  py_timer->Start();

  // def cfs_set_spaghetti(id, px, py, qx, qy, a_list, p):
  PyObject* func = PyObject_GetAttrString(module_, "cfs_set_spaghetti");
  PythonKernel::CheckPythonFunction(func, "cfs_set_spaghetti");

  for(Noodle& s : spaghetti)
  {
    // create a tuple for points, minimum two points
    PyObject* pyp = PyTuple_New(s.points.GetSize());
    for(unsigned pi = 0; pi < s.points.GetSize(); pi++)
    {
      StdVector<Variable>& p = s.points[pi];
      PyObject* lst = PyList_New(p.GetSize());
      for(unsigned int li = 0; li < p.GetSize(); li++)
        PyList_SetItem(lst, li, PyFloat_FromDouble(p[li].GetPlainDesignValue()));
      PyTuple_SetItem(pyp, pi, lst);
    }

    // create a tuple for a or r - is allowed to be empty
    const StdVector<Variable>& vec = s.type == Noodle::NORMAL ? s.a : s.r;
    PyObject* pyvec = PyTuple_New(vec.GetSize());
    for(unsigned int vi = 0; vi < vec.GetSize(); vi++)
      PyTuple_SetItem(pyvec, vi, PyFloat_FromDouble(vec[vi].GetPlainDesignValue()));

    PyObject* arg = PyTuple_New(s.HasAlpha() ? 5 : 4);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(s.idx));
    PyTuple_SetItem(arg, 1, pyp);
    PyTuple_SetItem(arg, 2, pyvec);
    PyTuple_SetItem(arg, 3, PyFloat_FromDouble(s.p.GetPlainDesignValue()));
    if(s.HasAlpha())
      PyTuple_SetItem(arg, 4, PyFloat_FromDouble(s.alpha.GetPlainDesignValue()));

    LOG_DBG2(pasta) << "PUS: s=" << s.idx << " vec=" << vec.GetSize() << " func=" << func << " arg=" << PyTuple_Size(arg) << " pyvec=" << pyvec << " alpha? " << s.HasAlpha();

    PyObject* ret = PyObject_CallObject(func, arg);
    PythonKernel::CheckPythonReturn(ret);

    Py_XDECREF(ret);
    Py_XDECREF(arg); // do NOT also Py_XDECREF(pya) as PyTuple_SetItem() steals the reference: https://docs.python.org/3/c-api/intro.html

    LOG_DBG(pasta) << "PUS sent " << s.ToString();
  }
  Py_XDECREF(func);

  py_timer->Stop();
}

void SpaghettiDesign::MapFeatureToDensity()
{
  mapping_timer_->Start();
  // this is the pure python implementation
  assert(module_);

  LOG_DBG(pasta) << "MFTD called";

  if(design_id != mapped_design_)
    PythonUpdateSpaghetti(); // has it's own timer, it would top our

  py_timer->Start();

  // up to now we only can do python
  PyObject* func = PyObject_GetAttrString(module_, "cfs_map_to_design");
  PythonKernel::CheckPythonFunction(func, "cfs_map_to_design");

  PyObject* ret = PyObject_CallObject(func, NULL);
  PythonKernel::CheckPythonReturn(ret);

  // more convenient to use a temporary object than to iterate directly on the numpy object
  Vector<double> tmp(ret,true); // decref ret
  Py_XDECREF(func);

  assert(!map_.IsEmpty());
  assert(tmp.GetSize() >= map_.GetSize()); // map_ can be smaller if the domain is not convex

  for(auto& item : map_)
    item.elemval->SetDesign(tmp[item.lexicographic_pos]);

  LOG_DBG(pasta)  << "MFTD: old mid=" << mapped_design_ << " current did=" << design_id << " size=" << tmp.GetSize();
  LOG_DBG2(pasta) << "MFTD: rho=" << tmp.ToString();

  mapped_design_ = design_id;

  py_timer->Stop();
  mapping_timer_->Stop();
}

void SpaghettiDesign::MapFeatureGradient(const Function* f)
{
  gradient_timer_->Start();

  assert(module_);

  assert(mapped_design_ == design_id);

  py_timer->Start();

  // python needs to be smart by itself to speedup gradient calculation because we will call this
  // for several functions for the same design set (compliance, volume, ...)

  // we cannot create a numpy array in C - so get it from python to be filled with the rho-density to be given to python
  PyObject* afunc = PyObject_GetAttrString(module_, "cfs_get_drho_vector");
  PythonKernel::CheckPythonFunction(afunc, "cfs_get_drho_vector");

  PyObject* np_rho = PyObject_CallObject(afunc, NULL);
  PythonKernel::CheckPythonReturn(np_rho);

  // we create a temporary vector for the gradient values such that we can use Export() - just for convenience.
  // With Python we expect no hpc performance anyway
  //Vector<double> drho(nx_*ny_*nz_);
  Vector<double> drho = Vector<double>(nx_*ny_*nz_*design.GetSize(), 0);
  for(auto& item : map_)
    drho[item.lexicographic_pos] = item.elemval->GetPlainGradient(f);

  // copy data to numpy array - does extensive range checks with good error message
  drho.Export(np_rho);

  // now get the gradient from python. We get all gradients and ignore efficiency optimization by skipping fixed variables
  PyObject* gfunc = PyObject_GetAttrString(module_, "cfs_get_gradient");
  PythonKernel::CheckPythonFunction(gfunc, "cfs_get_gradient");

  PyObject* arg = PyTuple_New(2);
  Py_INCREF(np_rho); // this is important, as PyTuple_SetItem() steals the reference: https://docs.python.org/3/c-api/intro.html
  PyTuple_SetItem(arg, 0, np_rho);
  PyTuple_SetItem(arg, 1, PyUnicode_FromString(f->ToString().c_str()));

  PyObject* np_grad = PyObject_CallObject(gfunc, arg);
  PythonKernel::CheckPythonReturn(np_grad);

  // create temporary vector - ordering is all nodes, profiles, normals as in density.xml and SetupDesign()
  // we set all to shape_param_ and opt_shape_param_ "gets updated automatically" - clearly its all pointer referencing
  Vector<double> dshape(np_grad, true); // decref!

  LOG_DBG2(pasta) << "MFG: " << f->ToString() << " -> " << dshape.ToString();

  assert(shape_param_.First()->GetIndex() == 0);
  assert(opt_shape_param_.First()->GetPlainGradient(f) == 0.0); // we expect it to be reset before

  if(dshape.GetSize() != shape_param_.GetSize())
    EXCEPTION("Got gradient of size " << dshape.GetSize() << " from Python, expected " << shape_param_.GetSize())

  for(ShapeParamElement* spe : shape_param_)
    spe->AddGradient(f, dshape[spe->GetIndex()]);

  Py_XDECREF(arg); // np_grad is already decref
  Py_XDECREF(np_rho);
  Py_XDECREF(gfunc);
  Py_XDECREF(afunc);

  py_timer->Stop();
  gradient_timer_->Stop();

  // on each iteration we might have new data to be reported
  StdVector<string> keys = PythonGetInfoFieldKeys();
  sp_info_->Get("info_field/genericElem")->SetValue(keys.ToString());
}


void SpaghettiDesign::PrepareSpecialResults()
{
  StdVector<const ResultDescription*> res = GetGenericResults();
  for(const ResultDescription* rs : res)
  {
    PyObject* afunc = PyObject_GetAttrString(module_, "cfs_get_info_field");
    PythonKernel::CheckPythonFunction(afunc, "cfs_get_info_field");

    PyObject* arg = PyTuple_New(1);
    PyTuple_SetItem(arg, 0, PyUnicode_FromString(rs->generic.c_str()));

    PyObject* np = PyObject_CallObject(afunc, arg);
    PythonKernel::CheckPythonReturn(np);

    Vector<double> sr(np, true); // decref
    assert(sr.GetSize() == data.GetSize()); // will fail once we have angle!!

    int ri = DesignElement::GetOptResultIndex(rs->solutionType);

    for(unsigned int i = 0; i < sr.GetSize(); i++)
      data[i].specialResult[ri] = sr[i];

    Py_XDECREF(arg); // np_grad is already decref
    Py_XDECREF(afunc);
  }
}


void SpaghettiDesign::AddVariable(Variable* var)
{
  // we do exact space "estimation"

  assert(shape_param_.HasSpace());
  var->SetIndex(shape_param_.GetSize());
  shape_param_.Push_back(var);

  if(!var->fixed)
  {
    assert(opt_shape_param_.HasSpace());
    opt_indices[opt_shape_param_.GetSize()] = shape_param_.GetSize()-1;
    var->SetOptIndex(opt_shape_param_.GetSize());
    opt_shape_param_.Push_back(var);
  }
}


void SpaghettiDesign::SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);
  assert(f->GetType() == Function::DISTANCE || f->GetType() == Function::BENDING);

  // we assume fixed only for profile if at all
  assert(locality == Function::Local::FUNCTION_SPECIFIC || locality == Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS);

  vem.Reserve(spaghetti.GetSize()); // assume nothing fixed

  assert(dim_ == 2);
  StdVector<BaseDesignElement*> nodes;

  bool two_signs = locality == Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS;

  int sign_1 = two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = -1;

  for(Noodle& s : spaghetti)
  {
    nodes.Clear();
    StdVector<Variable>& P = s.points.First();
    StdVector<Variable>& Q = s.points.Last();

    // assume nothing fixed
    if(P[0].fixed || P[1].fixed || Q[0].fixed || Q[1].fixed)
    {
      if (f->GetType() == Function::DISTANCE)
        throw Exception("distance constraints currently only for non-fixed nodes");
      // else: Bending
      if (s.IsFixed(P) && s.IsFixed(Q) && (s.a.GetSize() == 0))
        continue; // won't add empty constraint if all points are fixed -> next noodle
    }

    // px is element, then py, then qx then qy
    nodes.Push_back(&P[1]);
    nodes.Push_back(&Q[0]);
    nodes.Push_back(&Q[1]);

    if(f->GetType() == Function::DISTANCE)
    {
      if(f->GetDesignType() != BaseDesignElement::NODE)
        throw Exception("Configuration error: distance needs to be defined for node.");

      vem.Push_back(Function::Local::Identifier(&P[0], nodes));
      LOG_DBG(pasta) << "SVSEM: f=" << f->ToString() << " s=" << s.idx << " id=" << vem.Last().ToString();
    }
    else // Bending
    {
      // with  many segments, the typical case would be ai-1,ai,ai+1 (N N N)
      // but we have the cases for only two segments (0 N 0)
      // or for start ( 0 N N) or end ( N N 0)

      for(unsigned ai = 0; ai < s.a.GetSize(); ai++)
      {
        StdVector<BaseDesignElement*> buddies(nodes);
        assert(buddies.GetSize() == nodes.GetSize());
        assert(buddies[0] == nodes[0]); // no empty nodes for bending

        if(s.a.GetSize() == 1)
        {
          buddies.Push_back(&s.a[0]);
          vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
          vem.Last().bending = Function::Local::Identifier::ZNZ;
          if(two_signs){
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
            vem.Last().bending = Function::Local::Identifier::ZNZ;
          }
        }
        else
        {
          if(ai == 0)
          {
            buddies.Push_back(&s.a[0]);
            buddies.Push_back(&s.a[1]);
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::ZNN;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::ZNN;
            }
          }
          else if(ai == s.a.GetSize()-1)
          {
            buddies.Push_back(&s.a[ai-1]);
            buddies.Push_back(&s.a[ai]);
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::NNZ;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::NNZ;
            }
          }
          else
          {
            assert(ai > 0 && ai < s.a.GetSize()-1);
            buddies.Push_back(&s.a[ai-1]);
            buddies.Push_back(&s.a[ai]);
            buddies.Push_back(&s.a[ai+1]);
           // vem.Push_back(Function::Local::Identifier(&s.px, buddies));
            //vem.Last().bending = Function::Local::Identifier::NNN;
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::NNN;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::NNN;
            }
          }
        }
        LOG_DBG(pasta) << "SVSEM: f=" << f->ToString() << " s=" << s.idx << " id=" << vem.Last().ToString();
      } // end ai
    } // end bending
  } // end noodle
}


int SpaghettiDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  // practically updates the Variable objects within spaghetti. Set's design_id
  FeaturedDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent);

  if(mapped_design_ != design_id){
    PythonUpdateSpaghetti();
  }

  LOG_DBG(pasta) << "RDFE mid=" << mapped_design_ << " did=" << design_id;

  if(mapped_design_ != design_id)
    MapFeatureToDensity();

  return design_id;
}

void SpaghettiDesign::ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation)
{
  // we expect all nodes (start x, start y, end x, end y, ...) then all profiles then all normals (if any).
  // by the shape attribute we could alll collect from random order but for easy we expect this ordering
  assert(spaghetti.GetSize() > 0);

  ParamNodeList list = set->GetList("shapeParamElement");
  if(list.GetSize() != shape_param_.GetSize())
    EXCEPTION("expect " << shape_param_.GetSize() << " 'shapeParamElement' in density.xml but got " << list.GetSize());

  for(unsigned int i = 0; i < shape_param_.GetSize(); i++)
  {
    Variable*    var = dynamic_cast<Variable*>(shape_param_[i]);
    PtrParamNode pn  = list[i];

    // do a lot of validation
    const std::string pnnr = pn->Get("nr")->As<string>();
    if(pn->Get("nr")->As<unsigned int>() != i)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'nr' value " << i);

    if(pn->Get("shape")->As<int>() != var->noodle)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'shape' value " << var->noodle);

    if(DesignElement::type.Parse(pn->Get("type")->As<string>()) != var->GetType())
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'type' value " << DesignElement::type.ToString(var->GetType()));

    if(var->GetType() == Variable::NODE)
    {
      if(ShapeParamElement::dof.Parse(pn->Get("dof")->As<string>()) != var->dof_)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'dof' value " << ShapeParamElement::dof.ToString(var->dof_));
      if(SpaghettiDesign::tip.Parse(pn->Get("tip")->As<string>()) != var->tip)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'tip' value " << SpaghettiDesign::tip.ToString(var->tip));
    }

    var->SetDesign(pn->Get("design")->As<double>());
  }

  // it shall be save to update the design_id, because we changed stuff
  design_id++;

  MapFeatureToDensity();
}

void SpaghettiDesign::AddToDensityHeader(PtrParamNode header)
{
  PtrParamNode pn = header->Get("spaghetti");
  pn->Get("radius")->SetValue(radius);
}

void SpaghettiDesign::SetupDesign(PtrParamNode base)
{
  StdVector<PtrParamNode> pnl = base->GetList("noodle");
  spaghetti.Resize(pnl.GetSize());
  int total = 0;
  int opt   = 0;
  for(unsigned int i = 0; i < pnl.GetSize(); i++)
  {
    Noodle& noodle = spaghetti[i];
    noodle.Parse(pnl[i], i);
    total += noodle.GetTotalVariables();
    opt   += noodle.GetOptVariables();
  }
  // spaghetti contain the variable instances, now flatten the design

  shape_param_.Reserve(total);
  opt_shape_param_.Reserve(opt);
  opt_indices.Resize(opt);

  // the order is by noodles: nodes, profile, [alpha|, normals/radii
  // reading shall be transparent due to shape (noodle) idx.
  for(Noodle& noodle : spaghetti)
  {
    for(auto& p : noodle.points) // P, [inner], Q
      for(auto& v : p)
        AddVariable(&v);

    assert(noodle.a.IsEmpty() || noodle.r.IsEmpty()); // or both are empty
    for(Variable& var : noodle.r)
      AddVariable(&var);

    AddVariable(&(noodle.p));

    if(noodle.HasAlpha())
      AddVariable(&(noodle.alpha));

    for(Variable& var : noodle.a)
      AddVariable(&var);
  }

  // we did exact space "estimation"
  assert(shape_param_.GetCapacity() == shape_param_.GetSize());
  assert(opt_shape_param_.GetCapacity() == opt_shape_param_.GetSize());
}

SpaghettiDesign::Noodle::Noodle()
{
  type = dim_ == 2 ? NORMAL : POINTS;
}


void SpaghettiDesign::Noodle::Parse(PtrParamNode noodle, int idx)
{
  this->idx = idx;
  opt_variables_ = 0;

  int n = noodle->Get("segments")->As<int>();
  if(n < 1)
    throw Exception("'noodle' requires at least 1 segment");

  assert(points.IsEmpty());
  points.Resize(dim_ == 2 ? 2 : n + 1);

  // add tip
  StdVector<Variable>& start = points.First();
  start.Resize(dim_);
  start[0].Parse(noodle->GetByVal("node", "dof", "x", "tip", "start"),idx);
  start[1].Parse(noodle->GetByVal("node", "dof", "y", "tip", "start"),idx);
  if(dim_ == 3)
    start[2].Parse(noodle->GetByVal("node", "dof", "z", "tip", "start"),idx);
  opt_variables_ += CountNonFixed(start);

  StdVector<Variable>& end = points.Last();
  end.Resize(dim_);
  end[0].Parse(noodle->GetByVal("node", "dof", "x", "tip", "end"),idx);
  end[1].Parse(noodle->GetByVal("node", "dof", "y", "tip", "end"),idx);
  if(dim_ == 3)
    end[2].Parse(noodle->GetByVal("node", "dof", "z", "tip", "end"),idx);
  opt_variables_ += CountNonFixed(end);

  if(n == 1 && noodle->Has("inner"))
    throw Exception("'noodle' has defined 1 segment but 'inner' is given");

  if(dim_ == 3)
  {
    Vector<double> P = Variable::AsVector(start);
    Vector<double> Q = Variable::AsVector(end);
    for(int pi = 1; pi < n; pi++) // n = 1 -> no inner point; n = 2 -> 1 inner point 1; n = 3 -> 2 inner points 1,2
    {
      StdVector<Variable>& point = points[pi];
      point.Resize(dim_);

      for(unsigned int d = 0; d < dim_; d++)
      {
        double inter = Interpolate(P[d],Q[d],pi,points.GetSize());
        string dof = ShapeParamElement::dof.ToString((ShapeParamElement::Dof) d);
        // either we have inner or require clear lower and upper from P and Q
        if(noodle->HasByVal("inner", "dof", dof)) // either interpolate for value or math parser stuff
        {
          point[d].Parse(noodle->GetByVal("inner", "dof", dof), this->idx, inter);
          point[d].tip = INNER;
        }
        else
        {
          if(start[d].GetLowerBound() != end[d].GetLowerBound() || start[d].fixed || end[d].fixed)
            throw Exception("As bounds for tips do not match or are fixed, 'inner' needs to be given for dof " + dof);
          point[d].InitInnerNode(idx, (ShapeParamElement::Dof) d, inter, start[d].GetLowerBound(), start[d].GetUpperBound());
        }
      }

      opt_variables_ += dim_;
    }
  }

  p.Parse(noodle->Get("profile"),idx); // exactly one
  opt_variables_ += p.fixed ? 0 : 1;

  // number of normals is segments-1
  if(type == NORMAL)
  {
    if(noodle->Has("radius"))
      throw Exception("no 'radius' for this 'noodle' type expected (to be implemented :))");
    if(n == 1 && noodle->Has("normal"))
      throw Exception("'noodle' has defined 1 segment but 'normal' is given");
    if(n >= 2 && ! noodle->Has("normal"))
      throw Exception("'noodle' has defined " + std::to_string(n) + " segments but no 'normal' given");

    a.Resize(n-1);
    for(int i = 0; i < n-1; i++)
      a[i].Parse(noodle->Get("normal"),this->idx);

    opt_variables_ += !a.IsEmpty() && !a[0].fixed ? a.GetSize() : 0;
  }
  else
  {
    assert(type == POINTS);
    if(noodle->Has("normal"))
      throw Exception("no 'normal' for this 'noodle' type expected");

    if(n == 1 && noodle->Has("radius"))
      throw Exception("'noodle' has defined 1 segment but 'radius' is given");
    if(n >= 2 && ! noodle->Has("radius"))
      throw Exception("'noodle' has defined " + std::to_string(n) + " segments but no 'radius' given");

    r.Resize(n-1);
    for(int i = 0; i < n-1; i++)
      r[i].Parse(noodle->Get("radius"),this->idx);

    opt_variables_ += !r.IsEmpty() && !r[0].fixed ? r.GetSize() : 0;

    alpha.Parse(noodle->Get("alpha"),idx); // exactly one
    opt_variables_ += alpha.fixed ? 0 : 1;
  }

  LOG_DBG(pasta) << "N:start ov=" << opt_variables_ << " " << ToString();
}

bool SpaghettiDesign::Noodle::IsFixed(const StdVector<Variable>& point) const
{
  assert(point.GetSize() == dim_);
  for(const Variable& v : point)
    if(v.fixed)
      return true;
  return false;
}

int SpaghettiDesign::Noodle::CountNonFixed(const StdVector<Variable>& point) const
{
  int sum = 0;
  assert(point.GetSize() == dim_);
  for(const Variable& v : point)
    sum += v.fixed ? 0 : 1;
  return sum;
}

double SpaghettiDesign::Noodle::Interpolate(double start, double end, unsigned int idx, unsigned int n_segments)
{
  return start + (end-start)/(n_segments-1) * idx;
}

int SpaghettiDesign::Noodle::GetTotalVariables() const
{
  int vars = 1; // profile p
  assert(a.IsEmpty() || r.IsEmpty());
  vars += dim_ * points.GetSize() + a.GetSize() + r.GetSize();
  if(HasAlpha())
    vars += 1;

  return vars;
}


string SpaghettiDesign::Noodle::ToString()
{
  std::stringstream ss;

  ss << "idx=" << idx << " type=" << (type == NORMAL ? "normal" : "points");
  if(points.GetSize() > 1)
  {
    ss << " start=" << Variable::ToString(points.First());
    ss << " end=" << Variable::ToString(points.Last());
    for(unsigned int i = 1; i < points.GetSize()-1; i++)
      ss << " num_" << i << "=" << Variable::ToString(points[i]);
  }
  ss << " p=" << (p.fixed ? "*" : "") << p.GetPlainDesignValue();
  ss << " a=" << Variable::ToString(a);
  ss << " r=" << Variable::ToString(r);
  return ss.str();
}

void SpaghettiDesign::Noodle::ToInfo(PtrParamNode in)
{
  in->Get("id")->SetValue(idx);
  in->Get("segments")->SetValue(a.GetSize() + 1);

  // one could combine P and Q as coordinates.
  auto& P = points.First();
  auto& Q = points.Last();

  for(unsigned int d = 0; d < dim_; d++)
    P[d].ToInfo(in->Get("p" + Dof(d), ParamNode::APPEND));

  if(HasInner())
    for(unsigned int d = 0; d < dim_; d++)
      ToInfo(in, Variable::NODE, "i" + Dof(d), d);

  for(unsigned int d = 0; d < dim_; d++)
    Q[d].ToInfo(in->Get("p" + Dof(d), ParamNode::APPEND));

  p.ToInfo(in->Get("p", ParamNode::APPEND));

  if(HasAlpha())
    alpha.ToInfo(in->Get("alpha", ParamNode::APPEND));

  // will do nothing for empty vectors
  ToInfo(in, Variable::NORMAL, "a");
  ToInfo(in, Variable::RADIUS, "r");
}


void SpaghettiDesign::Noodle::CompareToInfoHelper(const Variable& v0, const Variable& test, std::string& fixed, std::string& lower, std::string& upper) const
{
  assert(v0.fixed == test.fixed);

  if(v0.GetPlainDesignValue() == test.GetPlainDesignValue()) {
    if(fixed != "not constant")
      fixed = std::to_string(v0.GetPlainDesignValue());
  } else
    fixed = "not constant";

  if(v0.GetLowerBound() == test.GetLowerBound()) {
    if(lower != "not constant")
      lower = std::to_string(v0.GetLowerBound());
  } else
    lower = "not constant";

  if(v0.GetUpperBound() == test.GetUpperBound()) {
    if(upper != "not constant")
      upper = std::to_string(v0.GetUpperBound());
  } else
    upper = "not constant";
}

void SpaghettiDesign::Noodle::ToInfo(PtrParamNode in, Variable::Type type, const std::string& label, int dim) const
{
  std::string fixed;
  std::string lower;
  std::string upper;
  const Variable* v0 = nullptr;
  int size = -1;

  if(type == Variable::NORMAL || type == Variable::RADIUS)
  {
    const StdVector<Variable>& vec = type == Variable::NORMAL ? a : r;
    if(vec.IsEmpty())
      return;

    v0 = &(vec[0]);
    size = vec.GetSize();

    // check if the element values are valid for all (e.g. from mathparser)
    for(unsigned int i = 0; i < vec.GetSize(); i++)
      CompareToInfoHelper(*v0, vec[i], fixed, lower, upper);
  }
  else
  {
    assert(type == Variable::NODE);
    assert(dim >= 0 && dim <= 2);
    if(points.GetSize() <= 2)
      return;

    v0 = &(points[1][dim]);
    assert(v0->tip == INNER);
    size = points.GetSize() - 2; // minus P and Q

    for(unsigned int i = 1; i < points.GetSize()-1; i++)
      CompareToInfoHelper(*v0, points[i][dim], fixed, lower, upper);
  }

  v0->ToInfo(in->Get(label, ParamNode::APPEND));
  in->Get(label + "/size")->SetValue(size);
  if(v0->fixed)
    in->Get(label + "/fixed")->SetValue(fixed);
  else {
    in->Get(label + "/lower")->SetValue(lower);
    in->Get(label + "/upper")->SetValue(upper);
  }
}

std::string SpaghettiDesign::Variable::ToString(const StdVector<Variable>& vec, bool show_fixed)
{
  std::stringstream ss;
  ss << "[";
  for(unsigned int i = 0; i < vec.GetSize(); i++)
  {
    const Variable& v = vec[i];
    if(show_fixed && v.fixed)
      ss << "*";
    ss << v.GetPlainDesignValue();
    if(i < vec.GetSize()-1)
      ss << ", ";
  }
  ss << "]";
  return ss.str();
}

Vector<double> SpaghettiDesign::Variable::AsVector(const StdVector<Variable>& vec)
{
  Vector<double> res(vec.GetSize());
  for(unsigned int i = 0; i < vec.GetSize(); i++)
    res[i] = vec[i].GetPlainDesignValue();
  return res;
}

void SpaghettiDesign::Variable::InitInnerNode(int noodle_idx, Dof dof, double val, double lower, double upper)
{
  this->type_ = NODE;
  this->fixed = false;
  this->noodle = noodle_idx;
  this->dof_ = dof;
  this->tip = INNER;
  SetDesign(val);
  SetLowerBound(lower);
  SetUpperBound(upper);
  LOG_DBG(pasta) << "V:IIN " << ToString();
}

void SpaghettiDesign::Variable::Parse(PtrParamNode pn, int noodle_idx, double interpolate_value)
{
  this->noodle = noodle_idx;

  MathParser* mp = domain->GetMathParser();
  unsigned int handle = mp->GetNewHandle();
  // there are not much additional variables stored.

  string value;

  if(pn->Has("fixed")) {
    value = pn->Get("fixed")->As<string>();
    if(pn->Has("upper") || pn->Has("lower") || pn->Has("initial"))
      throw Exception("When a '" + pn->GetName() + "' is 'fixed', don't use 'initial', 'lower' and 'upper'");
    fixed = true;
  } else {
    if(!pn->Has("upper") || !pn->Has("lower") || !pn->Has("initial"))
      throw Exception("'" + pn->GetName() + "' needs 'initial', 'lower' and 'upper' when not 'fixed'");
    value = pn->Get("initial")->As<string>();
    mp->SetExpr(handle, pn->Get("lower")->As<string>());
    SetLowerBound(mp->Eval(handle));
    mp->SetExpr(handle, pn->Get("upper")->As<string>());
    SetUpperBound(mp->Eval(handle));

    fixed = false;
  }

  // this might contain a formula or simply a value
  if(value == "interpolated")
  {
    assert(interpolate_value != -12.34);
    SetDesign(interpolate_value);
  }
  else
  {
    mp->SetExpr(handle, value);
    SetDesign(mp->Eval(handle));
  }

  mp->ReleaseHandle(handle);

  if(pn->Has("tip")) // node requires tip
    tip = SpaghettiDesign::tip.Parse(pn->Get("tip")->As<string>());

  if(pn->Has("dof")) // node and inner require top
    dof_ = dof.Parse(pn->Get("dof")->As<string>());

  type_ = pn->GetName() == "inner" ? NODE : type.Parse(pn->GetName()); // we use inner only for the problem xml, it is a NODE int the end

  if(type_ != NODE && (tip != NO_TIP || dof_ != NOT_SET))
    throw Exception("Within 'noodle' use 'dof' and 'tip' only for 'node'");
}

void SpaghettiDesign::Variable::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  if(type_ == NODE)
    in->Get("dof")->SetValue(ShapeParamElement::dof.ToString(dof_));
  if(type_ == NODE && tip != NO_TIP) // inner has no top -> one could use "inner" ?!
    in->Get("tip")->SetValue(SpaghettiDesign::tip.ToString(tip));

  // in the NORMAL case the values below might be overwritten by Noode::ToInfo().
  if(fixed)
   in->Get("fixed")->SetValue(GetPlainDesignValue());
  else {
   in->Get("lower")->SetValue(lower_);
   in->Get("upper")->SetValue(upper_);
  }
}


std::string SpaghettiDesign::Variable::ToString() const
{
  std::stringstream ss;
  ss << ShapeParamElement::ToString() << " noodle=" << noodle << " type=" << type.ToString(type_) << " dof=" << dof.ToString(dof_) << " tip=" << tip;
  return ss.str();
}

std::string SpaghettiDesign::Variable::GetLabel() const
{
  std::stringstream ss;
  ss << "s" << noodle << "_"; // 0-based
  switch(type_)
  {
  case NODE:
    if(tip == START)
      ss << "p";
    else if(tip == END)
      ss << "q";
    else // We increment after x/y/z. As P is first, inner is 1-based as desired
      ss << "i" << (unsigned int) GetIndex() / dim_;
    ss << dof.ToString(dof_);
    break;
  case PROFILE:
    ss << "p";
    break;
  case ALPHA:
    ss << "alpha";
    break;
  case NORMAL: // intentionally both cases
  case RADIUS:
  {
    assert(domain->GetOptimization() != NULL);
    SpaghettiDesign* sd = dynamic_cast<SpaghettiDesign*>(domain->GetOptimization()->GetDesign());
    const Noodle& n = sd->spaghetti[noodle];
    ss << (type_ == NORMAL ? "a" : "r");
    assert(!(n.a.GetSize() > 0 && n.r.GetSize() > 0)); // we are in GetLable() and shall have it
    unsigned int pos = GetIndex() - (type_ == NORMAL ? n.a.First().GetIndex() : n.r.First().GetIndex());
    ss << (pos + 1); // 1-based
    break;
  }

  default:
    assert(false);
    break;
  }

  return ss.str();

}


} // end of namespace

