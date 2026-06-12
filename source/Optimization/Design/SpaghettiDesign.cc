#include "Utils/CfsNumpy.hh"

#include "Optimization/Design/SpaghettiDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Utils/ToolsFull.hh"
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

  OpenGradPlot(pn);

  // map
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
  FeaturedDesign::ToInfo(em);

  sp_info_ = info_->Get("features");

  if(radius > 0.0)
    sp_info_->Get("radius")->SetValue(radius);
  sp_info_->Get("boundary/transition")->SetValue(transition);
  sp_info_->Get("orientation")->SetValue(orientation.ToString(orientation_));

  // python stuff written with PythonInit()

  // get keys of python glob.info_field -> to be repeated e.g. in MapFeatureGradient() as we have initially not much info
  StdVector<string> keys = PythonGetInfoFieldKeys();
  sp_info_->Get("info_field/genericElem")->SetValue(keys.ToString()); // TODO use nice Vector::ToString()
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
      StdVector<FeatureVariable>& p = s.points[pi];
      PyObject* lst = PyList_New(p.GetSize());
      for(unsigned int li = 0; li < p.GetSize(); li++)
        PyList_SetItem(lst, li, PyFloat_FromDouble(p[li].GetPlainDesignValue()));
      PyTuple_SetItem(pyp, pi, lst);
    }

    // create a tuple for a or r - is allowed to be empty
    const StdVector<FeatureVariable>& vec = s.type == Noodle::NORMAL ? s.a : s.r;
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

  assert(!map.IsEmpty());
  assert(tmp.GetSize() >= map.GetSize()); // map can be smaller if the domain is not convex

  for(auto& item : map)
    item.elemval->SetDesign(tmp[item.lexicographic_pos]);

  LOG_DBG(pasta)  << "MFTD: old mid=" << mapped_design_ << " current did=" << design_id << " size=" << tmp.GetSize();
  LOG_DBG2(pasta) << "MFTD: rho=" << tmp.ToString();

  mapped_design_ = design_id;

  py_timer->Stop();
  mapping_timer_->Stop();
}

void SpaghettiDesign::MapFeatureGradient(Function* f)
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
  for(auto& item : map)
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




void SpaghettiDesign::SetupVirtualShapeElementMapBending(Feature* f, StdVector<Function::Local::Identifier>& vem, StdVector<BaseDesignElement*>& nodes, bool two_signs, int sign_1, int sign_2)
{
      // with  many segments, the typical case would be ai-1,ai,ai+1 (N N N)
      // but we have the cases for only two segments (0 N 0)
      // or for start ( 0 N N) or end ( N N 0)

      Noodle* s = dynamic_cast<Noodle*>(f);
      StdVector<FeatureVariable>& P = f->points.First();

      for(unsigned ai = 0; ai < s->a.GetSize(); ai++)
      {
        StdVector<BaseDesignElement*> buddies(nodes);
        assert(buddies.GetSize() == nodes.GetSize());
        assert(buddies[0] == nodes[0]); // no empty nodes for bending

        if(s->a.GetSize() == 1)
        {
          buddies.Push_back(&s->a[0]);
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
            buddies.Push_back(&s->a[0]);
            buddies.Push_back(&s->a[1]);
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::ZNN;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::ZNN;
            }
          }
          else if(ai == s->a.GetSize()-1)
          {
            buddies.Push_back(&s->a[ai-1]);
            buddies.Push_back(&s->a[ai]);
            vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::NNZ;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&P[0], buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::NNZ;
            }
          }
          else
          {
            assert(ai > 0 && ai < s->a.GetSize()-1);
            buddies.Push_back(&s->a[ai-1]);
            buddies.Push_back(&s->a[ai]);
            buddies.Push_back(&s->a[ai+1]);
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
        LOG_DBG(pasta) << "SVSEM: f=" << f->ToString() << " s=" << s->idx << " id=" << vem.Last().ToString();
      } // end ai
    }


int SpaghettiDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  // practically updates the FeatureVariable objects within spaghetti. Set's design_id
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
    FeatureVariable*    var = dynamic_cast<FeatureVariable*>(shape_param_[i]);
    PtrParamNode pn  = list[i];

    // do a lot of validation
    const std::string pnnr = pn->Get("nr")->As<string>();
    if(pn->Get("nr")->As<unsigned int>() != i)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'nr' value " << i);

    if(pn->Get("shape")->As<int>() != var->feature)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'shape' value " << var->feature);

    if(DesignElement::type.Parse(pn->Get("type")->As<string>()) != var->GetType())
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'type' value " << DesignElement::type.ToString(var->GetType()));

    if(var->GetType() == FeatureVariable::NODE)
    {
      if(ShapeParamElement::dof.Parse(pn->Get("dof")->As<string>()) != var->dof_)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'dof' value " << ShapeParamElement::dof.ToString(var->dof_));
      if(FeatureVariable::tip_enum.Parse(pn->Get("tip")->As<string>()) != var->tip)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'tip' value " << FeatureVariable::tip_enum.ToString(var->tip));
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

void SpaghettiDesign::SetupParsedFeatures(PtrParamNode base)
{
  StdVector<PtrParamNode> pnl = base->GetList("noodle");
  spaghetti.Resize(pnl.GetSize());
  features_.Resize(pnl.GetSize());
  for(unsigned int i = 0; i < pnl.GetSize(); i++)
  {
    spaghetti[i].Parse(pnl[i], i);
    features_[i] = &  spaghetti[i];
  }
}

void SpaghettiDesign::SetupDesign(PtrParamNode base)
{
  FeaturedDesign::SetupDesign(base); // makes use of SetupParsedFeatures()
   
  // the order is by noodles: nodes, profile, [alpha|, normals/radii
  // reading shall be transparent due to shape (noodle) idx.
  for(Noodle& noodle : spaghetti)
  {
    for(auto& p : noodle.points) // P, [inner], Q
      for(auto& v : p)
        AddVariable(&v);

    assert(noodle.a.IsEmpty() || noodle.r.IsEmpty()); // or both are empty
    for(FeatureVariable& var : noodle.r)
      AddVariable(&var);

    AddVariable(&(noodle.p));

    if(noodle.HasAlpha())
      AddVariable(&(noodle.alpha));

    for(FeatureVariable& var : noodle.a)
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
  int n = noodle->Get("segments")->As<int>();
  if(n < 1)
    throw Exception("'noodle' requires at least 1 segment");

  if(n == 1 && noodle->Has("inner"))
    throw Exception("'noodle' has defined 1 segment but 'inner' is given");

  assert(points.IsEmpty());  
  Feature::Parse(noodle, idx);
  assert(points.GetSize() == 2);  
  if(dim_ == 3)
  {
    // we need to add inner points and move the last point to the end
    points.Resize(n + 1);
    points.Last() = points[1]; 

    StdVector<FeatureVariable>& start = points.First();
    StdVector<FeatureVariable>& end   = points.Last();

    // add inner points for 3d noodle   
    Vector<double> P = FeatureVariable::AsVector(start);
    Vector<double> Q = FeatureVariable::AsVector(end);

    for(int pi = 1; pi < n; pi++) // n = 1 -> no inner point; n = 2 -> 1 inner point 1; n = 3 -> 2 inner points 1,2
    {
      StdVector<FeatureVariable>& point = points[pi];
      point.Resize(dim_);
    
      for(unsigned int d = 0; d < dim_; d++)
      {
        double inter = Interpolate(P[d],Q[d],pi,points.GetSize());
        string dof = ShapeParamElement::dof.ToString((ShapeParamElement::Dof) d);
        // either we have inner or require clear lower and upper from P and Q
        if(noodle->HasByVal("inner", "dof", dof)) // either interpolate for value or math parser stuff
        {
          // not really sure, if NODE is right here - but is has no use anyway here
          point[d].Parse(noodle->GetByVal("inner", "dof", dof), this->idx,  DesignElement::NODE, inter);
          point[d].tip = FeatureVariable::INNER;
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
      a[i].Parse(noodle->Get("normal"),this->idx, DesignElement::NORMAL);

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
      r[i].Parse(noodle->Get("radius"),this->idx, DesignElement::RADIUS);

    opt_variables_ += !r.IsEmpty() && !r[0].fixed ? r.GetSize() : 0;

    alpha.Parse(noodle->Get("alpha"),idx, DesignElement::ALPHA); // exactly one
    opt_variables_ += alpha.fixed ? 0 : 1;
  }

  LOG_DBG(pasta) << "N:start ov=" << opt_variables_ << " " << ToString();
}

int SpaghettiDesign::Noodle::GetTotalVariables() const
{
  int vars = Feature::GetTotalVariables();

  assert(a.IsEmpty() || r.IsEmpty());
  vars += a.GetSize() + r.GetSize();
  if(HasAlpha())
    vars += 1;
  return vars;
}


string SpaghettiDesign::Noodle::ToString() const
{
  std::stringstream ss(Feature::ToString());
  ss << " type=" << (type == NORMAL ? "normal" : "points");
  ss << " a=" << FeatureVariable::ToString(a);
  ss << " r=" << FeatureVariable::ToString(r);
  return ss.str();
}

void SpaghettiDesign::Noodle::ToInfo(PtrParamNode in)
{
  Feature::ToInfo(in);

  in->Get("segments")->SetValue(a.GetSize() + 1);

  if(HasInner())
    for(unsigned int d = 0; d < dim_; d++)
      ToInfo(in, FeatureVariable::NODE, "i" + Dof(d), d);

  if(HasAlpha())
    alpha.ToInfo(in->Get("alpha", ParamNode::APPEND));

  // will do nothing for empty vectors
  ToInfo(in, FeatureVariable::NORMAL, "a");
  ToInfo(in, FeatureVariable::RADIUS, "r");
}

void SpaghettiDesign::Noodle::ToInfo(PtrParamNode in, FeatureVariable::Type type, const std::string& label, int dim) const
{
  if(type != FeatureVariable::NORMAL && type != FeatureVariable::RADIUS)
    Feature::ToInfo(in, type, label, dim);
  else
  {
    std::string fixed;
    std::string lower;
    std::string upper;
    const FeatureVariable* v0 = nullptr;
    int size = -1;

    const StdVector<FeatureVariable>& vec = type == FeatureVariable::NORMAL ? a : r;
    if(vec.IsEmpty())
      return;

    v0 = &(vec[0]);
    size = vec.GetSize();

    // check if the element values are valid for all (e.g. from mathparser)
    for(unsigned int i = 0; i < vec.GetSize(); i++)
      FeatureVariable::CompareToInfoHelper(v0, &vec[i], fixed, lower, upper);

    v0->ToInfo(in->Get(label, ParamNode::APPEND));
    in->Get(label + "/size")->SetValue(size);
    if(v0->fixed)
      in->Get(label + "/fixed")->SetValue(fixed);
    else {
      in->Get(label + "/lower")->SetValue(lower);
      in->Get(label + "/upper")->SetValue(upper);
    }
  }
}

} // end of namespace

