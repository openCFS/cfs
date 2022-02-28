#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include "Optimization/Design/SpaghettiDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Utils/PythonKernel.hh"
#include "Optimization/Function.hh"
#include "Optimization/TransferFunction.hh"

using std::string;

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

  PtrParamNode pn = empn->Get(ErsatzMaterial::method.ToString(method)); // spaghetti or spaghettiParamMat

  combine_ = combine.Parse(pn->Get("combine")->As<string>());
  boundary_ = boundary.Parse(pn->Get("boundary")->As<string>());

  transition = pn->Get("transition")->As<double>(); // make optional
  radius = pn->Get("radius")->As<double>();

  if(pn->Get("gradplot")->As<bool>()) // todo: move to FeaturedDesign
   gradplot_.open((progOpts->GetSimName() + ".grad.dat").c_str()); // the auto destructor does the job.

  // n_ and nx_, ny_, nz_
  SetupMeshStructure();

  // map_
  SetupMapping();

  // setup shaphetti with noodles and shape_param_ and opt_shape_param_
  SetupDesign(pn);

  // get rhomin from first density element and check transfer function
  assert(!data.IsEmpty());
  const DesignElement& de = data.First();
  rhomin = de.GetLowerBound(); // at least in python use for all

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
  sp_info_->Get("radius")->SetValue(radius);
  sp_info_->Get("boundary/type")->SetValue(boundary.ToString(boundary_));
  sp_info_->Get("boundary/transition")->SetValue(transition);

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

  // def cfs_init(rhomin, radius, boundary, transition, combine, nx, ny, nz):
  PyObject* func = PyObject_GetAttrString(module_, "cfs_init");
  PythonKernel::CheckPythonFunction(func, "cfs_init");

  PyObject* arg = PyTuple_New(11);
  PyTuple_SetItem(arg, 0, PyFloat_FromDouble(rhomin));
  PyTuple_SetItem(arg, 1, PyFloat_FromDouble(radius));
  PyTuple_SetItem(arg, 2, PyUnicode_FromString(boundary.ToString(boundary_).c_str()));
  PyTuple_SetItem(arg, 3, PyFloat_FromDouble(transition));
  PyTuple_SetItem(arg, 4, PyUnicode_FromString(combine.ToString(combine_).c_str()));
  PyTuple_SetItem(arg, 5, PyLong_FromLong(nx_));
  PyTuple_SetItem(arg, 6, PyLong_FromLong(ny_));
  PyTuple_SetItem(arg, 7, PyLong_FromLong(nz_));
  PyTuple_SetItem(arg, 8, PyFloat_FromDouble(dx_));

  PyObject* des = PyTuple_New(design.GetSize());
  for(unsigned int i = 0; i < design.GetSize(); i++)
    PyTuple_SetItem(des, i,  PyUnicode_FromString(DesignElement::type.ToString(design[i].design).c_str()));
  PyTuple_SetItem(arg, 9, des); // steals the reference, so no need to decref

  PyTuple_SetItem(arg, 10, PythonKernel::CreatePythonDict(pyopts));
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
    // create a tuple for a, is allowed to be empty
    PyObject* pya = PyTuple_New(s.a.GetSize());
    for(unsigned int ai = 0; ai < s.a.GetSize(); ai++)
      PyTuple_SetItem(pya, ai, PyFloat_FromDouble(s.a[ai]));

    PyObject* arg = PyTuple_New(7);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(s.idx));
    PyTuple_SetItem(arg, 1, PyFloat_FromDouble(s.P[0]));
    PyTuple_SetItem(arg, 2, PyFloat_FromDouble(s.P[1]));
    PyTuple_SetItem(arg, 3, PyFloat_FromDouble(s.Q[0]));
    PyTuple_SetItem(arg, 4, PyFloat_FromDouble(s.Q[1]));
    PyTuple_SetItem(arg, 5, pya);
    PyTuple_SetItem(arg, 6, PyFloat_FromDouble(s.p));

    LOG_DBG2(pasta) << "PUS: s=" << s.idx << " a=" << s.a.GetSize() << " func=" << func << " arg=" << arg << " pya=" << pya;

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
  assert(tmp.GetSize() == map_.GetSize());

  for(unsigned int i = 0; i < tmp.GetSize(); i++)
    map_[i].elemval->SetDesign(tmp[i]);

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
  Vector<double> drho(data.GetSize());
  for(const DesignElement& de : data)
    drho[de.GetIndex()] = de.GetPlainGradient(f);

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

  vem.Reserve(spaghetti.GetSize()); // assume nothing fixec

  assert(dim_ == 2);
  StdVector<BaseDesignElement*> nodes;

  bool two_signs = locality == Function::Local::FUNCTION_SPECIFIC_TWO_SIGNS;

  int sign_1 = two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = -1;

  for(Noodle& s : spaghetti)
  {
    nodes.Clear();
    // assume nothing fixed
    if(s.px.fixed || s.py.fixed || s.qx.fixed || s.qy.fixed)
      throw Exception("distance constraints currently only for non-fixed nodes");

    // px is element, then py, then qx then qy
    nodes.Push_back(&s.py);
    nodes.Push_back(&s.qx);
    nodes.Push_back(&s.qy);

    if(f->GetType() == Function::DISTANCE)
    {
      if(f->GetDesignType() != BaseDesignElement::NODE)
        throw Exception("Configuration error: distance needs to be defined for node.");

      vem.Push_back(Function::Local::Identifier(&s.px, nodes));
      LOG_DBG(pasta) << "SVSEM: f=" << f->ToString() << " s=" << s.idx << " id=" << vem.Last().ToString();
    }
    else // Bending
    {
      // with  many segments, the typical case would be ai-1,ai,ai+1 (N N N)
      // but we have the cases for only two segments (0 N 0)
      // or for start ( 0 N N) or end ( N N 0)

      for(unsigned ai = 0; ai < s.a_var.GetSize(); ai++)
      {
        StdVector<BaseDesignElement*> buddies(nodes);
        assert(buddies.GetSize() == nodes.GetSize());
        assert(buddies[0] == nodes[0]); // no empty nodes for bending

        if(s.a_var.GetSize() == 1)
        {
          buddies.Push_back(&s.a_var[0]);
          vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_1));
          vem.Last().bending = Function::Local::Identifier::ZNZ;
          if(two_signs){
            vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_2));
            vem.Last().bending = Function::Local::Identifier::ZNZ;
          }
        }
        else
        {
          if(ai == 0)
          {
            buddies.Push_back(&s.a_var[0]);
            buddies.Push_back(&s.a_var[1]);
            vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::ZNN;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::ZNN;
            }
          }
          else if(ai == s.a_var.GetSize()-1)
          {
            buddies.Push_back(&s.a_var[ai-1]);
            buddies.Push_back(&s.a_var[ai]);
            vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::NNZ;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::NNZ;
            }
          }
          else
          {
            assert(ai > 0 && ai < s.a_var.GetSize()-1);
            buddies.Push_back(&s.a_var[ai-1]);
            buddies.Push_back(&s.a_var[ai]);
            buddies.Push_back(&s.a_var[ai+1]);
           // vem.Push_back(Function::Local::Identifier(&s.px, buddies));
            //vem.Last().bending = Function::Local::Identifier::NNN;
            vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_1));
            vem.Last().bending = Function::Local::Identifier::NNN;
            if(two_signs){
              vem.Push_back(Function::Local::Identifier(&s.px, buddies, sign_2));
              vem.Last().bending = Function::Local::Identifier::NNN;
            }
          }
        }
        LOG_DBG(pasta) << "SVSEM: f=" << f->ToString() << " s=" << s.idx << " id=" << vem.Last().ToString();
      } // end ai
    } // end bending
  } // end noodle
}


int SpaghettiDesign::ReadDesignFromExtern(const double* space_in)
{
  // practically updates the Variable objects within spaghetti. Set's design_id
  FeaturedDesign::ReadDesignFromExtern(space_in);

  // update also the convenience variables
  for(Noodle& s : spaghetti)
    s.Update();

  PythonUpdateSpaghetti();

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

  for(Noodle& s : spaghetti)
    s.Update();

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
    assert(dim_ == 2);
    total += noodle.GetTotalVariables();
    opt   += noodle.GetOptVariables();
  }
  // spaghetti contain the variable instances, now flatten the design

  shape_param_.Reserve(total);
  opt_shape_param_.Reserve(opt);

  // the order of 1) nodes, 2) profiles, 3) normals
  // could also be by noodle - this is just a decision. Change it, when it hurts.
  // reading shall be transparent die to shape (noodle) idx.
  for(Noodle& noodle : spaghetti)
  {
    AddVariable(&(noodle.px));
    AddVariable(&(noodle.py));
    AddVariable(&(noodle.qx));
    AddVariable(&(noodle.qy));
    AddVariable(&(noodle.p_var));
    for(Variable& var : noodle.a_var)
      AddVariable(&var);
  }


  // we did exact space "estimation"
  assert(shape_param_.GetCapacity() == shape_param_.GetSize());
  assert(opt_shape_param_.GetCapacity() == opt_shape_param_.GetSize());
}


void SpaghettiDesign::Noodle::Parse(PtrParamNode noodle, int idx)
{
  this->idx = idx;

  int n = noodle->Get("segments")->As<int>();

  assert(dim_ == 2);
  px.Parse(noodle->GetByVal("node", "dof", "x", "tip", "start"),idx);
  py.Parse(noodle->GetByVal("node", "dof", "y", "tip", "start"),idx);
  qx.Parse(noodle->GetByVal("node", "dof", "x", "tip", "end"),idx);
  qy.Parse(noodle->GetByVal("node", "dof", "y", "tip", "end"),idx);
  p_var.Parse(noodle->Get("profile"),idx); // exactly one

  // number of normals is segments-1
  if(n < 1)
    EXCEPTION("'noodle' requires at least 1 segment");
  if(n == 1 && noodle->Has("normal"))
    EXCEPTION("'noodle' has defined 1 segment but 'normal' is given");
  if(n >= 2 && ! noodle->Has("normal"))
    EXCEPTION("'noodle' has defined " << n << " segments but no 'normal' given");

  a_var.Resize(n-1);
  for(int i = 0; i < n-1; i++)
    a_var[i].Parse(noodle->Get("normal"),this->idx); // todo handle mathparser

  opt_variables_  = px.fixed ? 0 : 1;
  opt_variables_ += py.fixed ? 0 : 1;
  opt_variables_ += qx.fixed ? 0 : 1;
  opt_variables_ += qy.fixed ? 0 : 1;
  opt_variables_ += p_var.fixed  ? 0 : 1;
  opt_variables_ += !a_var.IsEmpty() && !a_var[0].fixed ? a_var.GetSize() : 0;

  LOG_DBG(pasta) << "N:P px=" << px.GetPlainDesignValue() << " py=" << py.GetPlainDesignValue();

  Update();
}

void SpaghettiDesign::Noodle::Update()
{
  P.Resize(2); // cheap so save doing it in constructor
  P[0] = px.GetPlainDesignValue();
  P[1] = py.GetPlainDesignValue();

  Q.Resize(2);
  Q[0] = qx.GetPlainDesignValue();
  Q[1] = qy.GetPlainDesignValue();

  p = p_var.GetPlainDesignValue();

  a.Resize(a_var.GetSize());
  for(unsigned int i = 0; i < a_var.GetSize(); i++)
    a[i] = a_var[i].GetPlainDesignValue();

  LOG_DBG(pasta) << "N:U " << ToString();
}

string SpaghettiDesign::Noodle::ToString()
{
  std::stringstream ss;
  ss << "idx=" << idx << " P=" << P.ToString() << " Q=" << Q.ToString() << " p=" << p << " a=" << a.ToString();
  return ss.str();
}

void SpaghettiDesign::Noodle::ToInfo(PtrParamNode in)
{
  in->Get("id")->SetValue(idx);
  in->Get("segments")->SetValue(a.GetSize() + 1);

  // one could combine P and Q as coordinates.
  px.ToInfo(in->Get("px", ParamNode::APPEND));
  py.ToInfo(in->Get("py", ParamNode::APPEND));
  qx.ToInfo(in->Get("qx", ParamNode::APPEND));
  qy.ToInfo(in->Get("qy", ParamNode::APPEND));
  p_var.ToInfo(in->Get("p", ParamNode::APPEND));
  if(!a_var.IsEmpty()) {
    const Variable& a0 = a_var.First();
    a0.ToInfo(in->Get("a", ParamNode::APPEND));
    in->Get("a/size")->SetValue(a_var.GetSize());
    // check if the first element values are valid for all (e.g. from mathparser)
    for(unsigned int i = 1; i < a_var.GetSize(); i++) {
      const Variable& t = a_var[i];
      assert(a0.fixed == t.fixed);
      if(a0.fixed) {
        if(a0.GetPlainDesignValue() != t.GetPlainDesignValue()) {
          in->Get("a/fixed")->SetValue("not constant");
          break;
        }
      } else {
        bool cancel = false;
        if(a0.GetLowerBound() != t.GetLowerBound()) {
          in->Get("a/lower")->SetValue("not constant");
          cancel = true;
        }
        if(a0.GetUpperBound() != t.GetUpperBound()) {
          in->Get("a/upper")->SetValue("not constant");
          cancel = true;
        }
        if(cancel)
          break;
      }
    }
  }

}



void SpaghettiDesign::Variable::Parse(PtrParamNode pn, int noodle_idx)
{
  this->noodle = noodle_idx;

  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();
  // there are not much additional variables stored.

  string value;

  // todo: add mathparser!!
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
  mp->SetExpr(handle, value);
  SetDesign(mp->Eval(handle));

  mp->ReleaseHandle(handle);

  if(pn->Has("tip"))
    tip = SpaghettiDesign::tip.Parse(pn->Get("tip")->As<string>());

  if(pn->Has("dof"))
    dof_ = dof.Parse(pn->Get("dof")->As<string>());

  type_ = type.Parse(pn->GetName());

  if(type_ == NODE && (tip == NO_TIP || dof_ == NOT_SET))
    throw Exception("A 'node' within 'noodle' requires 'dof' and 'tip'");
  if(type_ != NODE && (tip != NO_TIP || dof_ != NOT_SET))
    throw Exception("Within 'noodle' use 'dof' and 'tip' only for 'node'");
}

void SpaghettiDesign::Variable::ToInfo(PtrParamNode in) const
{
  in->Get("type")->SetValue(type.ToString(type_));
  if(type_ == NODE) {
    in->Get("dof")->SetValue(ShapeParamElement::dof.ToString(dof_));
    in->Get("tip")->SetValue(SpaghettiDesign::tip.ToString(tip));
  }
  // in the NORMAL case the values below might be overwritten by Noode::ToInfo().
  if(fixed)
   in->Get("fixed")->SetValue(GetPlainDesignValue());
  else {
   in->Get("lower")->SetValue(lower_);
   in->Get("upper")->SetValue(upper_);
  }
}


std::string SpaghettiDesign::Variable::GetLabel() const
{
  std::stringstream ss;
  ss << "s" << noodle << "_"; // 0-based
  switch(type_)
  {
  case NODE:
    ss << (tip == START ? "p" : "q");
    ss << dof.ToString(dof_);
    break;
  case PROFILE:
    ss << "p";
    break;
  case NORMAL: {
    ss << "a";
    assert(domain->GetOptimization() != NULL);
    SpaghettiDesign* sd = dynamic_cast<SpaghettiDesign*>(domain->GetOptimization()->GetDesign());
    const Noodle& n = sd->spaghetti[noodle];
    assert(!n.a_var.IsEmpty()); // we are in GetLable() and shall have it
    unsigned int pos = GetIndex() - n.a_var.First().GetIndex();
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

