#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include "Optimization/Design/SpaghettiDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/PythonTools.hh"
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

  PtrParamNode pn = empn->Get("spaghetti");

  combine_ = combine.Parse(pn->Get("combine")->As<string>());
  boundary_ = boundary.Parse(pn->Get("boundary")->As<string>());

  transition = pn->Get("transition")->As<double>(); // make optional
  radius = pn->Get("radius")->As<double>();

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
  const TransferFunction* tf = GetTransferFunction(&de);
  if(tf->GetType() != TransferFunction::IDENTITY)
    throw Exception("Up to now Spaghetti expects identity transfer function"); // extend!

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
   if(python)
     PythonDestructor();
}


void SpaghettiDesign::ToInfo(ErsatzMaterial* em)
{
  AuxDesign::ToInfo(em);

  PtrParamNode sm = info_->Get("spaghetti");

  sm->Get("combine")->SetValue(combine.ToString(combine_));
  sm->Get("radius")->SetValue(radius);
  sm->Get("boundary/type")->SetValue(boundary.ToString(boundary_));
  sm->Get("boundary/transition")->SetValue(transition);

  PtrParamNode py = sm->Get("python");
  if(python) {
    py->Get("value")->SetValue(pythonfile);
    for(auto& opt : pyopts) {
      PtrParamNode o = py->Get("option", ParamNode::APPEND);
      o->Get("key")->SetValue(opt.first);
      o->Get("value")->SetValue(opt.second);
    }
  } else {
    py->SetValue(false);
  }

  PtrParamNode msh = sm->Get("mesh");
  msh->Get("n")->SetValue(n_.ToString());
  PtrParamNode base = info_->Get("designVariables");
  for(Noodle& s : spaghetti)
    s.ToInfo(base->Get("noodle", ParamNode::APPEND));
}

void SpaghettiDesign::PythonDestructor()
{
  py_timer->Start();

  Py_XDECREF(python);
  python = NULL;
  Py_Finalize();

  py_timer->Stop();
}


void SpaghettiDesign::PythonInit(PtrParamNode pn)
{
  py_timer->Start();

  pythonfile = pn->Get("file")->As<string>();

  pyopts = ParseOptions(pn->GetList("option"));

  string version;
  python = InitializePythonModule(pythonfile, nullptr, &version);

  // def cfs_init(rhomin, radius, boundary, transition, combine, nx, ny, nz):
  PyObject* func = PyObject_GetAttrString(python, "cfs_init");
  CheckPythonFunction(func, "cfs_init");

  PyObject* arg = PyTuple_New(9);
  PyTuple_SetItem(arg, 0, PyFloat_FromDouble(rhomin));
  PyTuple_SetItem(arg, 1, PyFloat_FromDouble(radius));
  PyTuple_SetItem(arg, 2, PyUnicode_FromString(boundary.ToString(boundary_).c_str()));
  PyTuple_SetItem(arg, 3, PyFloat_FromDouble(transition));
  PyTuple_SetItem(arg, 4, PyUnicode_FromString(combine.ToString(combine_).c_str()));
  PyTuple_SetItem(arg, 5, PyLong_FromLong(nx_));
  PyTuple_SetItem(arg, 6, PyLong_FromLong(ny_));
  PyTuple_SetItem(arg, 7, PyLong_FromLong(nz_));
  PyTuple_SetItem(arg, 8, CreatePythonDict(pyopts));
  PyObject* ret = PyObject_CallObject(func, arg);
  CheckPythonReturn(ret);

  Py_XDECREF(ret);
  Py_XDECREF(arg);
  Py_XDECREF(func);

  py_timer->Stop();
}

void SpaghettiDesign::PythonUpdateSpaghetti()
{
  py_timer->Start();

  // def cfs_set_spaghetti(id, px, py, qx, qy, a_list, p):
  PyObject* func = PyObject_GetAttrString(python, "cfs_set_spaghetti");
  CheckPythonFunction(func, "cfs_set_spaghetti");

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
    CheckPythonReturn(ret);

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
  assert(python);

  LOG_DBG(pasta) << "MFTD called";

  if(design_id != mapped_design_)
    PythonUpdateSpaghetti(); // has it's own timer, it would top our

  py_timer->Start();

  // up to now we only can do python
  PyObject* func = PyObject_GetAttrString(python, "cfs_map_to_density");
  CheckPythonFunction(func, "cfs_map_to_density");

  PyObject* ret = PyObject_CallObject(func, NULL);
  CheckPythonReturn(ret);

  // more convenient to use a temporary object than to iterate directly on the numpy object
  Vector<double> tmp(ret,true); // decref ret
  Py_XDECREF(func);

  assert(!map_.IsEmpty());
  assert(tmp.GetSize() == map_.GetSize());

  for(unsigned int i = 0; i < tmp.GetSize(); i++)
    map_[i].rho->SetDesign(tmp[i]);

  LOG_DBG(pasta)  << "MFTD: old mid=" << mapped_design_ << " current did=" << design_id << " size=" << tmp.GetSize();
  LOG_DBG2(pasta) << "MFTD: rho=" << tmp.ToString();

  mapped_design_ = design_id;

  py_timer->Stop();
  mapping_timer_->Stop();
}

void SpaghettiDesign::MapFeatureGradient(const Function* f)
{
  gradient_timer_->Start();

  assert(python);

  assert(mapped_design_ == design_id);

  py_timer->Start();

  // python needs to be smart by itself to speedup gradient calculation because we will call this
  // for several functions for the same design set (compliance, volume, ...)

  // we cannot create a numpy array in C - so get it from python to be filled with the rho-density to be given to python
  PyObject* afunc = PyObject_GetAttrString(python, "cfs_get_drho_vector");
  CheckPythonFunction(afunc, "cfs_get_drho_vector");

  PyObject* np_rho = PyObject_CallObject(afunc, NULL);
  CheckPythonReturn(np_rho);

  // we create a temporary vector for the gradient values such that we can use Export() - just for convenience.
  // With Python we expect no hpc performance anyway
  Vector<double> drho(data.GetSize());
  for(const DesignElement& de : data)
    drho[de.GetIndex()] = de.GetPlainGradient(f);

  // copy data to numpy array - does extensive range checks with good error message
  drho.Export(np_rho);

  // now get the gradient form python. We get all gradients and ignore efficiency optimization by skipping fixed variables
  PyObject* gfunc = PyObject_GetAttrString(python, "cfs_get_gradient");
  CheckPythonFunction(gfunc, "cfs_get_gradient");

  PyObject* arg = PyTuple_New(1);
  Py_INCREF(np_rho); // this is important, as PyTuple_SetItem() steals the reference: https://docs.python.org/3/c-api/intro.html
  PyTuple_SetItem(arg, 0, np_rho);

  PyObject* np_grad = PyObject_CallObject(gfunc, arg);
  CheckPythonReturn(np_grad);

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
  }

  for(Noodle& noodle : spaghetti)
    AddVariable(&(noodle.p_var));

  for(Noodle& noodle : spaghetti)
    for(Variable& var : noodle.a_var)
      AddVariable(&var);

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

} // end of namespace

