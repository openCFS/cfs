#include "Optimization/Context.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/OptimizationMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/BucklingDriver.hh"
#include "Driver/EigenValueDriver.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/MultiSequenceDriver.hh"
#include "Driver/Assemble.hh"
#include "Driver/FormsContexts.hh"
#include "Forms/LinForms/LinearForm.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "PDE/LatticeBoltzmannPDE.hh"

using namespace CoupledField;

DEFINE_LOG(context, "context")

Context::Context()
{
  // to be replaces by any current Excitation::Apply()
  excitation_ = NULL;

  driver = NULL;
  mat = NULL;
  manager_ = NULL;
  me_ = NULL;
  basic_excitations_ = 0;
  sequence = -1;
  active = false;
  context_idx = -1;
  driver_steps_ = 0;
  analysis = BasePDE::NO_ANALYSIS;

  num_harm_freq = 0;
  num_bloch_wave_vectors = 0;
  num_eigenmodes = 0;


  complex_ = false;
  harmonic_ = false;
  eigenvalue_ = false;
  bloch_ = false;
  buckling_ = false;
  homogenization = false; // to be set by the functions
  pde = NULL;
  stt = NO_TENSOR;

}

Context::~Context()
{
  delete mat;
}

void Context::Setup(ContextManager* manager, BasePDE::AnalysisType analyis, PtrParamNode node, int sequence_step)
{
  this->manager_ = manager;
  assert(sequence_step >= 1);
  this->sequence = sequence_step;
  this->context_idx = sequence_step - 1;
  this->analysis = analyis;

  switch(analysis)
  {
  case BasePDE::HARMONIC:
    complex_ = true;
    harmonic_ = true;
    num_harm_freq = HarmonicDriver::GetNumFreq(node);
    break;

  case BasePDE::EIGENFREQUENCY:
    complex_ = true;
    eigenvalue_ = true;
    bloch_                 = EigenFrequencyDriver::DoBloch(node);
    num_bloch_wave_vectors = EigenFrequencyDriver::GetNumBlochWave(node); // 0 when not bloch
    num_eigenmodes         = EigenFrequencyDriver::GetNumModes(node); // also in the non-bloch case
    break;

  case BasePDE::BUCKLING:
    complex_ = true;
    eigenvalue_ = true;
    buckling_ = true;
    // we need num_eigenmodes to reserve enough space in StateContainer::Get
    // however "numModes" is only given for Arpack, so we choose a hopefully large enough number
    if (node->Has("numModes"))
      num_eigenmodes = EigenFrequencyDriver::GetNumModes(node);
    else
      num_eigenmodes = 20;
    break;

  case BasePDE::EIGENVALUE:
    complex_ = true;
    eigenvalue_ = true;
    break;

  case BasePDE::STATIC:
    // FIXME num_static_loads not that easy to determine w/o pde. allows 0 (explixit excitations or homogenization) or more
    break;

  case BasePDE::TRANSIENT:
    assert(false); // not yet implemented
    break;

  case BasePDE::INVERSESOURCE:
    EXCEPTION("Optimizaion with analysis type INVERSESOURCE makes no sense!");
    break;

  case BasePDE::MULTI_SEQUENCE:
  case BasePDE::MULTIHARMONIC:
  case BasePDE::HARMONIC25D:
  case BasePDE::NO_ANALYSIS:
    assert(false);
    break;
  }

  // this is too early to access domain->GetOptimization()->infoNode
}

void Context::Update()
{
  driver = domain->GetSingleDriver();
  SetPDEs();

  ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization()); // not set in first place

  // once is enough but do it only when the pdes are initialized (not yet done in Domain::PostInit()) for the multimaterial case)
  if(pde != NULL && mat == NULL && em != NULL) // might be
  {
    assert(Optimization::context == this); // we are already set as active -> necessary for OptimizationMaterial
    if(em->pn->Has("materials")) {
      ParamNodeList list = em->pn->Get("materials")->GetListByVal("material", "sequence", this->sequence);
      if(list.IsEmpty())
        EXCEPTION("'materials' is given but does not contain 'material' for sequence " << sequence);
      assert(list.GetSize() == 1);
      mat = OptimizationMaterial::CreateInstance(OptimizationMaterial::system.Parse(list[0]->Get("type")->As<std::string>()), em, this);
    }
    else
      mat = OptimizationMaterial::CreateInstance(OptimizationMaterial::system.Parse(em->pn->Get("material")->As<std::string>()), em, this);

    // would be too early in Setup()
    if(domain->GetOptimization() != NULL)
      infoNode = domain->GetOptimization()->optInfoNode->Get(ParamNode::HEADER)->GetByVal("sequence", "step", sequence);
    else
      infoNode = domain->GetInfoRoot()->Get("coreOptContext")->GetByVal("sequence", "step", sequence);

    infoNode->Get("complex")->SetValue(complex_);
    infoNode->Get("harmonic")->SetValue(harmonic_);
    infoNode->Get("eigenvalue")->SetValue(eigenvalue_);
    infoNode->Get("bloch")->SetValue(bloch_);
    infoNode->Get("material")->SetValue(OptimizationMaterial::system.ToString(mat->GetSystem()));

    // read paramMat parameters but not for anisotropic feature mapping which does not use additional parameters
    // PARAM_MAT, SHAPE_PARAM_MAT (not tested), SPAGHETTI_PARAM_MAT
    if(dm == NULL && em->IsParamMat(em->GetMethod()) && em->GetMethod() != ErsatzMaterial::FEATURE_MAPPING_PARAM_MAT) 
    {
      // either em->pn->paramMat/designMaterials/designMaterial or em->pn->spaghetti/designMaterial
      std::string base = em->GetMethod() == ErsatzMaterial::SPAGHETTI_PARAM_MAT ? "spaghettiParamMat" : "paramMat/designMaterials";
      assert(em->pn->Has(base));
      assert(sequence != -1);
      ParamNodeList list = em->pn->Get(base)->GetListByVal("designMaterial", "sequence", sequence);

      if(list.IsEmpty() && em->GetMethod() != ErsatzMaterial::SPAGHETTI_PARAM_MAT)
        EXCEPTION("No designMaterial given in paramMat for sequence step " + std::to_string(sequence));

      if(!list.IsEmpty())
      {
        assert(list.GetSize() == 1);
        assert(mat != NULL);
        domain->GetDesign()->SetDesignMaterial(list[0], mat->GetSystem()); // the single material for the specified sequence
        assert(dm != NULL);
        LOG_DBG3(context) << "CTXT: set design material '" << OptimizationMaterial::system.ToString(mat->GetSystem()) << "' for sequence=" << sequence;
      }
    }
  }
}

EigenFrequencyDriver* Context::GetEigenFrequencyDriver()
{
  return dynamic_cast<EigenFrequencyDriver*>(driver);
}

BucklingDriver* Context::GetBucklingDriver()
{
  return dynamic_cast<BucklingDriver*>(driver);
}

EigenValueDriver* Context::GetEigenValueDriver()
{
  return dynamic_cast<EigenValueDriver*>(driver);
}

HarmonicDriver* Context::GetHarmonicDriver()
{
  return dynamic_cast<HarmonicDriver*>(driver);
}

LatticeBoltzmannPDE* Context::GetLatticeBoltzmannPDE()
{
  return dynamic_cast<LatticeBoltzmannPDE*>(pde);
}

bool Context::DoMultiSequence() const
{
  return manager_->context.GetSize() > 1;
}

void Context::SetExcitation(Excitation* ex)
{
  // add check for sequence switch
  this->excitation_ = ex;
}

bool Context::IsStrainExcitedSystem() const
{
  // this shall not be called to often, hence we don't cache
  if(homogenization)
    return true;

  StdVector<LinearFormContext*>& lf = pde->GetAssemble()->GetLinForms();
  // ignore the regions!!
  for(unsigned int i = 0;i < lf.GetSize();i++)
    if (lf[i]->GetIntegrator()->GetName() == "AddStrainRHSInt")
    return true;
  return false;
}


void Context::SetMultiExcitations(MultipleExcitation* me, unsigned int basic_excitaions)
{
  this->me_ = me;
  this->basic_excitations_ = basic_excitaions;

  assert(excitations.GetSize() >= basic_excitaions);
  assert(!(me->DoMetaExcitation(this) && me->GetNumberMeta(this, true) * basic_excitaions != excitations.GetSize()));
  assert(!(!me->DoMetaExcitation(this) && excitations.GetSize() != basic_excitaions));
}

Excitation* Context::GetExcitation(unsigned int base, unsigned int meta)
{
  assert(base <= basic_excitations_);

  return excitations[basic_excitations_ * meta + base];
}

Excitation* Context::GetExcitation(unsigned int base, const std::string& meta)
{
  if(meta.empty() || !isdigit(meta[0]))
    throw Exception("excitation parameter '" + meta + "' is not a number");

  return GetExcitation(base, boost::lexical_cast<unsigned int>(meta));
}

Excitation* Context::GetExcitation(unsigned int base, Function* f)
{
  assert(f->ctxt == this);
  assert(base < excitations.GetSize());
  if(!me_->DoMetaExcitation(f->ctxt)){
    LOG_DBG3(context) << "C::GE  base=" << base << " excitation=" << excitations[base]->test_strain.ToString();
    return excitations[base];
  }
  else
    return excitations[basic_excitations_ * f->GetExcitation()->meta_index + base]; // * and + swapped??
}

App::Type Context::ToApp(const SinglePDE* pde)
{
  if(pde->GetName() == "electrostatic") return App::ELEC;
  if(pde->GetName() == "mechanic") 
  {
    if(pde->GetAnalysisType() == BasePDE::BUCKLING)
      return App::BUCKLING;
    else
      return App::MECH;
  }
  if(pde->GetName() == "heatConduction") return App::HEAT;
  if(pde->GetName() == "acoustic") return App::ACOUSTIC;
  if(pde->GetName() == "LatticeBoltzmann") return App::LBM;
  if(pde->GetName() == "magnetic") return App::MAG;
  if(pde->GetName() == "magneticEdge") return App::MAG;

  throw Exception("invalid");
}



SinglePDE* Context::ToPDE(App::Type app, bool throw_exception)
{
  std::map<App::Type, SinglePDE*>::const_iterator it = pdes.find(app);
  if(it != pdes.end())
    return it->second;

  // nothing found
  if(throw_exception)
    EXCEPTION("No PDE '" << app << "' stored");

  return NULL;
}

bool Context::SetPDEs()
{
  const StdVector<SinglePDE*> avail = domain->GetSinglePDEs();

  // might be empty for the first call within Optimization constructor. Here we set the design which we need
  // to setup the pdes with the proper optimization based material coefficients
  // -> call Context::Update() again in Optimization::PostInit()
  for(unsigned int i = 0; i < avail.GetSize(); i++)
  {
    const SinglePDE* sp = avail[i];

    if(sp->GetName() == "heatConduction") {
      pde = domain->GetSinglePDE("heatConduction", true);
      pdes[App::HEAT] = pde;
    }

    if(sp->GetName() == "magnetic") {
      pde = domain->GetSinglePDE("magnetic", true);
      pdes[App::MAG] = pde;
    }

    if(sp->GetName() == "magneticEdge") {
      pde = domain->GetSinglePDE("magneticEdge", true);
      pdes[App::MAG] = pde;
    }

    if(sp->GetName() == "acoustic") {
      pde = domain->GetSinglePDE("acoustic", true);
      pdes[App::ACOUSTIC] = pde;
    }

    // this sets elec for the piezo-case even if we want only to optimize mechanic
    if(sp->GetName() == "electrostatic") {
      pde = domain->GetSinglePDE("electrostatic", true);
      pdes[App::ELEC] = pde;
    }

    if(sp->GetName() == "LatticeBoltzmann") {
      pde = domain->GetSinglePDE("LatticeBoltzmann", true);
      pdes[App::LBM] = pde;
    }

    // mechanic last to have mechanic in the piezo-case as pde shortcut
    if(sp->GetName() == "mechanic") {
      assert(!(pde != NULL && domain->GetSinglePDE("mechanic", true) != pde));
      pde = domain->GetSinglePDE("mechanic", true);
      if (sp->GetAnalysisType() == BasePDE::BUCKLING)
        pdes[App::BUCKLING] = pde;
      else
        pdes[App::MECH] = pde;
      stt = pde->GetSubTensorType();
    }
  }

  assert(pdes.size() == avail.GetSize());

  return !avail.IsEmpty();
}

BiLinFormContext* Context::GetBiLinFormContext(const RegionIdType reg, App::Type app1, App::Type app2, bool throw_exception)
{
  App::Type a1, a2;
  string integrator = "";

  if(app1 == App::MECH && (app2 == App::MECH || app2 == App::NO_APP))
  {
    a1 = a2 = App::MECH;
    integrator = "LinElastInt";
  }
  if(app1 == App::BUCKLING && (app2 == App::BUCKLING || app2 == App::NO_APP))
  {
    a1 = a2 = App::BUCKLING;
    integrator = "PreStressInt";
  }
  if(app1 == App::ELEC && (app2 == App::ELEC || app2 == App::NO_APP))
  {
    a1 = a2 = App::ELEC;
    integrator = "LinGradBDBInt";
  }
  if((app1 == App::MECH && app2 == App::ELEC) || (app1 == App::ELEC && app2 == App::MECH) || (app1 == App::PIEZO_COUPLING && app2 == App::NO_APP))
  {
    a1 = App::MECH;
    a2 = App::ELEC;
    integrator = "LinPiezoCoupling";
  }
  if(app1 == App::MASS && (app2 == App::MASS || app2 == App::NO_APP))
  {
    a1 = a2 =App:: MASS;
    integrator = "MassInt";
  }

  assert(integrator != "");

  SinglePDE* pde1 = ToPDE(a1, throw_exception);
  SinglePDE* pde2 = ToPDE(a2, throw_exception);

  if(pde1 == NULL || pde2 == NULL)
  {
    if(!throw_exception)
      return NULL;
    else
      EXCEPTION("No PDE for application " << a1 << " resp. " << a2);
  }

  // in pre-FE-Space there was also the entry type!
  assert(pde2 != NULL); // otherwise it would be a LinearForm
  return pde1->GetAssemble()->GetBiLinForm(integrator, reg, pde1, pde2, !throw_exception);
}


ContextManager::ContextManager()
{
  this->initialized_ = false;
  this->default_excitation_ = new Excitation(); // for the load ersatz material case
  context.Resize(1); // default
  Optimization::context = &context[0];
  Optimization::context->SetExcitation(default_excitation_);
}

ContextManager::~ContextManager()
{
  delete default_excitation_;
}

void ContextManager::Init()
{
  if(domain->GetMultiSequenceDriver() == NULL)
  {
    SingleDriver* sd = domain->GetSingleDriver();
    assert(sd->GetDriverClass() == BaseDriver::SINGLE_DRIVER);
    assert(context.GetSize() == 1);
    context[0].Setup(this, sd->GetAnalysisType(), sd->GetParam(), 1);
  }
  else
  {
    MultiSequenceDriver* msd = domain->GetMultiSequenceDriver();
    assert(msd != NULL);
    std::map<UInt, BasePDE::AnalysisType> map = msd->GetAnalyisPerStep();
    std::map<unsigned int, PtrParamNode>  pns = msd->paramPerStep;

    unsigned int nms = msd->GetNumberOfSequenceSteps();

    assert(nms > 1);
    assert(map.size() == nms);
    assert(map.find(1) != map.end()); // the steps are one based and need to be consecutive
    assert(map.find(nms) != map.end()); // no idea why a map is used ?!
    assert(map.find(nms+1) == map.end()); // 1-based!

    context.Resize(nms);

    for(unsigned int i = 0; i < context.GetSize(); i++)
      context[i].Setup(this, map[i+1], pns[i+1]->Get("analysis")->GetChild(), i+1); // becomes an issue when analysis gets attributes
  }


  // SwitchContext(0) is too early for the multi sequence case. It triggers initialization of the pdes which triggers
  // creation of integrator which check for the existence of optimization with is too early as we are still in the
  // constructor of Optimization.
  Optimization::context = &context[0];
  Optimization::context->Update();

  // we have the default excitation to allow ersatz material only without optimization
  this->initialized_ = true;
}

ContextManager::AnyContext::AnyContext()
{
  this->bloch = false;
  this->complex = false;
  this->harmonic = false;
  this->eigenvalue = false;
  this->homogenization = false;
}

const ContextManager::AnyContext ContextManager::any() const
{
  AnyContext any;

  // because we don't know when harmonic is set we don't cache it

  for(unsigned int i = 0; i < context.GetSize(); i++)
  {
    const Context& c = context[i];
    any.bloch           = c.DoBloch() ? true : any.bloch;
    any.complex         = c.IsComplex() ? true : any.complex;
    any.harmonic        = c.IsHarmonic() ? true : any.harmonic;
    any.eigenvalue      = c.IsEigenvalue() ? true : any.eigenvalue;
    any.homogenization  = c.homogenization ? true : c.homogenization;
  }

  return any;
}

void ContextManager::SwitchContext(int index)
{
  if(domain->GetMultiSequenceDriver() != NULL)
  {
    // we keep the drivers and pdes of the older context
    for(unsigned int i = 0; i < context.GetSize(); i++)
      context[i].active = false;

    // detect if we already have the context. Special care for initial initialization, then no single driver is set yet
    if(domain->GetSingleDriver() == NULL || (int) domain->GetMultiSequenceDriver()->GetActSequenceStep() != index + 1)
      domain->GetMultiSequenceDriver()->SetSequenceStep(index + 1);
  }
  else
  {
    assert(domain->GetDriver()->GetDriverClass() == BaseDriver::SINGLE_DRIVER);
    assert(context.GetSize() == 1);
    assert(index == 0);
  }

  Optimization::context = &context[index];
  Optimization::context->Update();
  Optimization::context->active = true;
}

void ContextManager::SwitchContext(Excitation* excitation)
{
  assert(excitation->sequence >= 1); // 1-based
  SwitchContext(excitation->sequence -1); // 0-based
  Optimization::context->SetExcitation(excitation);
}

Context& ContextManager::GetContext(const Excitation* ex)
{
  return context[ex->sequence -1];
}

Context* ContextManager::GetHomogenization()
{
  for(unsigned int i = 0; i < context.GetSize(); i++)
    if(context[i].homogenization)
      return &(context[i]);

  // assume we have not more sequences with homogenization

  return NULL;
}
