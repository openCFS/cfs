#include "Optimization/Context.hh"
#include "Optimization/Excitation.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/MultiSequenceDriver.hh"
#include "Optimization/Optimization.hh"

using namespace CoupledField;

Context::Context()
{
  // to be replaces by any current Excitation::Apply()
  excitation_ = NULL;

  driver = NULL;
  manager_ = NULL;
  sequence = -1;
  context_idx = -1;
  driver_steps_ = 0;
  analysis = BasePDE::NO_ANALYSIS;

  num_harm_freq = 0;
  num_bloch_wave_vectors = 0;

  complex_ = false;
  harmonic_ = false;
  eigenvalue_ = false;
  bloch_ = false;
  pde = NULL;

}

Context::~Context()
{
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
    bloch_ = EigenFrequencyDriver::DoBloch(node);
    num_bloch_wave_vectors = EigenFrequencyDriver::GetNumBlochWave(node); // 0 when not bloch
    break;

  case BasePDE::STATIC:
    // FIXME num_static_loads not that easy to determine w/o pde. allows 0 (explixit excitions or homogenization) or more
    break;

  case BasePDE::TRANSIENT:
    assert(false); // not yet implemented
    break;

  case BasePDE::MULTI_SEQUENCE:
  case BasePDE::NO_ANALYSIS:
    assert(false);
    break;
  }
}

void Context::Update()
{
  driver = domain->GetSingleDriver();
  SetPDEs();
}

EigenFrequencyDriver* Context::GetEigenFrequencyDriver()
{
  return dynamic_cast<EigenFrequencyDriver*>(driver);
}

HarmonicDriver* Context::GetHarmonicDriver()
{
  return dynamic_cast<HarmonicDriver*>(driver);
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


App::Type Context::ToApp(const SinglePDE* pde)
{
  if(pde->GetName() == "electrostatic") return App::ELEC;
  if(pde->GetName() == "mechanic") return App::MECH;
  if(pde->GetName() == "heatConduction") return App::HEAT;
  if(pde->GetName() == "acoustic") return App::ACOUSTIC;
  if(pde->GetName() == "LatticeBoltzmann") return App::LBM;

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


void Context::SetPDEs()
{
  // re-read pdes. Might change for multiple sequence steps
  pdes.clear();

  const StdVector<SinglePDE*> avail = domain->GetSinglePDEs();
  assert(!avail.IsEmpty());

  for(unsigned int i = 0; i < avail.GetSize(); i++)
  {
    const SinglePDE* sp = avail[i];

    if(sp->GetName() == "heatConduction") {
      pde = domain->GetSinglePDE("heatConduction", true);
      pdes[App::HEAT] = pde;
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
      pde = domain->GetSinglePDE("mechanic", true);
      pdes[App::MECH] = pde;
    }
  }

  assert(pdes.size() == avail.GetSize());
}



ContextManager::ContextManager()
{
  this->initialized_ = false;
  context.Resize(1); // default
  Optimization::context = &context[0];
}

ContextManager::~ContextManager()
{
}

void ContextManager::Init()
{
  if(domain->GetMultiSequenceDriver() == NULL)
  {
    assert(domain->GetDriver()->GetDriverClass() == BaseDriver::SINGLE_DRIVER);
    assert(context.GetSize() == 1);
    context[0].Setup(this, domain->GetSingleDriver()->GetAnalysisType(), domain->GetDriver()->GetParam(), 1);
  }
  else
  {
    MultiSequenceDriver* msd = domain->GetMultiSequenceDriver();
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
      context[i].Setup(this, map[i+1], pns[i+1], i+1);
  }

  for(unsigned int i = 0; i < context.GetSize(); i++)
  {
    Context& c = context[i];
    any_.bloch      = c.DoBloch() ? true : any_.bloch;
    any_.complex    = c.IsComplex() ? true : any_.complex;
    any_.harmonic   = c.IsHarmonic() ? true : any_.harmonic;
    any_.eigenvalue = c.IsEigenvalue() ? true : any_.eigenvalue;
  }

  SwitchContext(0);

  // we have no Excitations yet in this phase
  assert(Optimization::context->GetExcitation() == NULL);

  this->initialized_ = true;
}

void ContextManager::SwitchContext(int index)
{
  if(domain->GetMultiSequenceDriver() != NULL)
  {
    // remove all drivers, we create a new one
    for(unsigned int i = 0; i < context.GetSize(); i++)
      context[i].driver = NULL; // we don't have the ownership

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
}

