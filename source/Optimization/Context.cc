#include "Optimization/Context.hh"
#include "Optimization/Excitation.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/HarmonicDriver.hh"
#include "Driver/MultiSequenceDriver.hh"
#include "Optimization/Optimization.hh"

using namespace CoupledField;

Context::Context()
{
  // no instance to avoid the include in Context.hh
  default_excitation_ = new Excitation(); // has no assemble yet!

  // to be replaces by any current Excitation::Apply()
  excitation_ = default_excitation_;

  // we are a globally initialized object, there is no domain yet!
  assert(domain == NULL);
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

}

Context::~Context()
{
  delete default_excitation_;
  default_excitation_ = NULL;
}

void Context::Init(ContextManager* manager, BasePDE::AnalysisType analyis, PtrParamNode node, int sequence_step)
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

ContextManager::ContextManager()
{
  this->initialized_ = false;
  context.Resize(1); // default
  Optimization::context = &context[0];
}


void ContextManager::Init()
{
  if(domain->GetDriver()->GetDriverClass() == BaseDriver::SINGLE_DRIVER)
  {
    assert(context.GetSize() == 1);
    context[0].Init(this, domain->GetSingleDriver()->GetAnalysisType(), domain->GetDriver()->GetParam(), 1);
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
    assert(map.find(nms+1) != map.end()); // no idea why a map is used ?!

    context.Resize(nms);

    for(unsigned int i = 0; i < context.GetSize(); i++)
      context[i].Init(this, map[i+1], pns[i+1], i+1);
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
  this->initialized_ = true;
}

void ContextManager::SwitchContext(int index)
{
  if(domain->GetDriver()->GetDriverClass() == BaseDriver::SINGLE_DRIVER)
  {
    assert(context.GetSize() == 1);
    context[0].driver = domain->GetSingleDriver();
  }
  else
  {
    // remove all drivers, we create a new one
    for(unsigned int i = 0; i < context.GetSize(); i++)
      context[i].driver = NULL; // we don't have the ownership

    domain->GetMultiSequenceDriver()->SetSequenceStep(index + 1);
  }

  Optimization::context = &context[index];
}

