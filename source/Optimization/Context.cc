#include "Optimization/Context.hh"
#include "Optimization/Excitation.hh"

using namespace CoupledField;

Context::Context()
{
  // no instance to avoid the include in Context.hh
  default_excitation_ = new Excitation(); // has no assemble yet!

  // to be replaces by any current Excitation::Apply()
  excitation_ = default_excitation_;

  // we are a globally initialized object, there is no domain yet!
  assert(domain == NULL);
  driver_ = NULL;

}

Context::~Context()
{
  delete default_excitation_;
  default_excitation_ = NULL;
}

void Context::Init()
{
  assert(!domain->GetMultiSequenceDriver());
  driver_ = dynamic_cast<SingleDriver*>(domain->GetDriver());
  assert(driver_ != NULL);

  // we have a driver, therefore we can parse it
  ParseSequence();
}

void Context::SetExcitation(Excitation* ex)
{
  // add check for sequence switch
  this->excitation_ = ex;
}

void Context::ParseSequence()
{
  // even a standard EV problem is in CFS complex (with imag=0). For Bloch we are complex anyways
  complex_ = driver_->GetAnalysisType() == BasePDE::HARMONIC || driver_->GetAnalysisType() == BasePDE::EIGENFREQUENCY;
  harmonic_ = driver_->GetAnalysisType() == BasePDE::HARMONIC;
  eigenvalue_ = driver_->GetAnalysisType() == BasePDE::EIGENFREQUENCY;
  bloch_ = driver_->DoBlochModeEigenfrequency();

}

