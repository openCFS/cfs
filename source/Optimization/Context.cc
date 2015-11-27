#include "Optimization/Context.hh"
#include "Optimization/Excitation.hh"

using namespace CoupledField;

Context::Context()
{
  default_excitation_ = new Excitation();

  // to be replaces by any current Excitation::Apply()
  excitation = default_excitation_;
}

Context::~Context()
{
  delete default_excitation_;
}



