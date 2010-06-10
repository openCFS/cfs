#include "Optimization/Filter.hh"

using namespace CoupledField;

Filter::Filter()
{
  type_        = NO_FILTERING;
  sensitivity_ = SIGMUND;
  density_     = STANDARD;
  beta_        = -1.0;
}


