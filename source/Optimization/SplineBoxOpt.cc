#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/SplineBoxOpt.hh"

namespace CoupledField {

DEFINE_LOG(SBOpt, "splineBoxOpt")


SplineBoxOpt::SplineBoxOpt()
 : SIMP()
{

}

SplineBoxOpt::~SplineBoxOpt()
{

}


void SplineBoxOpt::PostInit()
{
  SIMP::PostInit();
}

}
