#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Optimization/SplineBoxOpt.hh"

namespace CoupledField {

DECLARE_LOG(SBOpt)
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
