#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Optimization/SplineBoxOpt.hh"

namespace CoupledField {

DECLARE_LOG(ShOpt)
DEFINE_LOG(ShOpt, "shapeOpt")


SplineBoxOpt::SplineBoxOpt()
 : ParamMat()
{

}

SplineBoxOpt::~SplineBoxOpt()
{

}


void SplineBoxOpt::PostInit()
{
  ParamMat::PostInit();
}

}
