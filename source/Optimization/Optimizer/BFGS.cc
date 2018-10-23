#include "Optimization/Optimizer/BFGS.hh"

#include "Utils/tools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ProgramOptions.hh"

DECLARE_LOG(BFGS)
DEFINE_LOG(BFGS, "bfgs")

using namespace CoupledField;
using std::pow;
using std::max;
using std::min;
using std::abs;
using std::string;

void BFGS::set_values (int x, int y) {
  a = x;
  b = y;
}

int BFGS::sum(){
  return a+b;
}
