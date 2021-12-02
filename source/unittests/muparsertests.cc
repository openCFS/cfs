#include <boost/test/unit_test.hpp>
#include "muParser.h"
#include "Utils/mathParser/mathParser.hh"

using namespace CoupledField;

BOOST_AUTO_TEST_CASE(muparser_tests)
{

  try
  {
    double var_rho = .5;
    mu::Parser mp;
    mp.DefineVar("rho", &var_rho);
    mp.DefineConst("p", 1);
    mp.SetExpr("1-(1-rho)");

    for(double rho=0; rho <= 1.0; rho += .1)
    {
      var_rho = rho;
      std::cout << "muparser(" << rho << ")=" << mp.Eval() << std::endl;
    }
  }
  catch (mu::Parser::exception_type &e)
  {
    std::cout << e.GetMsg() << std::endl;
  }
}

