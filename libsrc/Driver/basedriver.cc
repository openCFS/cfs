#include <fstream>
#include <iostream>
#include <string>

#include "basedriver.hh"

namespace CoupledField
{

BaseDriver :: BaseDriver(Domain<Point2D> * adomain)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::BaseDriver" << std::endl;
#endif
  ptdomain = adomain;

    // read info Should we save first,second derivatives or not
  ptdomain->GetInFile()->ReadOutputOptions(SaveDer1,SaveDer2);

}

BaseDriver :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif

}

/*
void BaseDriver :: WriteResultsInFile(BasePDE * ptPDE, const Integer step, const Double t)
{
 
  ptdomain->GetOutFile()->WriteSolution(ptPDE->getS(),step,t);

  if (SaveDer1)
  ptdomain->GetOutFile()->WriteFirstDerSolution(ptPDE->getS1(),step,t);

  if (SaveDer2)
  ptdomain->GetOutFile()->WriteSecondDerSolution(ptPDE->getS2(),step,t);
}
*/

}
