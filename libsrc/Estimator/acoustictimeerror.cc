#include <iostream>
#include <fstream>
#include <string>

#include "acoustictimeerror.hh"

namespace CoupledField
{

AcousticTimeErrorEstimator::AcousticTimeErrorEstimator(BasePDE * aptPDE)
:TimeErrorEstimator(aptPDE)
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimator::AcousticTimeErrorEstimator" << std::endl;
#endif
  
}

void AcousticTimeErrorEstimator::ChangeStep(Double & dt)
{
#ifdef TRACE
  (*trace) << "entering TimeErrorEstimator::ChangeStep" << std::endl;
#endif

  Double tol, theta;

  conf->get("tol",tol,"Acoustic");
  conf->get("theta", theta, "Acoustic");

 dt*=sqrt(theta*tol/relativeerror);

}

/*
void TimeErrorEstimator::ChangeStep(Double & dt);
{
#ifdef TRACE
  (*trace) << "entering TimeErrorEstimator::ChangeStep" << std::endl;
#endif
 
  conf->get("tol",tol,"Acoustic");
  conf->get("theta", theta, "Acoustic");

 dt*=sqrt(theta*tol/relativeerror); 
   
}

Boolean TimeErrorEstimator::TestError()
{
#ifdef TRACE
  (*trace) << "entering TimeErrorEstimator::TestError" << std::endl;
#endif

  

}
*/

} // end of namespace
