#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

namespace CoupledField
{

StaticDriver :: StaticDriver(Domain * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::StaticDriver" << std::endl;
#endif

}

StaticDriver :: ~StaticDriver()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::~StaticDriver" << std::endl;
#endif

}

void StaticDriver :: SetupMatricesPDE(const Integer pdenumber, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SetUpMatricesPDE" << std::endl;
#endif

  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);
}

void StaticDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;
  Integer matrixtype = 0;
  Integer nsys = 0;

  ptdomain_->GetPDE(pdenumber)->SolveStepStatic(ptdomain_->GetBCs(), level);

  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
}

}
