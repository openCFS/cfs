#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

namespace CoupledField
{

StaticDriver :: StaticDriver(Domain<Point2D> * adomain)
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


void StaticDriver :: SetupMatricesPDE(Integer pdenumber, Integer type)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SetUpMatricesPDE" << std::endl;
#endif
 
  ptdomain->GetPDE(pdenumber)->SetupMatrices(ptdomain->GetAlgSys(),type);

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

  ptdomain->GetPDE(pdenumber)->SolveStepStatic(ptdomain->GetAlgSys(), ptdomain->GetBCs(), level);

  //get solution
  Double *sol = ptdomain->GetAlgSys()->GetSolution(nsys);

  ptdomain->PrintSolution(sol,nsys);
}

}
