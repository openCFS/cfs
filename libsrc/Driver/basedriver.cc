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
}

BaseDriver :: ~BaseDriver()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::~BaseDriver" << std::endl;
#endif
  if (ptdomain) delete ptdomain;

}


void BaseDriver :: SetupMatricesPDE(Integer pdenumber)
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::SetUpMatricesPDE" << std::endl;
#endif
 
  ptdomain->ptpde[pdenumber]->AssembleGlobal(ptdomain->ptalgsys);

}

void BaseDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering BaseDriver::SolveProblem" << std::endl;
#endif
 
  Integer pdenumber = 0;
  SetupMatricesPDE(pdenumber);
}

}
