#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

#include "DataInOut/GMV/outGMV.hh"
#include <PDE/basepde.hh>
#include <CoupledPDE/basecoupledpde.hh>
#include <General/environment.hh>

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

void StaticDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SolveProblem" << std::endl;
#endif


  Integer level=0;
  Integer pdenumber  = 0;

  if (PrintGridOnly)
    {
      ptdomain_->PrintGrid(level);
      exit(0);
    }

  if (ptdomain_->GetNumPDE() <= 1) 
    {
      ptdomain_->GetPDE(pdenumber)->SolveStepStatic(level);   
      ptdomain_->GetPDE(pdenumber)->PostProcess(level);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
    }
  else
    {
      ptdomain_->GetCoupledPDE()->SolveStepStatic(level);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetCoupledPDE()->WriteResultsInFile();
    }
  
}

} // end of namespace
