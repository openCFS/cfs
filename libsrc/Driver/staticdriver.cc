#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

#include "DataInOut/GMV/outGMV.hh"
#include <CoupledPDE/basecoupledpde.hh>
#include <General/environment.hh>

#include <PDE/basePDE.hh>

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

  const Integer nstep = 1;
  Double  steptime = 0;
  Boolean reset = FALSE;
  
  if (ptdomain_->GetNumPDE() <= 1) 
    {
      ptdomain_->GetPDE(pdenumber)->WriteGeneralPDEdefines();
      ptdomain_->GetPDE(pdenumber)->SolveStepStatic(nstep, steptime, level, reset);   
      ptdomain_->GetPDE(pdenumber)->PostProcess(level);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
    }
  else
    {
      ptdomain_->GetCoupledPDE()->WriteGeneralPDEdefines();
      ptdomain_->GetCoupledPDE()->SolveStepStatic(nstep, steptime, level, reset);
      ptdomain_->PrintGrid(level);
      ptdomain_->GetCoupledPDE()->WriteResultsInFile();
    }
}

} // end of namespace
