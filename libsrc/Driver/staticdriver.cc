#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"
#include "spaceerror.hh"

#include "outGMV.hh"

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

  ptdomain_->GetPDE(pdenumber)->SetMatrixFactors();
  ptdomain_->GetPDE(pdenumber)->SolveStepStatic(ptdomain_->GetBCs(), level);

  ptdomain_->PrintGrid(level);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
}

void StaticDriver :: SolveProblemAdaptSpace()
{
#ifdef TRACE
  (*trace) << "entering SolveProblemAdapt::SolveProblemAdapt " << std::endl;
#endif

  Integer level = 0;
  Integer pdenumber = 0;

  SpaceErrorEstimator * ptSpaceError;
  ptSpaceError = ptdomain_->GetPDE(pdenumber)->CreatePtSpaceError();

  ptdomain_->GetPDE(pdenumber)->SolveStepStatic(ptdomain_->GetBCs(), level);

  Integer maxnumrepeat, numrepeat = 0;
  conf->get("maxnumrepeat",maxnumrepeat,"SpaceAdaptivity");

  while (ptSpaceError->TestError() && numrepeat !=maxnumrepeat) {
    ptSpaceError->RefineMesh();

//   Char * name="laplace_test";
//   WriteResults * ptOut=new WriteResultsGMV(name);
//   ptOut->Init(ptdomain_->GetGrid());
//   ptOut->WriteGrid(0);
//   if (ptOut) delete ptOut; 

    ptdomain_->Update(level);
  
    ptdomain_->GetPDE(pdenumber)->SolveStepStatic(ptdomain_->GetBCs(), level);

    numrepeat++;
  }
 
  ptdomain_->PrintGrid(level);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
}

}
