#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"

#include "DataInOut/GMV/outGMV.hh"

#ifndef NEWBASEPDE
#include <PDE/basepde.hh>
#else
#include <PDE/newBasePDE.hh>
#endif //#ifndef NEWBASEPDE


namespace CoupledField
{

HarmonicDriver :: HarmonicDriver(Domain * adomain)
:BaseDriver(adomain)
{
#ifdef TRACE
  (*trace) << "entering HarmonicDriver::HarmonicDriver" << std::endl;
#endif

}

HarmonicDriver :: ~HarmonicDriver()
{
#ifdef TRACE
  (*trace) << "entering HarmonicDriver::~HarmonicDriver" << std::endl;
#endif

}

void HarmonicDriver :: SetupMatricesPDE(const Integer pdenumber, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering HarmonicDriver::SetUpMatricesPDE" << std::endl;
#endif

  // ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Commented out for testing purposes!!!!!!!!!!!!!
  //  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);
}

void HarmonicDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering HarmonicDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;
  ptdomain_->PrintGrid(level);
  if (PrintGridOnly)
      exit(0);

  ptdomain_->GetPDE(pdenumber)->SolveStepHarmonic(level);
  ptdomain_->GetPDE(pdenumber)->PostProcess(level);
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

}


} // end of namespace
