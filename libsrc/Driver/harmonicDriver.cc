#include <fstream>
#include <iostream>
#include <string>

#include "harmonicDriver.hh"

#include "DataInOut/GMV/outGMV.hh"

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

  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);
}

void HarmonicDriver :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering HarmonicDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;

  ptdomain_->GetPDE(pdenumber)->SolveStepHarmonic(level);
  //  std::cout << "Solve Step ok" << std::endl;
  ptdomain_->GetPDE(pdenumber)->PostProcess(level);
  //  std::cout << "Post ok" << std::endl;
  ptdomain_->PrintGrid(level);
  //  std::cout << "Print Grid ok" << std::endl;
  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();

}


} // end of namespace
