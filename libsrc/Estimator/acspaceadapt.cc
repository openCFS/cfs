#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "acspaceadapt.hh"
#include "grid.hh"

namespace CoupledField
{

AcousticSpaceErrorEstimator::AcousticSpaceErrorEstimator(BasePDE * aptPDE, Grid * aptGrid)
:SpaceErrorEstimator(aptPDE,aptGrid)
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::AcousticTimeErrorEstimator" << std::endl;
#endif

  // read tolerance 
  conf->get("tol_sp",tol_,"Acoustic");
  conf->get("theta_sp", theta_, "Acoustic");
}

void AcousticSpaceErrorEstimator::CalcError()
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::CalcError" << std::endl;
#endif

 ;
}

Boolean AcousticSpaceErrorEstimator::TestError()
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::TestError" << std::endl;
#endif

  return TRUE;
}

void  AcousticSpaceErrorEstimator::RefineMesh()
{
 Integer levelGrid=ptGrid_->GetLastLevel();
 Integer numelems=ptGrid_->GetMaxnumElem(levelGrid);

 Integer i;
 for (i=0; i<numelems; i++) 
 {
   ptGrid_->SetRefinementFlag(i);  
 }

  ptGrid_->Refine();
}

} // end of namespace
