#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "acspaceadapt.hh"
#include "grid.hh"
#include "acoustic2dPDE.hh"

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

  maxenergynormsol_=1.0;
}

void AcousticSpaceErrorEstimator::CalcError()
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::CalcError" << std::endl;
#endif

 error_=-1;
 ;
}

Double AcousticSpaceErrorEstimator::CalcLocError(const Integer iElem)
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::CalcLocError" << std::endl;
#endif

 Double locerr=0;
 return locerr;

}

Boolean AcousticSpaceErrorEstimator::TestError()
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimator::TestError" << std::endl;
#endif
   
   Double normsol=ptPDE_->CalcEnergyNorm();
   if (normsol>maxenergynormsol_) maxenergynormsol_=normsol;

   std::cout << " norm sol " << normsol << std::endl;

   CalcError();
   
   if (InfoPrint)
    (*infofile) << " error " << error_ << " " << tol_ << std::endl; 

    std::cout << " error " << error_ << " " << tol_ << std::endl;

    if (error_<= tol_) return FALSE;
    else return TRUE;
}

void  AcousticSpaceErrorEstimator::RefineMesh()
{

 Integer levelGrid=ptGrid_->GetLastLevel();

 std::vector<Integer> elems2refinement;
 DefineRefinedElems(levelGrid, elems2refinement);

 ptGrid_->SetRefinementFlag(elems2refinement);

 ptGrid_->Refine();
}

void AcousticSpaceErrorEstimator::DefineRefinedElems(const Integer level, std::vector<Integer> & elems2ref)
{
  Integer numelems=ptGrid_->GetMaxnumElem(level);

  Double locerr_tol=theta_*tol_*maxenergynormsol_/sqrt(numelems);
  Double locerr;

  Integer i;
  for (i=0; i<numelems; i++)
   {
     locerr=CalcLocError(i);
     if (locerr > locerr_tol) elems2ref.push_back(i);     
   }
}

} // end of namespace
