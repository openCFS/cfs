#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include "elecst2dErr.hh"
#include "elecst2dPDE.hh"

#ifdef ADAPTGRID
#include "interface_adgrid.hh"
#endif

namespace CoupledField
{

Elecst2dErrEstimator::Elecst2dErrEstimator(BasePDE * aptPDE, Grid * aptGrid): SpaceErrorEstimator(aptPDE,aptGrid)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimator::Elecst2dErrEstimator" << std::endl;
#endif

  // read tolerance 
    conf->get("tol_sp",tol_,"Elecst2d");
    conf->get("theta_sp", theta_, "Elecst2d");

  maxenergynormsol_=0.0;
}

void Elecst2dErrEstimator::CalcError()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimator::CalcError" << std::endl;
#endif

 error_=-1;
 ;
}

Double Elecst2dErrEstimator::CalcLocError(const Integer iElem)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimator::CalcLocError" << std::endl;
#endif

 Double locerr=0;
 return locerr;

}

#ifdef ADAPTGRID
Boolean Elecst2dErrEstimator::TestLocError(grd::Element * t)
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimator::TestLocError" << std::endl;
#endif

  return TRUE;
}
#endif

Boolean Elecst2dErrEstimator::TestError()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimator::TestError" << std::endl;
#endif
   
  return TRUE;

  //  Double normsol=ptPDE_->CalcEnergyNorm();
//    if (normsol>maxenergynormsol_) maxenergynormsol_=normsol;

//    std::cout << " norm sol " << normsol << std::endl;

//    CalcError();
   
//    if (InfoPrint)
//     (*infofile) << " error " << error_ << " " << tol_ << std::endl; 

//     std::cout << " error " << error_ << " " << tol_ << std::endl;

//     if (error_<= tol_) return FALSE;
//     else return TRUE;
}

void  Elecst2dErrEstimator::RefineMesh()
{
#ifdef TRACE
  (*trace)<<" entering Elecst2dErrEstimator::RefineMesh " << std::endl;
#endif
  Integer i;

  std::vector<std::string> * listSDs=ptPDE_->getSDsPDE();

  //  SetRefFlag f(this);
  
#ifdef ADAPTGRID
  SetRefFlagTest f;  
  for(i=0; i<(*listSDs).size(); i++)
    ptGrid_->forEachElemSd(f,(*listSDs)[i]);
  ptGrid_->Refine();
#endif

}

void Elecst2dErrEstimator::DefineRefinedElems(const Integer level, std::vector<Integer> & elems2ref)
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
