#include "CoefFunctionSurfaceForceBalance.hh"

#include <limits>

namespace CoupledField {

  CoefFunctionSurfaceForceBalance::CoefFunctionSurfaceForceBalance( PtrCoefFct coef, 
                                                      shared_ptr<ResultFunctor> fnct,
                                                      shared_ptr<BaseFeFunction> feFct,
                                                      StdVector<std::string> surfList1,
                                                      StdVector<std::string> surfList2,
                                                      StdVector<std::string> volumeList,
                                                      StdVector<Vector<Double>> dirVec ) : CoefFunction() {

  // this type of coefficient is never analytic
  isAnalytic_ = false;
  dimType_ = VECTOR;
  coef_ = coef;
  // this is a bit weird but we only have one functor to all results using this functor
  fnct_ = fnct;
  fnctResName_ = fnct->GetResultInfo()->resultType;
  feFct_ = feFct;
  volumeList_ = volumeList;
  surfList1_ = surfList1;
  surfList2_ = surfList2;
  dirVec_ = dirVec;

  scalingFactor_ = 1.0;

  RegionIdType volRegion;
  RegionIdType surfRegion1;
  RegionIdType surfRegion2;
  for (UInt i = 0; i < surfList1_.size(); ++i) {
    volRegion = feFct_->GetGrid()->GetRegionId(volumeList_[i]);
    volRegions_.insert(volRegion);
    surfRegion1 = feFct_->GetGrid()->GetRegionId(surfList1_[i]);
    surfRegions1_.insert(surfRegion1);
    surfRegion2 = feFct_->GetGrid()->GetRegionId(surfList2_[i]);
    surfRegions2_.insert(surfRegion2);
  }

  // Consistency check: due to the nature of insert we would skip multiple occurences of the same region which will lead to different lengths than the original vectors
  // this won't work with the iterators used in the GetVector() function (and won't be needed), hence, we just perform a check here and throw an exception if inconsistencies occour
  if (volRegions_.size() != volumeList_.size() || surfRegions1_.size() != surfList1_.size() || surfRegions2_.size() != surfList2_.size()) {
    throw Exception("A surface occours in multiple surface pairs. Currently, each surface can be only used once in a single surface pair.");
  }
 }

//! \copydoc CoefFunction::IsZero
bool CoefFunctionSurfaceForceBalance::IsZero() const {
  EXCEPTION("Not implemented");
}

//! \copydoc CoefFunction::ToString
std::string CoefFunctionSurfaceForceBalance::ToString() const {
  EXCEPTION("Not implemented");
}

void CoefFunctionSurfaceForceBalance::GetVector(Vector<Double>& coefVec,
                const LocPointMapped& lpm ) {
  // This is basically a wrapper where we evaluate the original coefFct and scale it with a scaling factor based on the functor

  // since we can't define the resultFunctors in the constructor since some variables do not exist yet, we have to outsorce this to the GetVector function call
  // check if this function gets called the first time

  if( firstFuncCall_ ) {
    this->DefineResFunctor();
    firstFuncCall_= false;
  }

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( volRegions_);

  //Confirm that output variable is defined as a vector
  assert(this->dimType_ == VECTOR);

  // get the current region
  RegionIdType surfaceRegionId = surfLpm.ptEl->regionId;

  // search for the region we balance with
  std::string surfaceRegion;
  std::string balanceRegion;
  Vector<Double> surfaceRes;
  Vector<Double> balanceRes;
  Vector<Double> curDirVec;
  RegionIdType balanceRegionId = -2000000000; // hardcoded, but that should never pose a problem
  std::set<RegionIdType>::iterator surfRegions1It_ = surfRegions1_.begin();
  std::set<RegionIdType>::iterator surfRegions2It_ = surfRegions2_.begin();
  for (UInt i = 0; i < surfRegions1_.size(); ++i) {
    if (surfaceRegionId==*surfRegions1It_) {
      balanceRegionId = *surfRegions2It_;
      surfaceRegion = surfList1_[i];
      balanceRegion = surfList2_[i];
      
      //resHandler_->UpdateResultType(fnctResName_); // more efficient to let the coupling scheme do this
      StdVector<shared_ptr<EntityList>> resEnt = fnct_->GetResultEnt();
      StdVector<Vector<Double>> resVec = fnct_->GetResultVec();

      for(UInt u = 0; u < resEnt.size(); ++u) {
        if (resEnt[u]->GetName()==surfaceRegion) {
          surfaceRes = resVec[u];
        }
        if (resEnt[u]->GetName()==balanceRegion) {
          balanceRes = resVec[u];
        }
      }

      curDirVec = dirVec_[i];
      break;
    }
    if( surfaceRegionId==*surfRegions2It_ ) {
      balanceRegionId = *surfRegions1It_;
      surfaceRegion = surfList2_[i];
      balanceRegion = surfList1_[i];

      //resHandler_->UpdateResultType(fnctResName_); // more efficient to let the coupling scheme do this
      StdVector<shared_ptr<EntityList>> resEnt = fnct_->GetResultEnt();
      StdVector<Vector<Double>> resVec = fnct_->GetResultVec();

      for(UInt u = 0; u < resEnt.size(); ++u) {
        if (resEnt[u]->GetName()==surfaceRegion) {
          surfaceRes = resVec[u];
        }
        if (resEnt[u]->GetName()==balanceRegion) {
          balanceRes = resVec[u];
        }
      }
      
      curDirVec = dirVec_[i];
      break;
    }
    surfRegions1It_++;
    surfRegions2It_++;
  }
  if( balanceRegionId==-2000000000 ) {
    EXCEPTION("Did not find a suitable surface match for CoefFunctionForceBalance!")
  }

  // in the first iteration in the first step it can happen, that not everything is defined yet (e.g. isNeeded_)
  // this can cause problems in forward coupled problems, hence, we have to catch this 


  // Recalculate scaling factor
  scalingFactor_ = CalcScalingFactor(surfaceRes,balanceRes,curDirVec);

  // Get Vector of passed coef
  Vector<Double> temp;
  coef_->GetVector(temp, surfLpm);

  coefVec = temp*scalingFactor_;
}

double CoefFunctionSurfaceForceBalance::CalcScalingFactor(Vector<Double> surfaceRes, Vector<Double> balanceRes, Vector<Double> curDirVec) {
  double sf;
  // calculate the residual dF = F_1+F_2
  Vector<Double> res = surfaceRes;
  res.Add(balanceRes);
  // divide by - --> -dF/2
  Vector<Double> negHalfRes = res;
  negHalfRes.ScalarDiv(-2.0);
  
  // subtract halfRes from our result -- > (F_1-dF/2)
  Vector<Double> diff = surfaceRes;
  diff.Add(negHalfRes);

  // calc (F_1-dF/2) \cdot eta
  Vector<Double> newF = diff;
  double newScaledF;
  newF.Inner(curDirVec, newScaledF);

  // calc F_1 \cdot eta
  Vector<Double> F = surfaceRes;
  double scaledF;
  F.Inner(curDirVec, scaledF);

  // calc ((F_1-dF/2) \cdot eta) / (F_1 \cdot eta)
  sf = newScaledF/scaledF;

  // if it's the first run the scaling factor will be zero/nan
  if( !(abs(sf)>0.1 && abs(sf)<10.0) || std::isnan(sf) ) {
    sf = 1;
    //WARN("Overwriting scaling factor of CoefFunctionForceBalance to 1 (out of bound)");
  }

  return sf;
} 

void CoefFunctionSurfaceForceBalance::DefineResFunctor()
{
  resHandler_ = feFct_->GetPDE()->GetDomain()->GetResultHandler();
  /* // read in all required functors
  resHandler_ = feFct_->GetPDE()->GetDomain()->GetResultHandler();
  fnctList1_.Resize(volRegions_.size());
  fnctList2_.Resize(volRegions_.size());
  for (UInt i = 0; i < surfList1_.size(); ++i) {
    fnctList1_[i] = resHandler_->GetResFunctor(fnctResName_, surfList1_[i]);
    fnctList2_[i] = resHandler_->GetResFunctor(fnctResName_, surfList2_[i]);
  } */
}
} // end of namespace
