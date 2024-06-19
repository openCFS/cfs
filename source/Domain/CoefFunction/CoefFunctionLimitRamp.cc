/*
 * CoefFunctionLimitRamp.cc
 *
 *  Created on: 30 April 2024
 *      Author: Dominik Mayrhofer
 */

#include "CoefFunctionLimitRamp.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

DEFINE_LOG(coeffctlimitramp, "coeffctlimitramp")

namespace CoupledField {
  CoefFunctionLimitRamp::CoefFunctionLimitRamp(shared_ptr<BaseFeFunction> feFct,
                                                StdVector<std::string> surfList,
                                                StdVector<Double> primaryOffset, 
                                                StdVector<Double> primaryPeakVal, 
                                                StdVector<Double> coefFuncPeakVal, 
                                                StdVector<bool> useMeanPres) : CoefFunction() {
    dimType_ = SCALAR; // should be determined from input
    dependType_ = SOLUTION; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;

    feFct_ = feFct;
    surfList_ = surfList;
    primaryOffset_ = primaryOffset;
    primaryPeakVal_ = primaryPeakVal;
    coefFuncPeakVal_ = coefFuncPeakVal;
    useMeanPres_ = useMeanPres;

    RegionIdType surfRegion;
    for( UInt i = 0; i < surfList_.size(); ++i ) {
      surfRegion = feFct_->GetGrid()->GetRegionId(surfList_[i]);
      surfRegions_.insert(surfRegion);
    }
  }

  CoefFunctionLimitRamp::~CoefFunctionLimitRamp() {
    ;
  }

  void CoefFunctionLimitRamp::GetScalar(Double& scal, const LocPointMapped& surfLpm ) {
    LOG_DBG(coeffctlimitramp) << "+++++ CoefFunctionLimitRamp::GetScalar ++++++";
    
    // Check that we are actually dealing with a surface lpm
    //assert(surfLpm.isSurface); this check fails since the evaluation by the result handler doesn't set the isSurface flag

    RegionIdType surfaceRegion = surfLpm.ptEl->regionId;

    std::set<RegionIdType>::iterator surfRegionsIt_ = surfRegions_.begin();

    bool foundSurf = false;
    for (UInt i = 0; i < surfList_.size(); ++i) {
      if( surfaceRegion==*surfRegionsIt_ ) {
        foundSurf = true;
        if( useMeanPres_[i] ){
          EXCEPTION("CoefFunctionLimitRamp::GetScalar not implemented for useMeanPres_==true");
        } else {
          // obtain avergae surface value used to define the BC at that surface element
          Double surfValAvg;
          const SurfElem* srcElem;
          srcElem = surfLpm.GetShapeMap()->GetSurfElem();
          feFct_->GetAvgElemValue(surfValAvg,srcElem);

          // check which ramp behaviour we have
          if( primaryPeakVal_[i]>=primaryOffset_[i] ) {
            // positive ramp
            // now check the three possible cases for the limited ramp
            if( surfValAvg<primaryOffset_[i] ) {
              scal = 0.0;
            } else if( surfValAvg>primaryPeakVal_[i] ) {
              scal = coefFuncPeakVal_[i];
            } else {
              // linear interpolation
              scal = (surfValAvg-primaryOffset_[i])/(primaryPeakVal_[i]-primaryOffset_[i])*coefFuncPeakVal_[i];
            } 
          } else {
            // negative ramp
            // now check the three possible cases for the limited ramp
            if( surfValAvg>primaryOffset_[i] ) {
              scal = 0.0;
            } else if( surfValAvg<primaryPeakVal_[i] ) {
              scal = coefFuncPeakVal_[i];
            } else {
              // linear interpolation
              scal = (surfValAvg-primaryOffset_[i])/(primaryPeakVal_[i]-primaryOffset_[i])*coefFuncPeakVal_[i];
            } 
          }
        }
      } else {
        // increment iterators
        surfRegionsIt_++;
      }
    }
    if( !foundSurf ){
      EXCEPTION("CoefFunctionLimitRamp: Could not find surface to evaluate, check your definition in the XML!");
    }

    // For debugging purposes
    LOG_DBG(coeffctlimitramp) << "Calculated value: " << scal << std::endl;
  }

  void CoefFunctionLimitRamp::CalcMeanPressure() {
    EXCEPTION("CoefFunctionLimitRamp::CalcMeanPressure() not implemented yet");
  }

}
