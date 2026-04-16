/*
 * CoefFunctionContactForceDensity.cc
 *
 *  Created on: 18 Nov 2022
 *      Author: Dominik Mayrhofer
 */
#include "Domain/Domain.hh"
#include "Utils/mathParser/mathParser.hh"
#include "CoefFunctionContactForceDensity.hh"
#include "CoefFunctionDistance.hh"

namespace CoupledField {
  CoefFunctionContactForceDensity::CoefFunctionContactForceDensity(shared_ptr<BaseFeFunction> feFct,
                                                        StdVector<std::string> surfList1,
                                                        StdVector<std::string> surfList2,
                                                        StdVector<std::string> volumeList,
                                                        StdVector<std::string> contactLawList,
                                                        StdVector<bool> useSurfaceMidpointsList) : CoefFunction() {
    dimType_ = VECTOR; // should be determined from input
    dependType_ = SOLUTION; // should be determined from input
    isAnalytic_ = false;
    isComplex_  = false;

    // Functionality
    // We pass the feFunction, two surface regions as well as the analytic contact law to the coefFunction
    // Afterwards, we use NN search on elements to determine the actual penetration depth and to compute element based forceDensities
    feFct_ = feFct;
    surfList1_ = surfList1;
    surfList2_ = surfList2;
    volumeList_ = volumeList;
    contactLawList_ = contactLawList;
    useSurfaceMidpointsList_ = useSurfaceMidpointsList;

    RegionIdType volRegion;
    RegionIdType surfRegion1;
    RegionIdType surfRegion2;
    for (UInt i = 0; i < volumeList_.size(); ++i) {
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
      throw Exception("A surface or volume region occours in multiple contact-laws. Currently, only one contact-law per region is supported.");
    }
    

    // fetch math parser object
    mParser_ = domain->GetMathParser();


    // init function stuff (from postproc function SetResult)

    // intialize rVals and mathParser handles
    rVals_.Resize( contactLawList_.size() );
    rVals_.Init();
    rHandles_.Resize( contactLawList_.size() );
    coefs_.Resize( contactLawList_.size() );

#ifdef USE_OPENMP
    if(omp_get_num_threads()!=1){
      EXCEPTION("Call from parallel region which is not safe.")
    }
#endif

    for( UInt i = 0 ; i < contactLawList_.size(); i++ ) {
      // obtain new handle
      rHandles_[i] = mParser_->GetNewHandle(true);

      // register input 
      mParser_->SetValue( rHandles_[i], "u", 0 );      
        
      // register dof with dummy variable and set expression
      mParser_->SetExpr( rHandles_[i], contactLawList_[i] );

      shared_ptr<CoefFunctionDistance> coef;
      coef.reset(new CoefFunctionDistance(feFct,surfList1_[i],surfList2_[i],volumeList_[i],useSurfaceMidpointsList_[i]));
      coefs_[i] = coef;
    }

  }

  CoefFunctionContactForceDensity::~CoefFunctionContactForceDensity() {
    ;
  }

  void CoefFunctionContactForceDensity::GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) {
    
    
    // create local point for surface
    LocPointMapped surfLpm(lpm);
    surfLpm.SetSurfInfo( volRegions_);
    RegionIdType region = surfLpm.lpmVol->ptEl->regionId;

    RegionIdType surfaceRegion = surfLpm.ptEl->regionId;

    //Confirm that output variable is defined as a vector
    assert(this->dimType_ == VECTOR);

    // get the correct index for the evaluation of the contact law
    std::set<RegionIdType>::iterator volRegionsIt_ = volRegions_.begin();
    std::set<RegionIdType>::iterator surfRegions1It_ = surfRegions1_.begin();
    std::set<RegionIdType>::iterator surfRegions2It_ = surfRegions2_.begin();
    for (UInt i = 0; i < volumeList_.size(); ++i) {

      if( region==*volRegionsIt_ && (surfaceRegion==*surfRegions1It_ || surfaceRegion==*surfRegions2It_) ) {
        // get the scaling factor based on the contact law from the math parser

        Vector<Double> globMidPoint;
        surfLpm.GetShapeMap()->GetGlobMidPoint(globMidPoint);

        mParser_->SetCoordinates( rHandles_[i],*(domain->GetCoordSystem()),globMidPoint );

        // Register current input values for dofNames
        Double actVal;
        coefs_[i]->GetScalar(actVal, surfLpm);
        mParser_->SetValue( rHandles_[i], "u", actVal );
        
        // Apply function
        Double fac;
        fac = mParser_->Eval( rHandles_[i] );

        // Return force in surface normal direction
        vec = surfLpm.normal*fac;

      } else {
        EXCEPTION("Did not find a suitable volume/surface match for CoefFunctionContactForceDensity!")
      }

      // increment iterators
      volRegionsIt_++;
      surfRegions1It_++;
      surfRegions2It_++;
    }
  }
}
