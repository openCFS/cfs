#include "CoefFunctionVectorFromScalar.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

#include <limits>

namespace CoupledField {

  CoefFunctionVectorFromScalar::CoefFunctionVectorFromScalar(shared_ptr<BaseFeFunction> feFct,
                                                              StdVector<std::string> surfList,
                                                              StdVector<PtrCoefFct>& scalVals, 
                                                              StdVector<UInt> dofInd, 
                                                              UInt size): CoefFunction() {

    // this type of coefficient is never analytic
    isAnalytic_ = false;
    dofInd_ = dofInd;
    // we always have a vector
    dimType_ = VECTOR;

    // currently we only use this in combination with coefFunctions that explicitly depend on the solution, so we set this by hand
    // since the dependType is scalar and we get vectorial input, this would be problematic otherwise
    dependType_ = SOLUTION; 
    isComplex_ = false; // currently we only use this in transient simulations, not sure if it even makes sense for harmonic simulations
    
    surfList_ = surfList;

    // convert the string to a RegionIdType used for the evaluation in GetVector()
    RegionIdType surfRegion;
    for (UInt i = 0; i < surfList_.size(); ++i) {
      surfRegion = feFct->GetGrid()->GetRegionId(surfList_[i]);
      surfRegions_.insert(surfRegion);
    }
    
    size_ = size;

    mp_ = domain->GetMathParser();
    constZero_ = CoefFunction::Generate( mp_, Global::REAL, "0.0");

    scalVals_ = scalVals;
  }

//! \copydoc CoefFunction::GetVector
void CoefFunctionVectorFromScalar::GetVector( Vector<Double>& coefVec, const LocPointMapped& lpm ) {

  // currently we only allow surface evaluation since we only deal with a surface-based input (surfList)
  //assert(lpm.isSurface); this check fails since the evaluation by the result handler doesn't set the isSurface flag

  RegionIdType surfaceRegion = lpm.ptEl->regionId;

  //Confirm that output variable is defined as a vector
  assert(this->dimType_ == VECTOR);

  // get the correct index for the evaluation
  std::set<RegionIdType>::iterator surfRegionsIt_ = surfRegions_.begin();
  for (UInt i = 0; i < surfRegions_.size(); ++i) {

    if( surfaceRegion==*surfRegionsIt_ ) {
      // we have a match!

      // Resize and assign values
      coefVec.Resize(size_, 1);
      coefVec.Init();

      // get the actual non-zero value
      scalVals_[i]->GetScalar(coefVec[dofInd_[i]], lpm);
    } else {
      EXCEPTION("Did not find a suitable surface match for CoefFunctionVectorFromScalar::GetVector()!")
    }

    // increment iterator
    surfRegionsIt_++;
  }
}

//! \copydoc CoefFunction::IsZero
bool CoefFunctionVectorFromScalar::IsZero() const {
  EXCEPTION("CoefFunctionVectorFromScalar::IsZero() not implemented yet!");
}

//! \copydoc CoefFunction::ToString
std::string CoefFunctionVectorFromScalar::ToString() const {
  std::string ret = "";
  for( UInt ii = 0; ii < surfList_.size(); ++ii ) {
    ret += surfList_.ToString();
    ret += ": VectorFromScalar(";
    for( UInt i = 0; i < size_-1; ++i ) {
      if( dofInd_[ii]==i ) {
        ret += scalVals_[ii]->ToString();
      } else {
        ret += "0";
      }
      ret += ", ";
    }
    if( dofInd_[ii]==(size_-1) ) {
      ret += scalVals_[size_-1]->ToString();
    } else {
      ret += "0";
    }
    ret += ")\n";
  }
  return ret;
}
} // end of namespace
