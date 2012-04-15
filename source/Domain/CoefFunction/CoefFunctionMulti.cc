#include "CoefFunctionMulti.hh"

namespace CoupledField  {


CoefFunctionMulti::CoefFunctionMulti() {

  // set global data
  dimType_ = CoefFunction::NO_DIM;
  dependType_ = CoefFunction::GENERAL;
  
  // a ditributed coefficient function can never be analytic
  isAnalytic_ = false;
  isComplex_ = false;
}

CoefFunctionMulti::~CoefFunctionMulti() {
  // clear map
  regionCoefs_.clear();
}

void CoefFunctionMulti::AddRegion( RegionIdType region, PtrCoefFct coef ) {
  // check, if this is the first entry
  if( regionCoefs_.size() == 0 ) {
    dimType_ = coef->GetDimType();
    isComplex_ = coef->IsComplex();
  } else {
    PtrCoefFct first = regionCoefs_.begin()->second;
    if( coef->GetDimType() != first->GetDimType() ) {
      EXCEPTION( "The dimensionality of the coefficient functions "
          << "is not the same");
    }
    if( isComplex_ != first->IsComplex() ) {
      EXCEPTION( "All coefficient functions must have the same complexType");
    }
  }
  
  // in the end, check if there was already a coefficient function 
  // for this region
  if( regionCoefs_.find( region ) != regionCoefs_.end() ) {
    EXCEPTION( "There was already a coefficient function defined for "
               << "the region width id " << region );
  }
  
  regionCoefs_[region] = coef;
}

void CoefFunctionMulti::GetTensor(Matrix<Complex>& coefMat,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetTensor(coefMat, lpm);
}

void CoefFunctionMulti::GetVector(Vector<Complex>& coefVec,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetVector(coefVec, lpm);
}

void CoefFunctionMulti::GetScalar(Complex& coef,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetScalar(coef, lpm);

}


void CoefFunctionMulti::GetTensor(Matrix<Double>& coefMat,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetTensor(coefMat, lpm);
}

void CoefFunctionMulti::GetVector(Vector<Double>& coefVec,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetVector(coefVec, lpm);

}

void CoefFunctionMulti::GetScalar(Double& coef,
                                  const LocPointMapped& lpm ) {
  RegionIdType curRegion = lpm.ptEl->regionId;
  return GetRegionCoef(curRegion)->GetScalar(coef, lpm);
}

UInt CoefFunctionMulti::GetVecSize() const {
  assert(regionCoefs_.size());
  return regionCoefs_.begin()->second->GetVecSize();
}

void CoefFunctionMulti::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  assert(regionCoefs_.size());
  return regionCoefs_.begin()->second->GetTensorSize(numRows, numCols);
}

std::string CoefFunctionMulti::ToString() const {
  std::stringstream ret;
  ret << "CoefFunctionMulti defined on:\n";
  std::map<RegionIdType,PtrCoefFct >::const_iterator it = regionCoefs_.begin();
  for( ; it != regionCoefs_.end(); ++it ) {
    ret << "regionId " << it->first << ", value:" << it->second->ToString() << std::endl;
  }
  return ret.str();
}

bool CoefFunctionMulti::IsComplex() {
  if( regionCoefs_.size() == 0 ) {
    WARN("No entries set yet");
    return false;
  }
  return isComplex_;
}



} // end of namespace
