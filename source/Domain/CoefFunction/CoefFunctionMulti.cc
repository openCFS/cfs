#include "CoefFunctionMulti.hh"
#include "CoefFunctionConst.hh"

namespace CoupledField  {


CoefFunctionMulti::CoefFunctionMulti( CoefDimType dimType,
                                      UInt dim1, UInt dim2,
                                      bool isComplex,
                                      bool zeroEmptyRegions ) {

  // set global data
  dimType_ = dimType;
  dependType_ = CoefFunction::GENERAL;
  
  // a distributed coefficient function can never be analytic
  isAnalytic_ = false;
  isComplex_ = isComplex;
  rowSize_ = dim1;
  colSize_ = dim2;
  
  zeroEmptyRegions_ = zeroEmptyRegions;
}

CoefFunctionMulti::~CoefFunctionMulti() {
  // clear map
  regionCoefs_.clear();
}

void CoefFunctionMulti::AddRegion( RegionIdType region, PtrCoefFct coef, bool allowReplacement ) {
	// check, if this is the first entry
  if( regionCoefs_.size() == 0 ) {
    shared_ptr<CoefFunctionConst<Complex> > cFct(new CoefFunctionConst<Complex>());
    shared_ptr<CoefFunctionConst<Double> > rFct(new CoefFunctionConst<Double>());
    // generate empty coefficient functions
    if( dimType_ == CoefFunction::SCALAR) {
      cFct->SetScalar(0.0);
      rFct->SetScalar(0.0);
    } else if( dimType_ == CoefFunction::VECTOR ) {
      Vector<Complex> cZvec(coef->GetVecSize());
      Vector<Double> rZvec(coef->GetVecSize());
      cZvec.Init(0.0);
      rZvec.Init(0.0);
      cFct->SetVector(cZvec);
      rFct->SetVector(rZvec);
    } else if( dimType_ == CoefFunction::TENSOR ) {
      UInt numRows, numCols;
      coef->GetTensorSize( numRows, numCols);
      Matrix<Complex> cTens(numRows,numCols);
      Matrix<Double> rTens(numRows,numCols);
      cTens.Init();
      rTens.Init();
      cFct->SetTensor(cTens);
      rFct->SetTensor(rTens);
    } else  {
      EXCEPTION( "Unknown dimension type" );
    }

    if(isComplex_) {
      zeroCoef_ = cFct;
    } else {
      zeroCoef_ = rFct;
    }

  } else {
    PtrCoefFct first = regionCoefs_.begin()->second;

    if( coef->GetDimType() != dimType_ ) {
      EXCEPTION( "The dimensionality of the coefficient functions "
          << "is not the same");
    }

    // check size of tensor
    UInt nRows, nCols;
    switch (coef->GetDimType() ) {
    case VECTOR:
      /* if(coef->GetVecSize() != rowSize_ ) {
        EXCEPTION( "Vector size inconsistent: "<<"coef->GetVecSize() = "<<coef->GetVecSize()<<"; rowSize_="<<rowSize_ );
      } */
      break;
    case TENSOR:
      coef->GetTensorSize(nRows, nCols);
      if( nRows != rowSize_ || nCols != colSize_) {
        EXCEPTION( "Tensor size inconsistent" );
      }
      break;
    default:
      break;
    }
    if( isComplex_ != coef->IsComplex() ) {
      EXCEPTION("Developer info: You specified a coefficient function of complexType: " << coef->IsComplex() <<
                " for a CoefFunctionMulti of complexType: " << isComplex_ << "." << std::endl <<
                "Make sure all coefficient have the same complexType!");
    }
  }

  // in the end, check if there was already a coefficient function 
  // for this region
  if( regionCoefs_.find( region ) != regionCoefs_.end() ) {
    if(allowReplacement == false){
      EXCEPTION( "Multiple definition of region with id " << region << "in a CoefFunctionMulti.");
    }
  }

  // adjust dependency of this coeffunction
  dependType_ = GetMaxCoefDependType(this->GetDependency(), coef->GetDependency());
  
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
  assert(dimType_ == CoefFunction::VECTOR);
  return rowSize_;
}

void CoefFunctionMulti::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  assert(dimType_ == CoefFunction::TENSOR);
  numRows = rowSize_;
  numCols = colSize_;

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

void  CoefFunctionMulti::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                  StdVector<Double >  & vals,
                                                  Grid* ptGrid,
                                                  const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                  bool updatedGeo ){
  EXCEPTION("CoefFunctionMulti::GetVectorValuesAtCoords: not implemented")

}

void  CoefFunctionMulti::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                  StdVector<Vector<Double> >  & vals,
                                                  Grid* ptGrid,
                                                  const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                  bool updatedGeo ){
  EXCEPTION("CoefFunctionMulti::GetScalarValuesAtPoints: not implemented")
}


void  CoefFunctionMulti::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                  StdVector<Complex >  & vals, 
                                                  Grid* ptGrid,
                                                  const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                  bool updatedGeo ){
  EXCEPTION("CoefFunctionMulti::GetVectorValuesAtCoords: not implemented")
}

void  CoefFunctionMulti::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                  StdVector<Vector<Complex> >  & vals,
                                                  Grid* ptGrid,
                                                  const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                  bool updatedGeo ){
  EXCEPTION("CoefFunctionMulti::GetVectorValuesAtCoords: not implemented")
}



} // end of namespace
