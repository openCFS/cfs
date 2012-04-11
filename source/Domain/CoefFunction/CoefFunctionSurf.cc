#include "CoefFunctionSurf.hh"

#include <boost/tr1/type_traits.hpp>
namespace CoupledField {


template<typename T> CoefFunctionSurf<T>::
CoefFunctionSurf( bool mapNormal ) 
: CoefFunction() {


  // not sure about the following one
  dependType_ = GENERAL;
  mapNormal_ = mapNormal;
  
}


template<typename T> void
CoefFunctionSurf<T>::SetVolumeCoefs( std::map<RegionIdType, PtrCoefFct> coefs ) {
  // ensure that at least one entry is present
    assert( coefs.size() != 0);

    std::map<RegionIdType, PtrCoefFct>::iterator it = coefs.begin();

    // loop over all entries and set following variables
    dimType_ = it->second->GetDimType();
    regions_.insert(it->first);
    it++;
    for( ; it != coefs.end(); ++it ) {
      if( it->second->GetDimType() != dimType_ ) {
        EXCEPTION( "CoefFunctions for surface CoefFunction have " <<
                    "inconsistent dimension type" );
      }
      regions_.insert(it->first);
    }
    
    coefs_ = coefs;
    
    // adjust dimensionality: In case of normal mapping,
    // we perform a vector multiplication.
    if( mapNormal_ ) {
      switch(dimType_) {
        case SCALAR:
          dimType_ = VECTOR;
          break;
        case VECTOR:
          dimType_ = SCALAR;
          break;
        case TENSOR:
          dimType_ = VECTOR;
          break;
        default:
          EXCEPTION( "Unknown dimensionality" ); 
          break;     
      }
    }
}

template<typename T>
CoefFunctionSurf<T>::~CoefFunctionSurf() {

}


template<typename T> void  
CoefFunctionSurf<T>::GetTensor(Matrix<T>& coefMat, 
                               const LocPointMapped& lpm ){
  assert(this->dimType_ == TENSOR);
}

template<typename T> void 
CoefFunctionSurf<T>::GetVector(Vector<T>& coefVec, 
                               const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR);
  EXCEPTION("not implemented");

}

template<typename T> void
CoefFunctionSurf<T>::GetScalar(T& coefScalar, 
                               const LocPointMapped& lpm ){
  assert(this->dimType_ == SCALAR);

  // create local point for surface
  LocPointMapped surfLpm(lpm);
  surfLpm.SetSurfInfo( regions_);
  RegionIdType region = surfLpm.lpmVol->ptEl->regionId;
  if( mapNormal_ ) {
    Vector<T> coefVec;
    coefs_[region]->GetVector(coefVec, *surfLpm.lpmVol );
    coefScalar = coefVec * surfLpm.normal;
  } else {
    coefs_[region]->GetScalar(coefScalar, *surfLpm.lpmVol );
  }
}

template<typename T> UInt
CoefFunctionSurf<T>::GetVecSize() const {
  return 0;
}

template<typename T> void
CoefFunctionSurf<T>::GetTensorSize( UInt& numRows, UInt& numCols ) const {

}

template<typename T> bool
CoefFunctionSurf<T>::IsComplex() {
  return std::tr1::is_same<T,Complex>::value; 
}


template<typename T> std::string 
CoefFunctionSurf<T>::ToString() const {
  return "";
}

template class CoefFunctionSurf<Double>;
template class CoefFunctionSurf<Complex>;
}
