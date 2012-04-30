#include "CoefFunctionConst.hh"

namespace CoupledField {
template<>
PtrCoefFct CoefFunctionConst<Double>::GetComplexPart( Global::ComplexPart part ) {
  // the only meaningful value here is part == Global::REAL, 
  if( part != Global::REAL ) {
    EXCEPTION( "Can only return the real-part of a real-valued coefFunction" );
  }
  return shared_from_this();
}

template<>
PtrCoefFct CoefFunctionConst<Complex>::GetComplexPart( Global::ComplexPart part ) {
  PtrCoefFct ret;
  if( part  == Global::COMPLEX ) {
    ret = shared_from_this();
  } else if ( part == Global::REAL || part == Global::IMAG ) {
    shared_ptr<CoefFunctionConst<Double> > real(new CoefFunctionConst<Double>());
    switch(dimType_) {
      case SCALAR:
        real->SetScalar( coefScalar_.real() );
        break;
      case VECTOR:
        real->SetVector( coefVec_.GetPart( part ) );
        break;
      case TENSOR:
        real->SetTensor( constCoefMat_.GetPart( part ) );
        break;
      default:
        EXCEPTION( "Unknown dimType" );
        break;
    }
    ret = real;
  }
  return ret;
}

template<>
void CoefFunctionConst<Double>::GetStrScalar( std::string& real, std::string& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
  real = lexical_cast<std::string>(coefScalar_);
  imag = "0.0";
}

template<>
void CoefFunctionConst<Complex>::GetStrScalar( std::string& real, std::string& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
  real = lexical_cast<std::string>(coefScalar_.real() );
  imag = lexical_cast<std::string>(coefScalar_.imag() );
}

template<>
void CoefFunctionConst<Double>::GetStrVector( StdVector<std::string>& real, 
                                              StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  UInt numEntries = coefVec_.GetSize();
  real.Resize( numEntries );
  imag.Resize( numEntries );
  imag.Init("0.0");
  for( UInt i = 0; i < numEntries; ++i ) {
    real[i] = lexical_cast<std::string>( coefVec_[i] );
  }
}

template<>
void CoefFunctionConst<Complex>::GetStrVector( StdVector<std::string>& real, 
                                               StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  UInt numEntries = coefVec_.GetSize();
   real.Resize( numEntries );
   imag.Resize( numEntries );
   for( UInt i = 0; i < numEntries; ++i ) {
     real[i] = lexical_cast<std::string>( coefVec_[i].real() );
     imag[i] = lexical_cast<std::string>( coefVec_[i].imag() );
   }
}

template<>
void CoefFunctionConst<Double>::GetStrTensor( UInt& numRows, UInt& numCols,
                                              StdVector<std::string>& real, 
                                              StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  numRows = constCoefMat_.GetNumRows();
  numCols = constCoefMat_.GetNumCols();
  UInt numEntries = numRows * numCols;
  real.Resize( numEntries );
  imag.Resize( numEntries );
  
  imag.Init("0.0");
  UInt pos = 0;
  for( UInt i = 0; i < numRows; ++i ) {
    for( UInt j = 0; j < numCols; ++j ) {
      real[pos++] = lexical_cast<std::string>( constCoefMat_[i][j] );
    }
  }
}

template<>
void CoefFunctionConst<Complex>::GetStrTensor( UInt& numRows, UInt& numCols,
                                               StdVector<std::string>& real, 
                                               StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  numRows = constCoefMat_.GetNumRows();
  numCols = constCoefMat_.GetNumCols();
  UInt numEntries = numRows * numCols;
  real.Resize( numEntries );
  imag.Resize( numEntries );

  UInt pos = 0;
  for( UInt i = 0; i < numRows; ++i ) {
    for( UInt j = 0; j < numCols; ++j ) {
      real[pos] = lexical_cast<std::string>( constCoefMat_[i][j].real() );
      imag[pos++] = lexical_cast<std::string>( constCoefMat_[i][j].imag() );
    }
  }
}

}
