#include "CoefFunctionTimeFreq.hh"
#include "boost/bind.hpp" 
#include "boost/lexical_cast.hpp"

namespace CoupledField{


// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================

CoefFunctionTimeFreq<Double>::
CoefFunctionTimeFreq() {
  
  dependType_ = TIMEFREQ;
  
  // obtain handle from internal variable coefficient function
  mp_ = domain->GetMathParser();
  mHandle_ = mp_->GetNewHandle(true);
  
  // register callback mechanism if expression changes
  mp_->AddExpChangeCallBack(
      boost::bind(&CoefFunctionTimeFreq<Double>::Recalculate, this ),
      mHandle_ );
}

 CoefFunctionTimeFreq<Double>::
~CoefFunctionTimeFreq() {
  
}
 
 void CoefFunctionTimeFreq<Double>::
 SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols ) {

   assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );

   // ensure that expression really depends only on time / freq
   for( UInt i = 0; i < val.GetSize(); ++i ) {
     if( ExprDependsOnSpace(val[i]) ) {
       EXCEPTION("Expression '" << val.ToString() 
                 << "' depends on spatial coordinates. "
                 << "Refusing to create CoefFunctionTimeFreq." );
     }
   }

   // set value at coef and recalculate  
   // security check: make sure that nRows x nCols = val.size()
   if( nCols * nRows != val.GetSize() ) {
     EXCEPTION("Tensor is supposed to have size " << nRows << " x " 
               << nCols << ", but only " << val.GetSize() 
               << " entries have been provided!");
   }
   
   this->coefMat_ = val;
   this->dimType_ = CoefFunction::TENSOR;
   this->numRows_ = nRows;
   this->numCols_ = nCols;

   // Set expression at math parser, where each entry is separated by ", "
   mp_->SetExpr(mHandle_, val.Serialize(','));

   // Important: recalculate 
   Recalculate();
}

void CoefFunctionTimeFreq<Double>::
SetVector(const StdVector<std::string>& val) {
  
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  
  // ensure that expression really depends only on time / freq
  for( UInt i = 0; i < val.GetSize(); ++i ) {
    if( ExprDependsOnSpace(val[i]) ) {
      EXCEPTION("Expression '" << val.ToString() 
                << "' depends on spatial coordinates. "
                << "Refusing to create CoefFunctionTimeFreq." );
    }
  }

  this->coefVec_ = val;
  this->dimType_ = CoefFunction::VECTOR;

  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandle_, val.Serialize(','));

  // Important: recalculate 
  Recalculate();
}

void CoefFunctionTimeFreq<Double>::
SetScalar(const std::string& val) {

  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );

  // ensure that expression really depends only on time / freq
  if( ExprDependsOnSpace(val) ) {
    EXCEPTION("Expression '" << val << "' depends on spatial coordinates. "
              << "Refusing to create CoefFunctionTimeFreq." );
  }

  // set value at coef and recalculate  
  this->coefScalar_ = val;
  this->dimType_ = CoefFunction::SCALAR;

  // Register expression at math parser
  mp_->SetExpr(mHandle_, val);
  
  // Important: recalculate 
  Recalculate();
}

std::string  CoefFunctionTimeFreq<Double>::ToString() const {
  switch( this->dimType_) {
    case NO_DIM:
      return "";
      break;
    case SCALAR:
      return lexical_cast<std::string>(this->coefScalar_);
      break;
    case VECTOR:
      return this->coefVec_.ToString();
      break;
    case TENSOR:
      WARN("Not clean implemented");
      return this->coefMat_.ToString();
      break;
    default:
      EXCEPTION("Missing case");
  }
}

void CoefFunctionTimeFreq<Double>::
Recalculate() {

  if( dimType_ == TENSOR )  {
    this->mp_->EvalMatrix(mHandle_, constCoefMat_, 
                          this->numRows_, this->numCols_ ); 
  } else if( dimType_ == VECTOR ) {
    mp_->EvalVector(mHandle_, constCoefVec_ );
  } else if( dimType_ == SCALAR ) {
    constCoefScalar_ = this->mp_->Eval(mHandle_);
  } else {
    EXCEPTION("Not handled case")
  }
}  

// ===========================================================================
//  COMPLEX VALUED COEFFICIENT FUNCTION
// ===========================================================================
CoefFunctionTimeFreq<Complex>::
CoefFunctionTimeFreq() {
  
  dependType_ = TIMEFREQ;
  
  // obtain handle from internal variable coefficient function
  mp_ = domain->GetMathParser();
  mHandleReal_ = mp_->GetNewHandle(true);
  mHandleImag_ = mp_->GetNewHandle(true);
  
  // register callback mechanism if expression changes
  mp_->AddExpChangeCallBack(
      boost::bind(&CoefFunctionTimeFreq<Complex>::Recalculate, this ),
      mHandleReal_ );

  mp_->AddExpChangeCallBack(
      boost::bind(&CoefFunctionTimeFreq<Complex>::Recalculate, this ),
      mHandleImag_ );
}

CoefFunctionTimeFreq<Complex>::
~CoefFunctionTimeFreq() {
  
}

void CoefFunctionTimeFreq<Complex>::
SetTensor(const StdVector<std::string>& realVal, 
          const StdVector<std::string>& imagVal,
          UInt nRows, UInt nCols ) {

  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );

  // ensure that expression really depends only on time / freq
   for( UInt i = 0; i < realVal.GetSize(); ++i ) {
     if( ExprDependsOnSpace(realVal[i]) ) {
       EXCEPTION("Expression '" << realVal.ToString() 
                 << "' depends on spatial coordinates. "
                 << "Refusing to create CoefFunctionTimeFreq." );
     }
   }
   for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
     if( ExprDependsOnSpace(imagVal[i]) ) {
       EXCEPTION("Expression '" << imagVal.ToString() 
                 << "' depends on spatial coordinates. "
                 << "Refusing to create CoefFunctionTimeFreq." );
     }
   }
   
   // set value at coef and recalculate  
   // security check: make sure that nRows x nCols = val.size()
   if( nCols * nRows != realVal.GetSize() ) {
     EXCEPTION("Tensor is supposed to have size " << nRows << " x " 
               << nCols << ", but only " << realVal.GetSize() 
               << " entries have been provided!");
   }
   
   if( nCols * nRows != imagVal.GetSize() ) {
     EXCEPTION("Tensor is supposed to have size " << nRows << " x " 
               << nCols << ", but only " << imagVal.GetSize() 
               << " entries have been provided!");
   }
  
   this->coefMatReal_ = realVal;
   this->coefMatImag_ = imagVal;
   this->dimType_ = CoefFunction::TENSOR;

   this->numRows_ = nRows;
   this->numCols_ = nCols;
   // Set expression at math parser, where each entry is separated by ", "
   mp_->SetExpr(mHandleReal_, realVal.Serialize(','));
   mp_->SetExpr(mHandleImag_, imagVal.Serialize(','));

  // Important: recalculate 
  Recalculate();

}

void CoefFunctionTimeFreq<Complex>::
SetVector(const StdVector<std::string>& realVal,
                 const StdVector<std::string>& imagVal ) {
  
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  
  // ensure that expression really depends only on time / freq
  for( UInt i = 0; i < realVal.GetSize(); ++i ) {
    if( ExprDependsOnSpace(realVal[i]) ) {
      EXCEPTION("Expression '" << realVal.ToString() 
                << "' depends on spatial coordinates. "
                << "Refusing to create CoefFunctionTimeFreq." );
    }
  }
  for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
    if( ExprDependsOnSpace(imagVal[i]) ) {
      EXCEPTION("Expression '" << imagVal.ToString() 
                << "' depends on spatial coordinates. "
                << "Refusing to create CoefFunctionTimeFreq." );
    }
  }

  this->coefVecReal_ = realVal;
  this->coefVecImag_ = imagVal;
  this->dimType_ = CoefFunction::VECTOR;

  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandleReal_, realVal.Serialize(','));
  mp_->SetExpr(mHandleImag_, imagVal.Serialize(','));

  // Important: recalculate 
  Recalculate();
}

void CoefFunctionTimeFreq<Complex>::
SetScalar(const std::string& realVal,
          const std::string& imagVal) {

  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );

  // ensure that expression really depends only on time / freq
  if( ExprDependsOnSpace(realVal)  || 
      ExprDependsOnSpace(imagVal)) {
    EXCEPTION("Expression '(" << realVal  << ", " << imagVal 
              << ")' depends on spatial coordinates. "
              << "Refusing to create CoefFunctionTimeFreq." );
  }

  // Register expression at math parser
  mp_->SetExpr(mHandleReal_, realVal);
  mp_->SetExpr(mHandleImag_, imagVal);

  // Important: recalculate 
  Recalculate();
}

std::string CoefFunctionTimeFreq<Complex>::ToString() const {
  std::string ret;
  switch( this->dimType_) {
    case NO_DIM:
      return "";
      break;
    case SCALAR:
      ret = std::string("(") + this->coefScalarReal_;
      ret + ", ";
      ret += this->coefScalarImag_;
      ret +" )";
      return ret;
      break;
    case VECTOR:
      ret  = "Real-Part: " + this->coefVecReal_.ToString() + "\n";
      ret += "Imag-Part: " + this->coefVecImag_.ToString();
      return ret;
      break;
    case TENSOR:
      WARN("Not clean implemented");
      ret  = "Real-Part: " + this->coefMatReal_.ToString() + "\n";
      ret += "Imag-Part: " + this->coefMatImag_.ToString();
      return ret;
      break;
    default:
      EXCEPTION("Missing case");
  }
  return "";
}
void CoefFunctionTimeFreq<Complex>::
Recalculate() {
  
  
  std::cerr << "Re-calculate!\n";
  if( dimType_ == TENSOR )  {
    Matrix<Double> temp;
    constCoefMat_.Resize(this->numRows_, this->numCols_);
    this->mp_->EvalMatrix(mHandleReal_, temp, 
                          this->numRows_, this->numCols_ );
    constCoefMat_.SetPart( Global::REAL, temp );
    this->mp_->EvalMatrix(mHandleImag_, temp, 
                          this->numRows_, this->numCols_ );
    constCoefMat_.SetPart( Global::IMAG, temp );
    
  } else if( dimType_ == VECTOR ) {
    
    Vector<Double> temp;
    this->mp_->EvalVector(mHandleReal_, temp );
    constCoefVec_.Resize(temp.GetSize());
    constCoefVec_.SetPart(Global::REAL, temp);
    this->mp_->EvalVector(mHandleImag_, temp );
    constCoefVec_.SetPart(Global::IMAG, temp);
    
  } else if( dimType_ == SCALAR ) {
    Double real, imag;
    real = this->mp_->Eval(mHandleReal_);
    imag = this->mp_->Eval(mHandleImag_);
    constCoefScalar_ = Complex(real, imag);
  } else {
    EXCEPTION("Not handled case")
  }
}
} // namespace
