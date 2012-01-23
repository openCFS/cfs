#ifndef COEFFUNCTIONTIMEFREQ_HH
#define COEFFUNCTIONTIMEFREQ_HH

#include "CoefFunction.hh"
#include "CoefFunctionExpression.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField{

//! Coefficient function, depending on time / frequency

//! This coefficient function is just depending on the current time /
//! frequency, but not on space. Thus we can use an efficient evaluation scheme
//! which is based on a callback function. Every instance of the CoefFunction
//! registers itself with the MathParser class and if the time / frequency
//! changes, the coefficients are re-calculated. Otherwise we just return
//! the internally stored data. 
//! To parse the string representation, we call 
template<typename T>
class CoefFunctionTimeFreq : public CoefFunction {
  
};

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================
//! Real-valued coefficient function, depending on time / frequency
template<>
class CoefFunctionTimeFreq<Double> : public CoefFunction {
  public:
  //! Constructor
  CoefFunctionTimeFreq();
  
  //! Destructor
  virtual ~CoefFunctionTimeFreq();
  
  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, const LocPointMapped& lpm ) {
    assert(this->dimType_ == TENSOR);
    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      coefMat =  constCoefMat_;
    } else {
      EXCEPTION(
          "The rotation is not fully finished ':-(\n" << 
          "Here we have to add a call to the method BaseMaterial::PerformRotation "
          "This method should be moved to the base class of the CoefFunction"
          "In addition the initial rotation of the material must be incorporated"
          "somewhere in string-notation, as we are generally dealing with string"
          "parameters."
          "Thus we should treat the case, where rotation angles are multiples of "
          "90 degree separately, where the entries are just interchanged");
    }
  }

  //! \see CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == VECTOR);

    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      coefVec = constCoefVec_;

    } else {
      // otherwise, perform local -> global mapping
      Vector<Double> pointCoord;
      lpm.shapeMap->Local2Global(pointCoord,lpm.lp);

      // Afterwards rotate the local vector back to global coordinates
      this->coordSys_->Local2GlobalVector( coefVec, 
                                           constCoefVec_, pointCoord );
    }
  }

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, const LocPointMapped& lpm ) {
    assert(this->dimType_ == SCALAR);
    coefScalar =  constCoefScalar_;
  }

  //! Set string tensor representation
  void SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols );

  //! Set string vector representation
  void SetVector(const StdVector<std::string>& val);

  //! Set string scalar representation
  void SetScalar(const std::string& val);

  //! \see CoefFunction::IsComplex
  bool IsComplex(){ 
    return false; 
  }

  //! \see CoefFunction::ToString
  std::string ToString() const;

  protected:

    //! Call-back method for re-calculation
    void Recalculate();

    // =====================================
    //  STRING PARAMETER REPRESENTATION
    // =====================================
    //@{
    //! \name String Parameter Representation
    
    //! Coefficients for tensor
    StdVector<std::string > coefMat_;

    //! Number of rows of tensor
    UInt numRows_;

    //! Number of columns of tensor
    UInt numCols_;

    //! Coefficients for vector
    StdVector<std::string> coefVec_;

    //! Scalar coefficient
    std::string coefScalar_;

    //! Pointer to math parser instance
    MathParser* mp_;

    //! Handle for expression
    MathParser::HandleType mHandle_;
    //@}
    
    
    // =====================================
    //  CONSTANT PARAMETER REPRESENTATION
    // =====================================
    //@{
    //! \name Constant Parameter Representation

    //! Constant coefficient tensor
    Matrix<Double> constCoefMat_;

    //! Constant coefficient vector
    Vector<Double> constCoefVec_;

    //! Constant scalar value
    Double constCoefScalar_;
    //@}
};

// ===========================================================================
//  COMPLEX VALUED COEFFICIENT FUNCTION
// ===========================================================================
//! Complex-valued coefficient function, depending on time / frequency
template<>
class CoefFunctionTimeFreq<Complex> : public CoefFunction {
public:
  //! Constructor
  CoefFunctionTimeFreq();

  //! Destructor
  virtual ~CoefFunctionTimeFreq();

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat, const LocPointMapped& lpm ) {
    assert(this->dimType_ == TENSOR);
    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      coefMat =  constCoefMat_;
    } else {
      EXCEPTION(
          "The rotation is not fully finished ':-(\n" << 
          "Here we have to add a call to the method BaseMaterial::PerformRotation "
          "This method should be moved to the base class of the CoefFunction"
          "In addition the initial rotation of the material must be incorporated"
          "somewhere in string-notation, as we are generally dealing with string"
          "parameters."
          "Thus we should treat the case, where rotation angles are multiples of "
          "90 degree separately, where the entries are just interchanged");
    }
  }

  //! \see CoefFunction::GetVector
  void GetVector(Vector<Complex>& coefVec, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == VECTOR);

    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      coefVec = constCoefVec_;
    } else {
      // otherwise, perform local -> global mapping
      Vector<Double> pointCoord;
      lpm.shapeMap->Local2Global(pointCoord,lpm.lp);

      // Afterwards rotate the local vector back to global coordinates
      this->coordSys_->Local2GlobalVector( coefVec, constCoefVec_,
                                           pointCoord );
    }
  }

  //! \see CoefFunction::GetScalar
  void GetScalar(Complex& coefScalar, const LocPointMapped& lpm ) {
    assert(this->dimType_ == SCALAR);
    coefScalar =  constCoefScalar_;
  }

  //! Set string tensor representation
  void SetTensor(const StdVector<std::string>& realVal, 
                 const StdVector<std::string>& imagVal,
                 UInt nRows, UInt nCols );

  //! Set string vector representation
  void SetVector(const StdVector<std::string>& realVal,
                 const StdVector<std::string>& imagVal );

  //! Set string scalar representation
  void SetScalar(const std::string& realVal,
                 const std::string& imagVal);

  //! \see CoefFunction::IsComplex
  bool IsComplex(){ 
    return true; 
  }

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:

  //! Call-back method for re-calculation
  void Recalculate();

  // =====================================
  //  STRING PARAMETER REPRESENTATION
  // =====================================
  //@{
  //! \name String Parameter Representation
  
  //! Coefficients for tensor (real )
  StdVector<std::string > coefMatReal_;
  
  //! Coefficients for tensor (imag )
  StdVector<std::string > coefMatImag_;

  //! Number of rows of tensor
  UInt numRows_;

  //! Number of columns of tensor
  UInt numCols_;

  //! Coefficients for vector (real)
  StdVector<std::string> coefVecReal_;
  
  //! Coefficients for vector (imag)
  StdVector<std::string> coefVecImag_;

  //! Scalar coefficient (real)
  std::string coefScalarReal_;
  
  //! Scalar coefficient (imag)
  std::string coefScalarImag_;

  //! Pointer to math parser instance
  MathParser* mp_;

  //! Handle for expression (real)
  MathParser::HandleType mHandleReal_;
  
  //! Handle for expression (imag)
  MathParser::HandleType mHandleImag_;
  //@}
  
  // =====================================
  //  CONSTANT PARAMETER REPRESENTATION
  // =====================================
  //@{
  //! \name Constant Parameter Representation

  //! Constant coefficient tensor
  Matrix<Complex> constCoefMat_;

  //! Constant coefficient vector
  Vector<Complex> constCoefVec_;

  //! Constant scalar value
  Complex constCoefScalar_;
   //@}
};

}
#endif
