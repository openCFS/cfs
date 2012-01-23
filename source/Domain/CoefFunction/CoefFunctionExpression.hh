#ifndef COEFFUNCTIONEXPRESSION_HH
#define COEFFUNCTIONEXPRESSION_HH

#include "CoefFunction.hh"

#include "Utils/mathParser/mathParser.hh"

namespace CoupledField{

//! Base class for real- and complex valued coefficient functions
template<class DATA_TYPE>
class CoefFunctionExpression : public CoefFunction{
};

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================

//! Coefficient function defined by real-valued mathematical expression
template<>
class CoefFunctionExpression<Double> : public CoefFunction {

public:

  //! Constructor
  CoefFunctionExpression();

  virtual ~CoefFunctionExpression();

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& CoefMat, const LocPointMapped& lpm );

  //! \see CoefFunction::GetVector
  void GetVector(Vector<Double>& CoefVec, const LocPointMapped& lpm );

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& CoefScalar, const LocPointMapped& lpm );

  //! Set string tensor representation
  void SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols );

  //! Set string vector representation
  void SetVector(const StdVector<std::string>& val);

  //! Set scalar vector representation
  void SetScalar(const std::string& val);

  //! \see CoefFunction::IsComplex
  bool IsComplex(){ 
    return false; 
  }
  
  //! \see CoefFunction::ToString
  std::string ToString() const;

  //! Return MathParser handle
  MathParser::HandleType GetHandle() { 
    return mHandle_;
  }

protected:

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
};

// ===========================================================================
//  COMPLEX VALUED COEFFICIENT FUNCTION
// ===========================================================================

//! Coefficient function defined by complex-valued mathematical expression
template<>
class CoefFunctionExpression<Complex> : public CoefFunction {

public:

  //! Constructor
  CoefFunctionExpression();

  //! Destructor
  virtual ~CoefFunctionExpression();

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& CoefMat, const LocPointMapped& lpm );

  //! \see CoefFunction::GetVector
  void GetVector(Vector<Complex>& CoefVec, const LocPointMapped& lpm );

  //! \see CoefFunction::GetScalar
  void GetScalar(Complex& CoefScalar, const LocPointMapped& lpm );

  //! Set string tensor representation
  void SetTensor(const StdVector<std::string>& realVal, 
                 const StdVector<std::string>& imagVal,
                 UInt nRows, UInt nCols );

  //! Set string vector representation
  void SetVector(const StdVector<std::string>& realVal,
                 const StdVector<std::string>& imagVal );

  //! Set scalar vector representation
  void SetScalar(const std::string& realVal,
                 const std::string& imagVal);

  //! \see CoefFunction::IsComplex
  bool IsComplex(){ 
    return true; 
  }
  
  //! \see CoefFunction::ToString
  std::string ToString() const;

  //! Return MathParser handle for real-part
  MathParser::HandleType GetHandleReal() { 
    return mHandleReal_;
  }

  //! Return MathParser handle for imag-part
  MathParser::HandleType GetHandleImag() { 
    return mHandleImag_;
  }

protected:

  //@{
  //! Coefficients for tensor
  StdVector<std::string > coefMatReal_;
  StdVector<std::string > coefMatImag_;
  //@}

  //! Number of rows of tensor
  UInt numRows_;

  //! Number of columns of tensor
  UInt numCols_;

  //@{
  //! Coefficients for vector
  StdVector<std::string> coefVecReal_;
  StdVector<std::string> coefVecImag_;
  //@}

  //@{
  //! Scalar coefficient
  std::string coefScalarReal_;
  std::string coefScalarImag_;
  //@}

  //! Pointer to math parser instance
  MathParser* mp_;

  //@{
  //! Handle for expression
  MathParser::HandleType mHandleReal_;
  MathParser::HandleType mHandleImag_;
  //@}
};

}
#endif
