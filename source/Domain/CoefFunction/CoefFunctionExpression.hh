// =====================================================================================
// 
//       Filename:  coefFunctionExpression.hh
// 
//    Description:  This coefficient function handles coefficients based on math parser
//                  expressions. Therefore, the set Matrix/scalar/vector function expect
//                  a string rather then a Double or Complex
// 
//        Version:  1.0
//        Created:  11/03/2011 12:17:37 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================


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
    CoefFunctionExpression();

    ~CoefFunctionExpression();

    void GetTensor(Matrix<Double>& CoefMat, const LocPointMapped& lpm );

    void GetVector(Vector<Double>& CoefVec, const LocPointMapped& lpm );

    void GetScalar(Double& CoefScalar, const LocPointMapped& lpm );

    void SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols );

    void SetVector(const StdVector<std::string>& val);

    void SetScalar(const std::string& val);
    
    std::string ToString() const;

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
    CoefFunctionExpression();

    ~CoefFunctionExpression();

    void GetTensor(Matrix<Complex>& CoefMat, const LocPointMapped& lpm );

    void GetVector(Vector<Complex>& CoefVec, const LocPointMapped& lpm );

    void GetScalar(Complex& CoefScalar, const LocPointMapped& lpm );

    void SetTensor(const StdVector<std::string>& realVal, 
                   const StdVector<std::string>& imagVal,
                   UInt nRows, UInt nCols );

    void SetVector(const StdVector<std::string>& realVal,
                   const StdVector<std::string>& imagVal );

    void SetScalar(const std::string& realVal,
                   const std::string& imagVal);
    
    std::string ToString() const;

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
