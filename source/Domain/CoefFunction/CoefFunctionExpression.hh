// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionExpression.hh
 *       \brief    {This coefficient function handles coefficients based on math parser
 *                 expressions. Therefore, the set Matrix/scalar/vector function expect
 *                 a string rather then a Double or Complex}
 *
 *       \date     Nov. 03, 2011
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONEXPRESSION_HH
#define COEFFUNCTIONEXPRESSION_HH

#include "CoefFunction.hh"

namespace CoupledField{

class MathParser;

//! Base class for real- and complex valued coefficient functions
template<class DATA_TYPE>
class CoefFunctionExpression : public CoefFunctionAnalytic {
};

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================

//! Coefficient function defined by real-valued mathematical expression
template<>
class CoefFunctionExpression<Double> : public CoefFunctionAnalytic,
                                       public enable_shared_from_this<CoefFunctionExpression<Double> >{
  
  public:
    CoefFunctionExpression(MathParser * mp);

    virtual ~CoefFunctionExpression();

    virtual string GetName() const { return "CoefFunctionExpression"; }


    //! \copydoc CoefFunction::GetComplexPart
    virtual PtrCoefFct GetComplexPart( Global::ComplexPart part );
    
    void GetTensor(Matrix<Double>& CoefMat, const LocPointMapped& lpm );

    void GetVector(Vector<Double>& CoefVec, const LocPointMapped& lpm );

    void GetScalar(Double& CoefScalar, const LocPointMapped& lpm );

    void SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols );

    void SetVector(const StdVector<std::string>& val);

    void SetScalar(const std::string& val);
    
    //! \copydoc CoefFunction::GetVecSize
    virtual UInt GetVecSize() const {
      assert(this->dimType_ == VECTOR );
      return coefVec_.GetSize();
    }
      
    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      assert(this->dimType_ == TENSOR );
      numRows = numRows_;
      numCols = numCols_;
    }


    std::string ToString() const;
    
    // =========================================================================
    // STRING REPRESENTATION 
    // =========================================================================
    //! \copydoc CoefFunctionAnalytic::GetStrScalar
    virtual void GetStrScalar( std::string& real, std::string& imag );

    //! \copydoc CoefFunctionAnalytic::GetStrVector
    virtual void GetStrVector( StdVector<std::string>& real, 
                               StdVector<std::string>& imag );

    //! \copydoc CoefFunctionAnalytic::GetStrTensor
    virtual void GetStrTensor( UInt& numRows, UInt& numCols,
                               StdVector<std::string>& real, 
                               StdVector<std::string>& imag );

    // COLLECTION ACCESS
    virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector< Double >  & vals,
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                          StdVector<shared_ptr<EntityList> >() );

    virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector<Vector< Double> >  & vals, 
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                          StdVector<shared_ptr<EntityList> >(),
                                          bool updatedGeo = false);
    
    virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector<Matrix<Double> >  & vals,
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                          StdVector<shared_ptr<EntityList> >() );



  protected:

    //! Coefficients for tensor
    StdVector<std::string > coefMat_;

    //! Number of rows of tensor
    UInt numRows_ = 0;
    
    //! Number of columns of tensor
    UInt numCols_ = 0;
    
    //! Coefficients for vector
    StdVector<std::string> coefVec_;

    //! Scalar coefficient
    std::string coefScalar_;

    //! Pointer to math parser instance
    MathParser* mp_ = NULL;
    
    //! Default coordinate system 
    CoordSystem* coordSysDefault_ = NULL;
    
    //! Handle for expression
    unsigned int mHandle_;
};

// ===========================================================================
//  COMPLEX VALUED COEFFICIENT FUNCTION
// ===========================================================================

//! Coefficient function defined by complex-valued mathematical expression
template<>
class CoefFunctionExpression<Complex> : public CoefFunctionAnalytic,
                                        public enable_shared_from_this<CoefFunctionExpression<Complex> >{
  
  public:
    CoefFunctionExpression(MathParser * mp);

    virtual ~CoefFunctionExpression();

    //! \copydoc CoefFunction::GetComplexPart
    virtual PtrCoefFct GetComplexPart( Global::ComplexPart part );
    
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
    
    //! \copydoc CoefFunction::GetVecSize
    virtual UInt GetVecSize() const {
      assert(this->dimType_ == VECTOR );
      return coefVecReal_.GetSize();
    }

    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      assert(this->dimType_ == TENSOR );
      numRows = numRows_;
      numCols = numCols_;
    }
    
    std::string ToString() const;
    
    
    // =========================================================================
    // STRING REPRESENTATION 
    // =========================================================================
    //! \copydoc CoefFunctionAnalytic::GetStrScalar
    virtual void GetStrScalar( std::string& real, std::string& imag );

    //! \copydoc CoefFunctionAnalytic::GetStrVector
    virtual void GetStrVector( StdVector<std::string>& real, 
                               StdVector<std::string>& imag );

    //! \copydoc CoefFunctionAnalytic::GetStrTensor
    virtual void GetStrTensor( UInt& numRows, UInt& numCols,
                               StdVector<std::string>& real, 
                               StdVector<std::string>& imag );


    // COLLECTION ACCESS
    virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector< Complex >  & vals,
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                          StdVector<shared_ptr<EntityList> >() );

    virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector< Complex >  & vals, 
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                              StdVector<shared_ptr<EntityList> >(),
                                              bool updatedGeo = false);

    virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector<Vector< Complex> >  & vals, 
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                              StdVector<shared_ptr<EntityList> >(),
                                              bool updatedGeo = false);

    virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                          StdVector<Matrix<Complex> >  & vals,
                                          Grid* ptGrid,
                                          const StdVector<shared_ptr<EntityList> >& srcEntities =
                                              StdVector<shared_ptr<EntityList> >()  );


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
    
    //! Default coordinate system 
    CoordSystem* coordSysDefault_;
    
    //@{
    //! Handle for expression
    unsigned int mHandleReal_;
    unsigned int mHandleImag_;
    //@}
};

}
#endif
