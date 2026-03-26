#ifndef COEFFUNCTIONTIMEFREQ_HH
#define COEFFUNCTIONTIMEFREQ_HH

#include <boost/signals2.hpp>
#include <memory>

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
  
  virtual string GetName() const { return "CoefFunctionTimeFreq"; }

};

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================
//! Real-valued coefficient function, depending on time / frequency
template<>
class CoefFunctionTimeFreq<Double> : public CoefFunctionAnalytic,
                                     public enable_shared_from_this<CoefFunctionTimeFreq<Double> > {
  public:

  //! Constructor
  CoefFunctionTimeFreq(MathParser * mp);

  //! Destructor
  virtual ~CoefFunctionTimeFreq();
  
  virtual string GetName() const { return "CoefFunctionTimeFreq<double>"; }


  //! \copydoc CoefFunction::GetComplexPart
  virtual PtrCoefFct GetComplexPart( Global::ComplexPart part );
  
  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, const LocPointMapped& lpm ) {
    assert(this->dimType_ == TENSOR);
    TransformTensorByCoordSys(coefMat, constCoefMat_, lpm);
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Double>& coefVec, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == VECTOR ||
           this->dimType_ == SCALAR );

    // in case of scalars, just set one entry in the vector
    if( this->dimType_ == SCALAR ) {
      coefVec.Resize(1);
      coefVec[0] = constCoefScalar_;
    } else {
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
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, const LocPointMapped& lpm ) {
    assert(this->dimType_ == SCALAR);
    coefScalar =  constCoefScalar_;
  }

  //! \copydoc CoefFunction::GetVecSize
  UInt GetVecSize() const {
    assert(this->dimType_ == VECTOR );
    return coefVec_.GetSize();
  }
    
  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    assert(this->dimType_ == TENSOR );
    numRows = constCoefMat_.GetNumRows();
    numCols = constCoefMat_.GetNumCols();
  }
  
  //! Set string tensor representation
  void SetTensor(const StdVector<std::string>& val, UInt nRows, UInt nCols );

  //! Set string vector representation
  void SetVector(const StdVector<std::string>& val);

  //! Set string scalar representation
  void SetScalar(const std::string& val);

  //! \copydoc CoefFunction::ToString
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
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Double >  & vals, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                            StdVector<shared_ptr<EntityList> >(),
                                            bool updatedGeo = false);

  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Vector<Double> >  & vals,
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

    //! Call-back method for re-calculation
    void Recalculate();

    //! Connection to math parser instance 
    boost::signals2::connection conn_;
    
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
    unsigned int mHandle_;
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
class CoefFunctionTimeFreq<Complex> : public CoefFunctionAnalytic,
                                      public enable_shared_from_this<CoefFunctionTimeFreq<Complex> > {
public:
  
  //! Constructor
  CoefFunctionTimeFreq(MathParser * mp);

  //! Destructor
  virtual ~CoefFunctionTimeFreq();

  virtual string GetName() const { return "CoefFunctionTimeFreq<complex>"; }

  //! \copydoc CoefFunction::GetComplexPart
  virtual PtrCoefFct GetComplexPart( Global::ComplexPart part );
  
  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<Complex>& coefMat, const LocPointMapped& lpm ) {
    assert(this->dimType_ == TENSOR);
    TransformTensorByCoordSys(coefMat, constCoefMat_, lpm);
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<Complex>& coefVec, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == VECTOR ||
           this->dimType_ == SCALAR );

    // in case of scalars, just set one entry in the vector
    if( this->dimType_ == SCALAR ) {
      coefVec.Resize(1);
      coefVec[0] = constCoefScalar_;
    } else {

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
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(Complex& coefScalar, const LocPointMapped& lpm ) {
    assert(this->dimType_ == SCALAR);
    coefScalar =  constCoefScalar_;
  }

  //! \copydoc CoefFunction::GetVecSize
  UInt GetVecSize() const {
    assert(this->dimType_ == VECTOR );
    return coefVecReal_.GetSize();
  }
    
  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    assert(this->dimType_ == TENSOR );
    numRows = constCoefMat_.GetNumRows();
    numCols = constCoefMat_.GetNumCols();
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

  //! \copydoc CoefFunction::ToString
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
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Complex >  & vals, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                            StdVector<shared_ptr<EntityList> >(),
                                            bool updatedGeo = false);

  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Vector<Complex> >  & vals,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                            StdVector<shared_ptr<EntityList> >(),
                                            bool updatedGeo = false);
  virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Matrix<Complex> >  & vals,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                            StdVector<shared_ptr<EntityList> >() );


protected:

  //! Call-back method for re-calculation
  void Recalculate();
  
  //! Connection to math parser instance (real) 
  boost::signals2::connection connReal_;

  //! Connection to math parser instance (imag) 
  boost::signals2::connection connImag_;

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
  unsigned int mHandleReal_;
  
  //! Handle for expression (imag)
  unsigned int mHandleImag_;
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
