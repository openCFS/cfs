#ifndef COEFFUNCTIONCONST_HH
#define COEFFUNCTIONCONST_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField {

//! Provide a coefficient function with constant expressions
template<typename T>
class CoefFunctionConst : public CoefFunctionAnalytic{
public:

  //! Constructor
  CoefFunctionConst() 
  : CoefFunctionAnalytic() {
    // this type of coefficient is always constant
    dependType_ = CONST;
    isAnalytic_ = true;
  }

  //! Destructor
  virtual ~CoefFunctionConst(){
    ;
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<T>& coefMat, 
                 const LocPointMapped& lpm ) {
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

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<T>& coefVec, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == VECTOR);

    // if no coordinate system is set, just
    // use internal vector
    if( !coordSys_ ) {
      coefVec = coefVec_;

    } else {
      // otherwise, perform local -> global mapping
      Vector<Double> pointCoord;
      lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
      
      // Afterwards rotate the local vector back to global coordinates
      this->coordSys_->Local2GlobalVector( coefVec, coefVec_, pointCoord );
    }
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(T& coefScalar, 
                 const LocPointMapped& lpm ) {
    assert(this->dimType_ == SCALAR);
    coefScalar =  coefScalar_;
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

  //! Set tensor value
  void SetTensor(Matrix<T>& CoefMat){
    assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
    this->dimType_ = TENSOR;
    this->constCoefMat_ = CoefMat;
  }

  //! Set vector values
  void SetVector(Vector<T> val){
    assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
    this->coefVec_ = val;
    this->dimType_ = VECTOR;
  }
  //! Set scalar value
  void SetScalar(T val){
    assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
    this->coefScalar_ = val;
    this->dimType_ = SCALAR;
  }

  //! \copydoc CoefFunction::IsComplex
  bool IsComplex(){
    return std::tr1::is_same<T,Complex>::value;
  }
  
  //! \copydoc CoefFunction::ToString
  std::string ToString() const {
    switch( dimType_ ) {
      case NO_DIM:
        return "";
        break;
      case SCALAR:
        return lexical_cast<std::string>(coefScalar_);
        break;
      case VECTOR:
        return coefVec_.ToString();
        break;
      case TENSOR:
        return constCoefMat_.ToString();
        break;
      default:
        EXCEPTION("Missing case");
        break;
    }
    return "";
  }

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
protected:
  
  //! Constant coefficient tensor
  Matrix<T> constCoefMat_;

  //! Constant coefficient vector
  Vector<T> coefVec_;

  //! Constant scalar value
  T coefScalar_;

};

}
#endif
