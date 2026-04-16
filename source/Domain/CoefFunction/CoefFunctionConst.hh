#ifndef COEFFUNCTIONCONST_HH
#define COEFFUNCTIONCONST_HH

#include <memory>

#include "CoefFunction.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField {

//! Provide a coefficient function with constant expressions
template<typename T>
class CoefFunctionConst : public CoefFunctionAnalytic, 
                          public enable_shared_from_this<CoefFunctionConst<T> >
{
public:

  //! Constructor
  CoefFunctionConst() 
  : CoefFunctionAnalytic() {
    // this type of coefficient is always constant
    dependType_ = CONSTANT;
    isAnalytic_ = true;
    isComplex_ = std::is_same<T,Complex>::value;
    supportDerivative_ = true;
  }

  //! Destructor
  virtual ~CoefFunctionConst(){
    ;
  }
 
  virtual string GetName() const { return "CoefFunctionConst"; }


  //! \copydoc CoefFunction::GetComplexPart
  virtual PtrCoefFct GetComplexPart( Global::ComplexPart part );

   /** This is the const special case which needs no lpm as long as we have no coordSys */
  Matrix<T>& GetTensor()
  {
     assert(this->dimType_ == TENSOR);
     return constCoefMat_;
   }

  /** Again a CoefFunctionConst special case */
  void GetTensor(Matrix<T>& coefMat) {
    coefMat = GetTensor();
  }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm) {
    // we are the const version
    Matrix<T> locMatrix;
    GetTensor(locMatrix);
    TransformTensorByCoordSys(coefMat, locMatrix, lpm);
  }

  const Matrix<T>& GetTensor() const
  {
    assert(dimType_ == TENSOR);
    return constCoefMat_;
  }

  //! \copydoc CoefFunction::GetVector
  void GetVector(Vector<T>& coefVec, const LocPointMapped& lpm )
  {
    assert(this->dimType_ == VECTOR || this->dimType_ == SCALAR );

    // in case of scalars, just set one entry in the vector
    if( this->dimType_ == SCALAR ) {
      coefVec.Resize(1);
      coefVec[0] = coefScalar_;
    } else {
      // if no coordinate system is set, just
      // use internal vector
      if( !coordSys_ ) {
        coefVec = coefVec_;
      } else {
        // otherwise, perform local -> global mapping
        Vector<Double> pointCoord;

        // in case we are interested in nodal points or have no shapeMap for other reasons
        // lpm.GetGlobal(pointCoord) would be an option - if one resolves the const issue
        lpm.shapeMap->Local2Global(pointCoord,lpm.lp);


        // Afterwards rotate the local vector back to global coordinates
        this->coordSys_->Local2GlobalVector( coefVec, coefVec_, pointCoord );
      }
    }
  }

  /** special CoefFunctionConst variant without lpm which is not used anyway!
   * Does not work with coordSys set!  */
  void GetVector(Vector<T>& coefVec)
  {
    assert(this->dimType_ == VECTOR || this->dimType_ == SCALAR );

    if(coordSys_)
      EXCEPTION("Call simplief CoefFunctionConst::GetVector() only without coordSys!");

    // in case of scalars, just set one entry in the vector. Therefore no Vector return!
    if( this->dimType_ == SCALAR ) {
      coefVec.Resize(1);
      coefVec[0] = coefScalar_; // TODO: why not check coordSys for scalar?!
    }
    else
      coefVec = coefVec_;
  }

  const Vector<T>& GetVector() const {
      assert(this->dimType_ == VECTOR);
      return coefVec_;
    }


  /** special CoefFunctionConst variant without lpm which is not used anyway! */
  T GetScalar() const {
    assert(this->dimType_ == SCALAR);
    return coefScalar_;
  }


  /** special CoefFunctionConst variant without lpm which is not used anyway! */
  void GetScalar(T& coefScalar) {
    coefScalar = GetScalar();
  }

  //! \copydoc CoefFunction::GetScalar
  void GetScalar(T& coefScalar, const LocPointMapped& lpm ) {
    GetScalar(coefScalar);
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
  
  //! Set associated coordinate system
  virtual void SetCoordinateSystem(CoordSystem* cSys) {
    if( cSys->GetName() != "default")
      coordSys_ = cSys;
  }

  //! Set tensor value
  void SetTensor(const Matrix<T>& CoefMat){
    assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
    this->dimType_ = TENSOR;
    this->constCoefMat_ = CoefMat;
  }

  //! Set vector values
  void SetVector(const Vector<T>& val){
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

  //! \copydoc CoefFunction::IsZero
  bool IsZero() const;
  
  //! \copydoc CoefFunction::ToString
  std::string ToString() const {
    switch( dimType_ ) {
      case NO_DIM:
        return "";
        break;
      case SCALAR:
        return boost::lexical_cast<std::string>(coefScalar_);
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

  //! \copydoc CoefFunction::SetDerivativeOperation
  virtual void SetDerivativeOperation(CoefDerivativeType type){
    this->derivType_ = type;
    //derivative of a constant is always zero...
    switch(dimType_){
    case SCALAR:
      this->coefScalar_ = 0.0;
      break;
    case VECTOR:
      this->coefVec_.Init();
      break;
    case TENSOR:
      this->constCoefMat_.Init();
      break;
    default:
      break;
    }
    return;
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


  // COLLECTION ACCESS
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<T >  & vals, Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );

  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Vector<T> >  & vals,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >(),
                                        bool updatedGeo = false);

  virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                        StdVector<Matrix<T> >  & vals,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );

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
