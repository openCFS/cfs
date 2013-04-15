#include "CoefFunctionExpression.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"
namespace CoupledField{


// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================
CoefFunctionExpression<Double>::CoefFunctionExpression(MathParser * mp) :
        CoefFunctionAnalytic(),
        mp_(mp),
        mHandle_(mp_->GetNewHandle(true))
        {
  
  // this type of coefficient is always variable
  dependType_ = GENERAL;
  
  // this coefficient function is still analytic, as it
  // is defined by analytical expressions.
  isAnalytic_ = true;
  
  isComplex_ = false;
  
  // do not use rotated coordsys initially
  this->coordSys_ = domain->GetCoordSystem();
}   

CoefFunctionExpression<Double>::~CoefFunctionExpression(){
  if( domain ) {
    mp_->ReleaseHandle(mHandle_);
  }
}

PtrCoefFct CoefFunctionExpression<Double>::GetComplexPart( Global::ComplexPart part ) {
  if( part != Global::REAL ) {
    EXCEPTION( "Can only return the real-part of a real-valued coefFunction" );
  }
  return shared_from_this();

}

void CoefFunctionExpression<Double>::GetTensor( Matrix<Double>& coefMat, 
                                                const LocPointMapped& lpm ){
  assert(this->dimType_ == CoefFunction::TENSOR);
  Vector<Double> pointCoord;
  Matrix<Double> locMatrix;

  // First, obtain global coordinates of current point, register it at the mathParser
  // and compute local coefficient vector
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), pointCoord);
  this->mp_->EvalMatrix(mHandle_, locMatrix, 
                        this->numRows_, this->numCols_ ); 

  // Rotate material, if coordinate system is not the global one
  if( coordSys_->GetName() != "default") {
    // Obtain rotation matrix
    Matrix<Double> rotMatrix;
    coordSys_->GetFullGlobRotationMatrix( rotMatrix, pointCoord );

    EXCEPTION("The rotation is not fully finished ':-(\n" << 
              "Here we have to add a call to the method BaseMaterial::PerformRotation "
              "This method should be moved to the base class of the CoefFunction"
              "In addition the initial rotation of the material must be incorporated"
              "somewehre in string-notation, as we are generally dealing with string"
              "parameters."
              "Thus we should treat the case, where rotation angles are multiples of "
              "90 degree separately, where the entries are just interchanged");
  } else {
    coefMat = locMatrix;
  }
}

void CoefFunctionExpression<Double>::GetVector( Vector<Double>& coefVec, 
                                                const LocPointMapped& lpm ){
  assert(this->dimType_ == CoefFunction::VECTOR ||
         this->dimType_ == CoefFunction::SCALAR);
  Vector<Double> pointCoord, locVector;

  // in case of scalars, just set one entry in the vector
  if( this->dimType_ == SCALAR ) {
    coefVec.Resize(1);
    GetScalar(coefVec[0], lpm);
  } else {
    // First, obtain global coordinates of current point, register it at the mathParser
    // and compute local coefficient vector
    lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
    this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), pointCoord);
    this->mp_->EvalVector(mHandle_, locVector );

    // Afterwards rotate the local vector back to global coordinates
    this->coordSys_->Local2GlobalVector( coefVec, locVector, pointCoord );
  }
}

void CoefFunctionExpression<Double>::GetScalar(Double& coefScalar, 
                                               const LocPointMapped& lpm){
  assert(this->dimType_ == CoefFunction::SCALAR);
  // First, obtain global coordinates of current point and  register it at the mathParser
  Vector<Double> pointCoord;;
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), pointCoord);
  coefScalar = this->mp_->Eval(mHandle_);
}

void CoefFunctionExpression<Double>::
SetTensor( const StdVector<std::string>& val, 
           UInt nRows, UInt nCols ){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  this->coefMat_ = val;
  this->dimType_ = CoefFunction::TENSOR;
  
  // security check: make sure that nRows x nCols = val.size()
  if( nCols * nRows != val.GetSize() ) {
    EXCEPTION("Tensor is supposed to have size " << nRows << " x " 
              << nCols << ", but only " << val.GetSize() 
              << " entries have been provided!");
  }
  this->numRows_ = nRows;
  this->numCols_ = nCols;
  
  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandle_, val.Serialize(','));
}

void CoefFunctionExpression<Double>::
SetVector( const StdVector<std::string>& val ){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  this->coefVec_ = val;
  this->dimType_ = CoefFunction::VECTOR;
  
  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandle_, val.Serialize(','));
}

void CoefFunctionExpression<Double>::SetScalar(const std::string& val){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
  this->coefScalar_ = val;
  this->dimType_ = CoefFunction::SCALAR;
  
  // Register expression at math parser
  mp_->SetExpr(mHandle_, val);
}

std::string CoefFunctionExpression<Double>::ToString() const {
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
      break;
  }
  return "";
}

void CoefFunctionExpression<Double>::
GetStrScalar( std::string& real, std::string& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
  real = coefScalar_;
  imag = "0.0";
  
}

void CoefFunctionExpression<Double>::
GetStrVector( StdVector<std::string>& real, 
              StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  real = coefVec_;
  imag = "0.0";

}

void CoefFunctionExpression<Double>::
GetStrTensor( UInt& numRows, UInt& numCols,
              StdVector<std::string>& real, 
              StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  numRows = numRows_;
  numCols = numCols_;
  real = coefMat_;
  imag.Resize(real.GetSize());
  imag.Init("0.0");
}

// ===========================================================================
//  COMPLEX VALUED COEFFICIENT FUNCTION
// ===========================================================================
CoefFunctionExpression<Complex>::CoefFunctionExpression(MathParser * mp) :
        CoefFunctionAnalytic(),
        mp_(mp),
        mHandleReal_(mp_->GetNewHandle(true)),
        mHandleImag_(mp_->GetNewHandle(true)){
  
  // this type of coefficient is always variable
  dependType_ = GENERAL;

  // this coefficient function is still analytic, as it
  // is defined by analytical expressions.
  isAnalytic_ = true;
  
  isComplex_ = true;

  // always store default coordinate system
  this->coordSys_ = domain->GetCoordSystem();
}   

CoefFunctionExpression<Complex>::~CoefFunctionExpression(){
  if( domain ) {
    mp_->ReleaseHandle(mHandleReal_);
    mp_->ReleaseHandle(mHandleImag_);
  }
}

PtrCoefFct CoefFunctionExpression<Complex>::GetComplexPart( Global::ComplexPart part ) {
  PtrCoefFct ret;
  if( part  == Global::COMPLEX ) {
    ret = shared_from_this();
  } else if ( part == Global::REAL || part == Global::IMAG ) {
    shared_ptr<CoefFunctionExpression<Double> > real(new CoefFunctionExpression<Double>(mp_));
    switch(dimType_) {
      case SCALAR:
        if( part == Global::REAL )
          real->SetScalar( coefScalarReal_ );
        else
          real->SetScalar( coefScalarImag_ );  
        break;
      case VECTOR:
        if( part == Global::REAL )
          real->SetVector( coefVecReal_ );
        else
          real->SetVector( coefVecImag_ );
        break;
      case TENSOR:
        if( part == Global::REAL )
          real->SetTensor( coefMatReal_, numRows_, numCols_ );
        else
          real->SetTensor( coefMatImag_, numRows_, numCols_ );
        break;
      default:
        EXCEPTION( "Unknown dimType" );
        break;
    }
    ret = real;
  }
  return ret;
}
    
void CoefFunctionExpression<Complex>::GetTensor( Matrix<Complex>& coefMat, 
                                                 const LocPointMapped& lpm ){
  assert(this->dimType_ == CoefFunction::TENSOR);
  Vector<Double> pointCoord;
  Matrix<Complex> locMatrix(this->numRows_, this->numCols_);
  Matrix<Double> temp;

  // First, obtain global coordinates of current point, register it at the mathParser
  // and compute local coefficient vector
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), pointCoord);
  this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), pointCoord);
  this->mp_->EvalMatrix(mHandleReal_, temp, 
                        this->numRows_, this->numCols_ );
  locMatrix.SetPart(Global::REAL, temp);
  this->mp_->EvalMatrix(mHandleImag_, temp, 
                        this->numRows_, this->numCols_ );
  locMatrix.SetPart(Global::IMAG, temp);
  
  // Rotate material, if coordinate system is not the global one
  if( coordSys_->GetName() != "default") {
    // Obtain rotation matrix
    Matrix<Double> rotMatrix;
    coordSys_->GetFullGlobRotationMatrix( rotMatrix, pointCoord );

    EXCEPTION("The rotation is not fully finished ':-(\n" << 
              "Here we have to add a call to the method BaseMaterial::PerformRotation "
              "This method should be moved to the base class of the CoefFunction"
              "In addition the initial rotation of the material must be incorporated"
              "somewehre in string-notation, as we are generally dealing with string"
              "parameters."
              "Thus we should treat the case, where rotation angles are multiples of "
              "90 degree separately, where the entries are just interchanged");
  } else {
    coefMat = locMatrix;
  }
}

void CoefFunctionExpression<Complex>::GetVector( Vector<Complex>& coefVec, 
                                                const LocPointMapped& lpm ){
  assert(this->dimType_ == CoefFunction::VECTOR);
  Vector<Double> temp, pointCoord;
  Vector<Complex> locVector;

  if( this->dimType_ == SCALAR ) {
     coefVec.Resize(1);
     GetScalar(coefVec[0], lpm);
   } else {

     // First, obtain global coordinates of current point, register it at the mathParser
     // and compute local coefficient vector
     lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
     this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), pointCoord);
     this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), pointCoord);
     this->mp_->EvalVector(mHandleReal_, temp );
     locVector.Resize(temp.GetSize());
     locVector.SetPart(Global::REAL, temp);
     this->mp_->EvalVector(mHandleImag_, temp );
     locVector.SetPart(Global::IMAG, temp);

     // Afterwards rotate the local vector back to global coordinates
     this->coordSys_->Local2GlobalVector( coefVec, locVector, pointCoord );
   }
}

void CoefFunctionExpression<Complex>::GetScalar( Complex& coefScalar, 
                                                 const LocPointMapped& lpm){
  Double real, imag;
  assert(this->dimType_ == CoefFunction::SCALAR);
  // First, obtain global coordinates of current point and  register it at the mathParser
  Vector<Double> pointCoord;;
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), pointCoord);
  this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), pointCoord);
  
  real = this->mp_->Eval(mHandleReal_);
  imag = this->mp_->Eval(mHandleImag_);
  coefScalar = Complex(real, imag);
}

void CoefFunctionExpression<Complex>::SetTensor( const StdVector<std::string>& realVal, 
                                                 const StdVector<std::string>& imagVal,
                                                 UInt nRows, UInt nCols ){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  this->coefMatReal_ = realVal;
  this->coefMatImag_ = imagVal;
  this->dimType_ = CoefFunction::TENSOR;
  
  // security check: make sure that nRows x nCols = val.size()
  if( nCols * nRows != realVal.GetSize() ) {
    EXCEPTION("Tensor is supposed to have size " << nRows << " x " 
              << nCols << ", but only " << realVal.GetSize() 
              << " entries have been provided!");
  }
  this->numRows_ = nRows;
  this->numCols_ = nCols;
  
  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandleReal_, realVal.Serialize(','));
  mp_->SetExpr(mHandleImag_, imagVal.Serialize(','));
}

void CoefFunctionExpression<Complex>::SetVector( const StdVector<std::string>& realVal,
                                                 const StdVector<std::string>& imagVal ){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  this->coefVecReal_ = realVal;
  this->coefVecImag_ = imagVal;
  this->dimType_ = CoefFunction::VECTOR;
  
  // Set expression at math parser, where each entry is separated by ", "
  mp_->SetExpr(mHandleReal_, realVal.Serialize(','));
  mp_->SetExpr(mHandleImag_, imagVal.Serialize(','));
}

void CoefFunctionExpression<Complex>::SetScalar( const std::string& realVal,
                                                 const std::string& imagVal ){
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
  this->coefScalarReal_ = realVal;
  this->coefScalarImag_ = imagVal;
  this->dimType_ = CoefFunction::SCALAR;
  
  // Register expression at math parser
  mp_->SetExpr(mHandleReal_, realVal);
  mp_->SetExpr(mHandleImag_, imagVal);
}

std::string CoefFunctionExpression<Complex>::ToString() const {
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

void CoefFunctionExpression<Complex>::
GetStrScalar( std::string& real, std::string& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == SCALAR) );
   real = coefScalarReal_;
   imag = coefScalarImag_;
}

void CoefFunctionExpression<Complex>::
GetStrVector( StdVector<std::string>& real, 
              StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == VECTOR) );
  real = coefVecReal_;
  imag = coefVecImag_;
}

void CoefFunctionExpression<Complex>::
GetStrTensor( UInt& numRows, UInt& numCols,
                            StdVector<std::string>& real, 
                            StdVector<std::string>& imag ) {
  assert((this->dimType_ == NO_DIM) || (this->dimType_ == TENSOR) );
  numRows = numRows_;
  numCols = numCols_;
  real = coefMatReal_;
  imag = coefMatImag_;
}

void CoefFunctionExpression<Double>::GetScalarValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<Double >  & vals){
  assert(this->dimType_ == CoefFunction::SCALAR);

  vals.Resize(points.GetSize());
  vals.Init();
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), points[curPoint]);
    vals[curPoint] = this->mp_->Eval(mHandle_);
  }
}

void CoefFunctionExpression<Double>::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                      StdVector<Vector<Double> >  & vals){
  Vector<Double> locVector;
  if( this->dimType_ == SCALAR ) {
    StdVector<Double> tmpVals;
    this->GetScalarValuesAtCoords(points,tmpVals);
    vals.Resize(tmpVals.GetSize(), Vector<Double>(1));
    for(UInt i=0;i<tmpVals.GetSize();i++){
      vals[i][0] = tmpVals[i];
    }
  }else{
    vals.Resize(points.GetSize());
    vals.Init();
    for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
      this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), points[curPoint]);
      this->mp_->EvalVector(mHandle_, locVector );
      this->coordSys_->Local2GlobalVector( vals[curPoint], locVector, points[curPoint] );
    }
  }

}

void CoefFunctionExpression<Double>::GetTensorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                       StdVector<Matrix<Double> >  & vals){

  assert(this->dimType_ == CoefFunction::TENSOR);

  Matrix<Double> locMatrix;
  vals.Resize(points.GetSize(),Matrix<Double>(this->numRows_,this->numCols_));
  vals.Init();
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    this->mp_->SetCoordinates(mHandle_, *(this->coordSys_), points[curPoint]);
    this->mp_->EvalMatrix(mHandle_, locMatrix,
                          this->numRows_, this->numCols_ );

    if( coordSys_->GetName() != "default") {
      // Obtain rotation matrix
      Matrix<Double> rotMatrix;
      coordSys_->GetFullGlobRotationMatrix( rotMatrix, points[curPoint] );

      EXCEPTION("The rotation is not fully finished ':-(\n" <<
                "Here we have to add a call to the method BaseMaterial::PerformRotation "
                "This method should be moved to the base class of the CoefFunction"
                "In addition the initial rotation of the material must be incorporated"
                "somewehre in string-notation, as we are generally dealing with string"
                "parameters."
                "Thus we should treat the case, where rotation angles are multiples of "
                "90 degree separately, where the entries are just interchanged");
    } else {
      vals[curPoint] = locMatrix;
    }
  }
}


void CoefFunctionExpression<Complex>::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<Complex >  & vals){
  Double real, imag;
  assert(this->dimType_ == CoefFunction::SCALAR);
  // First, obtain global coordinates of current point and  register it at the mathParser
  Vector<Double> pointCoord;;
  vals.Resize(points.GetSize());
  vals.Init();
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), points[curPoint]);
    this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), points[curPoint]);
    real = this->mp_->Eval(mHandleReal_);
    imag = this->mp_->Eval(mHandleImag_);
    vals[curPoint] = Complex(real, imag);
  }
}

void CoefFunctionExpression<Complex>::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                                      StdVector<Vector<Complex> >  & vals){
  assert(this->dimType_ == CoefFunction::VECTOR);
  Vector<Double> temp;
  Vector<Complex> locVector;

  if( this->dimType_ == SCALAR ) {
     StdVector<Complex> tmpVals;
     this->GetVectorValuesAtCoords(points,tmpVals);
     vals.Resize(tmpVals.GetSize(), Vector<Complex>(1));
     for(UInt i=0;i<tmpVals.GetSize();i++){
       vals[i][0] = tmpVals[i];
     }
   } else {
     vals.Resize(points.GetSize());
     vals.Init();
     for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
       this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), points[curPoint]);
       this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), points[curPoint]);
       this->mp_->EvalVector(mHandleReal_, temp );
       locVector.Resize(temp.GetSize());
       locVector.SetPart(Global::REAL, temp);
       this->mp_->EvalVector(mHandleImag_, temp );
       locVector.SetPart(Global::IMAG, temp);

       // Afterwards rotate the local vector back to global coordinates
       this->coordSys_->Local2GlobalVector( vals[curPoint], locVector, points[curPoint] );
     }
   }
}

void CoefFunctionExpression<Complex>::GetVectorValuesAtCoords( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<Matrix<Complex> >  & vals){

  assert(this->dimType_ == CoefFunction::TENSOR);
  Matrix<Complex> locMatrix(this->numRows_, this->numCols_);
  Matrix<Double> temp;

  vals.Resize(points.GetSize(),Matrix<Complex>(this->numRows_,this->numCols_));
  vals.Init();
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    this->mp_->SetCoordinates(mHandleReal_, *(this->coordSys_), points[curPoint]);
    this->mp_->SetCoordinates(mHandleImag_, *(this->coordSys_), points[curPoint]);
    this->mp_->EvalMatrix(mHandleReal_, temp,
                          this->numRows_, this->numCols_ );
    locMatrix.SetPart(Global::REAL, temp);
    this->mp_->EvalMatrix(mHandleImag_, temp,
                          this->numRows_, this->numCols_ );
    locMatrix.SetPart(Global::IMAG, temp);


    // Rotate material, if coordinate system is not the global one
    if( coordSys_->GetName() != "default") {
      // Obtain rotation matrix
      Matrix<Double> rotMatrix;
      coordSys_->GetFullGlobRotationMatrix( rotMatrix, points[curPoint] );

      EXCEPTION("The rotation is not fully finished ':-(\n" <<
                "Here we have to add a call to the method BaseMaterial::PerformRotation "
                "This method should be moved to the base class of the CoefFunction"
                "In addition the initial rotation of the material must be incorporated"
                "somewehre in string-notation, as we are generally dealing with string"
                "parameters."
                "Thus we should treat the case, where rotation angles are multiples of "
                "90 degree separately, where the entries are just interchanged");
    } else {
      vals[curPoint] = locMatrix;
    }
  }
}

}
