#include "BaseMaterial.hh"

#include "Utils/StdVector.hh"
#include "MatVec/Matrix.hh"

#include "Materials/Models/Preisach.hh"
//#include "Materials/Models/VectorPreisach.hh"
//#include "Materials/Models/VectorPreisachv7.hh"
#include "Materials/Models/VectorPreisachMayergoyz.hh"
#include "Materials/Models/VectorPreisachSutor.hh"
#include "Materials/Models/SimplePreisachInv.hh"
#include "Materials/Models/PiezoMicroModelHF.hh"
#include "Materials/Models/PiezoMicroModelBK.hh"

#include "Domain/Domain.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionApprox.hh"

#include "Domain/CoefFunction/CoefXpr.hh"

#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Utils/helperStructs.hh"

using std::string;
using std::map;
using std::set;


namespace CoupledField
{

  BaseMaterial::MatDescriptorNl::MatDescriptorNl() {
    measAccuracy = 0.0;
    maxVal = 0.0;
    angle = 0.0;
    temperature = 0.0;
    zScaling = 0.0;
    factor = 0.0;
    approxData = NULL;
    approxType = NO_APPROX_TYPE;
  }
  
  BaseMaterial::MatDescriptorNl::~MatDescriptorNl() {
    if (approxData) {
      delete approxData;
    }
  }



  // ***********************
  //   Default Constructor
  // ***********************
  BaseMaterial::BaseMaterial( MaterialClass matClass,
                              MathParser* mp,
                              CoordSystem * defaultCoosy )
  : class_(matClass)
  {
    name_ = "";
    model_ = MAT_MODEL_LINEAR;
    mp_ = mp;
    coosy_ = defaultCoosy;

    hyst_  = NULL;
    isHysteresis_  = false;
    isHystInverse_ = false;
    computeHystInverse_ = false;
    dim_ = 0;
    dimVecHyst_ = 0;
    XpreviousVEC_ = NULL;
    YpreviousVEC_ = NULL;
    vecHyst_ = NULL;
    hystY_ = NULL;
    hystZ_ = NULL;

    piezoMicroModel_   = NULL;
    isPiezoMicroModel_ = false;
  }

  BaseMaterial::~BaseMaterial() {
    mp_ = NULL;
    coosy_ = NULL;
  }


  // set the name of the material set
  void BaseMaterial::SetName(const std::string &name) {
    name_ = name;
  }

  // returns the name of the material
  std::string BaseMaterial::GetName() const {
    return name_;
  }

  // returns the class of the material
  MaterialClass BaseMaterial::GetClass() const {
    return class_;
  }

  // set the symmetry type
  void BaseMaterial::SetSymmetryType(MaterialType matType, SymmetryType symType) {
    symmetryType_[matType] = symType;
  }

  // get the symmetry type
  BaseMaterial::SymmetryType BaseMaterial::GetSymmetryType(MaterialType matType) const {
    SymmetryType ret = GENERAL;
    std::map<MaterialType, SymmetryType>::const_iterator it;
    it = symmetryType_.find(matType);
    if( it == symmetryType_.end() ) {
      EXCEPTION("Could not find symmetry type for material '"
          << MaterialTypeEnum.ToString(matType) << "'");
    } else {
      ret = it->second;
    }

    return ret;
  }

  // Returns the material model selected for this material
  BaseMaterial::MaterialModel BaseMaterial::GetModel() const {
    return model_;
  }

  // Set the material model
  void BaseMaterial::SetModel(MaterialModel model) {
    if (model == MAT_MODEL_LINEAR) {
      model_ = model;
    }
    else if (class_ == MECHANIC && model == MAT_MODEL_LINEARVISCOELASTIC) {
      model_ = model;
    }
    else if (class_ == SMOOTH && model == MAT_MODEL_LINEARVISCOELASTIC) {
      model_ = model;
    }
    else {
      EXCEPTION("Cannot use " << MaterialModelEnum.ToString(model)
          << " model for material '" << name_ /*<< "', because it is a "
          << MaterialClassEnum.ToString(class_) << " material."*/);
    }
  }


  void BaseMaterial::SetString(const std::string& param, MaterialType matType) {
    //check, if allowed
    if ( isAllowed_.find( matType ) == isAllowed_.end() ) {
      matTypeNotAllowed( matType, "string" );
    }
    else {
      stringParams_[matType] = param;
      isSet_.insert( matType );
    }
  }
  
  void BaseMaterial::GetString( std::string& param, MaterialType matType) const {
    StringMap::const_iterator pos = stringParams_.find( matType );

    if ( pos == stringParams_.end() ) {
      matTypeNotInDataBase( matType, "string" );
    }
    else {
      param = pos->second;
    }
  }    

  void BaseMaterial::SetScalar(Integer param, MaterialType matType) {
    //check, if allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
      matTypeNotAllowed( matType, "scalar" );
    }
    else {
      isSet_.insert( matType );
      integerParams_[matType] = param;
    }
  }
  
  // Set a scalar real constant material parameter
  void BaseMaterial::SetScalar(Double param, MaterialType matType,
                               Global::ComplexPart dataType )
  {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "scalar");
    }
    else if (dataType == Global::REAL) {
      shared_ptr< CoefFunctionConst<Double> > fct(new CoefFunctionConst<Double>());
      fct->SetScalar(param);
      scalarCoef_[matType] = fct;
      isSet_.insert(matType);
    }
    else {
      dataTypeNotAllowed4SetGet(dataType, "SetScalar-Double");
    }
  }

  // set a scalar complex constant material parameter
  void BaseMaterial::SetScalar(Complex param, MaterialType matType,
                               Global::ComplexPart dataType )
  {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "scalar");
    }
    else if (dataType == Global::INTEGER) {
      dataTypeNotAllowed4SetGet(dataType, "SetScalar-Double");
    }
    else {
      shared_ptr< CoefFunctionConst<Complex> > fct(new CoefFunctionConst<Complex>());
      if (dataType == Global::COMPLEX) {
        fct->SetScalar(param);
      }
      else if (dataType == Global::REAL) {
        fct->SetScalar(Complex(param.real(), 0.0));
      }
      else {
        fct->SetScalar(Complex(0.0, param.imag()));
      }
      scalarCoef_[matType] = fct;
      isSet_.insert(matType);
    }
  }

  // set a real constant material vector
  void BaseMaterial::SetVector(const Vector<Double>& param, MaterialType matType,
                               Global::ComplexPart dataType) {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "vector");
    }
    else if (dataType == Global::REAL) {
      shared_ptr< CoefFunctionConst<Double> > fct(new CoefFunctionConst<Double>());
      fct->SetVector(param);
      vectorCoef_[matType] = fct;
      isSet_.insert(matType);
    }
    else {
      dataTypeNotAllowed4SetGet(dataType, "SetVector-Double");
    }
  }

  // set a complex constant material vector
  void BaseMaterial::SetVector(const Vector<Complex>& param, MaterialType matType,
                               Global::ComplexPart dataType)
  {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "vector");
    }
    else if (dataType == Global::INTEGER) {
      dataTypeNotAllowed4SetGet(dataType, "SetVector-Double");
    }
    else {
      shared_ptr< CoefFunctionConst<Complex> > fct(new CoefFunctionConst<Complex>());
      if (dataType == Global::COMPLEX) {
        fct->SetVector(param);
      }
      else{
        Vector<Complex> v(param.GetSize());
        v.SetPart(dataType, param.GetPart(dataType), true);
        fct->SetVector(v);
      }
      vectorCoef_[matType] = fct;
      isSet_.insert(matType);
    }
  }

  // set a real constant material tensor
  void BaseMaterial::SetTensor(const Matrix<Double>& param, MaterialType matType,
                               Global::ComplexPart dataType)
  {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "tensor");
    }
    else if (dataType == Global::REAL) {
      shared_ptr< CoefFunctionConst<Double> > fct(new CoefFunctionConst<Double>());
      fct->SetTensor(param);
      tensorCoef_[matType] = fct;
      tensorOrigCoef_[matType] = fct;
      isSet_.insert(matType);
    }
    else {
      dataTypeNotAllowed4SetGet(dataType, "SetTensor-Double");
    }
  }

  // set a complex constant material tensor
  void BaseMaterial::SetTensor(const Matrix<Complex>& param, MaterialType matType,
                               Global::ComplexPart dataType )
  {
    if (isAllowed_.find(matType) == isAllowed_.end()) {
      matTypeNotAllowed(matType, "tensor");
    }
    else if (dataType == Global::INTEGER) {
      dataTypeNotAllowed4SetGet(dataType, "SetTensor-Double");
    }
    else {
      shared_ptr< CoefFunctionConst<Complex> > fct(new CoefFunctionConst<Complex>());
      if (dataType == Global::COMPLEX) {
        fct->SetTensor(param);
      }
      else{
        Matrix<Complex> t(param.GetNumRows(), param.GetNumCols());
        t.SetPart(dataType, param.GetPart(dataType), true);
        fct->SetTensor(t);
      }
      tensorCoef_[matType] = fct;
      tensorOrigCoef_[matType] = fct;
      isSet_.insert(matType);
    }
  }

  void BaseMaterial::GetScalar(Integer& param, MaterialType matType) const {
    IntegerMap::const_iterator pos = integerParams_.find( matType );

    if ( pos == integerParams_.end() ) {
      matTypeNotInDataBase( matType, "scalar" );
    }
    else {
      param = pos->second;
    }
  }

  void BaseMaterial::GetScalar( Double& param, MaterialType matType,
                                Global::ComplexPart matDataType) const
  {
    PtrCoefFct scal = GetScalCoefFnc( matType, matDataType );
    if( scal->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    LocPointMapped lpm;
    scal->GetScalar(param, lpm);
  }

  void BaseMaterial::GetScalar( Complex& param, MaterialType matType,
                                Global::ComplexPart matDataType) const
  {
    PtrCoefFct scal = GetScalCoefFnc( matType, matDataType);
    if( scal->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    LocPointMapped lpm;
    scal->GetScalar(param, lpm);
  }

  // Get a constant vector
  void BaseMaterial::GetVector( Vector<Double>& param, MaterialType matType,
                                Global::ComplexPart matDataType ) const
  {
    PtrCoefFct vec = GetVectorCoefFnc( matType, matDataType );
    if (vec->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    LocPointMapped lpm;
    vec->GetVector(param, lpm);
  }

  void BaseMaterial::GetVector( Vector<Complex>& param,  MaterialType matType,
                                Global::ComplexPart matDataType ) const
  {
    PtrCoefFct vec = GetVectorCoefFnc( matType, matDataType );
    if (vec->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    LocPointMapped lpm;
    vec->GetVector(param, lpm);
  }

  // get a real material tensor
  void BaseMaterial::GetTensor( Matrix<Double>& param, MaterialType matType,
                                Global::ComplexPart dataType,
                                SubTensorType subTensor ) const
  {
    PtrCoefFct tens = GetTensorCoefFnc(matType, subTensor, dataType);
    if (tens->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    if ( dataType == Global::REAL || dataType == Global::IMAG) {
      Matrix<Complex> matTensor;
      LocPointMapped lpm;
      tens->GetTensor(matTensor, lpm);
      param = matTensor.GetPart( dataType );
    }
    else {
      dataTypeNotAllowed4SetGet( dataType, "GetTensor-Double" );
    }
  }

  // get a complex material tensor
  void BaseMaterial::GetTensor( Matrix<Complex>& param, MaterialType matType,
                                Global::ComplexPart dataType,
                                SubTensorType subTensor ) const
  {
    PtrCoefFct tens = GetTensorCoefFnc(matType, subTensor, dataType);
    if (tens->GetDependency() != CoefFunction::CONSTANT) {
      matDataNotConstant(matType);
    }

    Matrix<Complex> matTensor;
    LocPointMapped lpm;
    tens->GetTensor(matTensor, lpm);

    if ( dataType == Global::COMPLEX ) {
      param = matTensor;
    }
    if ( dataType == Global::REAL || dataType == Global::IMAG) {
      Matrix<Double> help;
      help = matTensor.GetPart( dataType );
      param.Resize( matTensor.GetNumRows(), matTensor.GetNumCols() );
      param.Init();
      param.SetPart( dataType, help );
    }
    else {
      dataTypeNotAllowed4SetGet( dataType, "GetTensor-Complex" );
    }
  }

  void BaseMaterial::SetNonLinMatIso( MaterialType matType, MatDescriptorNl& data ) {
    
    if( nonlinIsoParams_.find(matType) != nonlinIsoParams_.end() ) {
      EXCEPTION( "Nonlinear material parameter '" << MaterialTypeEnum.ToString(matType)
                 << "' was already set");
    }
    nonlinIsoParams_[matType] = data;
  }

  void BaseMaterial::SetNonLinMatAniso( MaterialType matType, StdVector<MatDescriptorNl>& data ) {
    if( nonlinAnisoParams_.find(matType) != nonlinAnisoParams_.end() ) {
      EXCEPTION( "Nonlinear material parameter '" << MaterialTypeEnum.ToString(matType)
                 << "' was already set");
    }
    nonlinAnisoParams_[matType] = data;
  }

  void BaseMaterial::SetNonLinMatIsoTempDependBH( MaterialType matType, StdVector<MatDescriptorNl>& data ) {
    if( nonlinIsoTempDependBHParams_.find(matType) != nonlinIsoTempDependBHParams_.end() ) {
      EXCEPTION( "Nonlinear material parameter '" << MaterialTypeEnum.ToString(matType)
                 << "' was already set");
    }
    nonlinIsoTempDependBHParams_[matType] = data;
  }

  void BaseMaterial::SetAnhystMagModel( const std::string name ){
    anhystereticModel_ = name;
  }
  
  void BaseMaterial::SetCoefFct( MaterialType matType, PtrCoefFct coef ) {
    
    // check, if material type is allowed
    if (  isAllowed_.find( matType ) == isAllowed_.end() ) {
       std::string dim = CoefFunction::coefDimType.ToString(coef->GetDimType());
       matTypeNotAllowed( matType, dim );
     }
    
    // switch depending on dimType of coefficient function
    switch(coef->GetDimType()) {
      case CoefFunction::SCALAR:
        scalarCoef_[matType] = coef;
        break;
      case CoefFunction::VECTOR:
        vectorCoef_[matType] = coef;
        break;
      case CoefFunction::TENSOR:
        tensorCoef_[matType] = coef;
        tensorOrigCoef_[matType] = coef;
        break;
      default:
        EXCEPTION("Unknown entry type of coefficient function");
        break;
    }

    isSet_.insert( matType ); // enables printing to info.xml
  }

  bool BaseMaterial::IsSet( MaterialType matType ) const {

    StringMap::const_iterator stringIt = stringParams_.find( matType );
    IntegerMap::const_iterator intIt = integerParams_.find( matType );
    CoefMap::const_iterator coefScalIt = scalarCoef_.find( matType );
    CoefMap::const_iterator coefVecIt = vectorCoef_.find( matType );
    CoefMap::const_iterator coefTensIt = tensorCoef_.find( matType );
        
    return ( stringIt != stringParams_.end() ||
             intIt != integerParams_.end()   ||
             coefScalIt != scalarCoef_.end() ||
             coefVecIt != vectorCoef_.end()  ||
             coefTensIt != tensorCoef_.end() );
  }

  void  BaseMaterial::matTypeNotAllowed(MaterialType matType, const string& dim ) const {

    string help = MaterialTypeEnum.ToString( matType );
    EXCEPTION( "Material type (" <<  dim <<  ") " << help << 
               " is not available for " << MaterialClassEnum.ToString(class_)
               << " Database" );
  }
  

  void  BaseMaterial::dataTypeNotAllowed4SetGet(Global::ComplexPart dataType, 
						                                    const string& msg ) const
  {
    string help;
    help = Global::complexPart.ToString( dataType );
    EXCEPTION( "Datatype " << help << " is not allowed in function " << msg );
  }


  void BaseMaterial::dataTypeNotAllowed(Global::ComplexPart dataType, MaterialType matType ) const {

    string help1, help2;
    help1 = Global::complexPart.ToString( dataType );
    help2 = MaterialTypeEnum.ToString( matType );
    EXCEPTION( "Datatype " << help1 << " is not allowed for material type " 
               << help2 << " in material data base "
               << MaterialClassEnum.ToString(class_) );
  }

  void BaseMaterial::matTypeNotInDataBase(MaterialType matType, const string& dim ) const {

    string help = MaterialTypeEnum.ToString( matType );

    EXCEPTION( "Material type (" << dim << ") " << help << " of material "
               << name_ << " was not read from/defined in material file" );
  }


  void  BaseMaterial::setMakesNoSense(Global::ComplexPart dataType, const string& msg ) const {
    EXCEPTION( "Set of " << msg << " makes no sense with datatype "
               << Global::complexPart.ToString( dataType ) );
  }


  void BaseMaterial::subTensorNotAvailable(MaterialType matType, SubTensorType subTensor )
  {
    string help1, help2;
    help1 = MaterialTypeEnum.ToString( matType );
    Enum2String(subTensor, help2);
    EXCEPTION("Sub-tensor " << help2 <<" not available for material type " << help1);
  }

  // Error for material data is not constant, but was queried via Get{Scalar,Vector,Tensor}
  void BaseMaterial::matDataNotConstant(MaterialType matType) const {
    EXCEPTION( "Material data ' " << MaterialTypeEnum.ToString(matType)
               << "' of material '" << name_ << "' is not constant!");
  }

  void BaseMaterial::StoreTensor(PtrParamNode in, PtrCoefFct tensorFunc)
  {
    LocPointMapped lpm;

    if (tensorFunc->GetDependency() == CoefFunction::CONSTANT) {
      // get the tensor by default as complex
      if (tensorFunc->IsComplex()) {
        Matrix<Complex> mat;
        tensorFunc->GetTensor(mat, lpm);
        if (mat.GetPart(Global::IMAG).NormL2() == 0.0)
          in->SetValue(mat.GetPart(Global::REAL));
        else
          in->SetValue(mat);
      }
      else {
        Matrix<Double> mat;
        tensorFunc->GetTensor(mat, lpm);
        in->SetValue(mat);
      }
    }
    else {
      shared_ptr<CoefFunctionAnalytic> analTensor =
          dynamic_pointer_cast<CoefFunctionAnalytic>(tensorFunc);
      if (analTensor) {
        UInt numRows, numCols;
        StdVector<std::string> real, imag;
        analTensor->GetStrTensor(numRows, numCols, real, imag);

        in->Get("dim1")->SetValue(numRows);
        in->Get("dim2")->SetValue(numCols);
        in->Get("real")->SetValue(real);
        in->Get("imag")->SetValue(imag);
      }
    }
  }


  void BaseMaterial::ToInfo(PtrParamNode in, SubTensorType stt, const Vector<Double>* rot)
  {
    set<MaterialType>::iterator iter;

    // in isSet the actually read material parameters are stored
    for (iter = isSet_.begin(); iter != isSet_.end(); ++iter) {
      MaterialType mt = *iter;

      //check for material data
      CoefMap::const_iterator posScal = scalarCoef_.find(mt);
      CoefMap::const_iterator posVec  = vectorCoef_.find(mt);
      CoefMap::const_iterator posTens = tensorCoef_.find(mt);
      map<MaterialType, string >::const_iterator  posStr = stringParams_.find(mt);
      map<MaterialType, Integer >::const_iterator posInt = integerParams_.find(mt);

      PtrParamNode in_ = in->Get("property", ParamNode::APPEND);
      in_->Get("name")->SetValue(MaterialTypeEnum.ToString(mt));

      if (posScal != scalarCoef_.end()) {
        in_->Get("value")->SetValue(posScal->second->ToString());
      }

      if (posVec != vectorCoef_.end()) {
        in_->Get("value")->SetValue(posVec->second->ToString());

        if (stt != FULL && stt != NO_TENSOR) {
          PtrCoefFct subVec = GetSubVectorCoefFnc(mt, stt, Global::COMPLEX);
          in_->Get("subvector")->SetValue(subVec->ToString());
        }
      }

      if (posTens != tensorCoef_.end()) {
        StoreTensor(in_->Get("tensor"), posTens->second);

        // e.g. in the flatShellPlateEV test case we have NO_TENSOR which cannot be
        // handled by ComputeSubTensor()
        if (stt != FULL && stt != NO_TENSOR) { // electrostatic is NO_TENSOR
          PtrCoefFct subTensor = GetSubTensorCoefFnc(mt, stt, Global::COMPLEX);
          StoreTensor(in_->Get("subtensor"), subTensor);
        }
      }

      if (posStr != stringParams_.end()) {
        in_->Get("value")->SetValue(posStr->second);
      }

      if (posInt != integerParams_.end()) {
        in_->Get("value")->SetValue(posInt->second);
      }

      if (rot != NULL && rot->NormMax() > 0) {
        assert(rot->GetSize() == 3);
        PtrParamNode rot_ = in_->Get("rotation");
        rot_->Get("alpha")->SetValue((*rot)[0]);
        rot_->Get("beta")->SetValue((*rot)[1]);
        rot_->Get("gamma")->SetValue((*rot)[2]);
      }
    }
  }


  void BaseMaterial::RotateTensorByRotationAngles( const Vector<Double> &rotAngle,
                                                   MaterialType matType,
                                                   bool persistent) {
    // Calculate rotation matrix (based on Kardan-Angles)
    // Ref.: C. Woernle, "Skript: Dynamik von Mehrkoerpersystemen,
    // Kapitel 2 "Grundlagen der Kinematik", S. 12, Univ. Rostock
    // http://iamserver.fms.uni-rostock.de/studium/mehrkoerpersysteme/unterlagen.htm
    // This rotates a tensor first around the z-axis by gamma, then around the
    // y-axis by beta and last around the x-axis by alpha
    // Kardan-Angles, x-y-z (extrinsic rotation)

    Double alpha = rotAngle[0] * M_PI / 180.0;
    Double beta  = rotAngle[1] * M_PI / 180.0;;
    Double gamma = rotAngle[2] * M_PI / 180.0;;
    Matrix<Double> R(3,3);
    R.Resize(3,3);

    R[0][0] =  cos(beta) * cos(gamma);
    R[0][1] = -cos(beta) * sin(gamma);
    R[0][2] =  sin(beta);
    R[1][0] =  cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);
    R[1][1] =  cos(alpha)*cos(gamma) - sin(alpha)*sin(beta)*sin(gamma);
    R[1][2] = -sin(alpha)*cos(beta);
    R[2][0] =  sin(alpha)*sin(gamma) - cos(alpha)*sin(beta)*cos(gamma);
    R[2][1] =  sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
    R[2][2] =  cos(alpha)*cos(beta);
   
    RotateTensorByRotationMatrix( R, matType, persistent );
  }
  
  void BaseMaterial::RotateAllTensorsByRotationAngles( const Vector<Double>& rotAngle, 
                                                       bool persistent )
  {
    // Rotate all coefficient tensors
    BaseMaterial::CoefMap::iterator cIt = tensorCoef_.begin();
    for( ; cIt != tensorCoef_.end(); cIt++ ) {
      RotateTensorByRotationAngles( rotAngle, cIt->first, persistent );
    }
  }
  
  void BaseMaterial::RotateTensorByRotationMatrix( const Matrix<Double>& rotMatrix, 
                                                   MaterialType matType,
                                                   bool persistent ) {
    if( tensorCoef_.find( matType) != tensorCoef_.end() ) {
      PtrCoefFct coef;
      PtrCoefFct coefOrig = tensorOrigCoef_[matType];

      // perform rotation. In case the matrix is constant (99% of all cases)
      // we resort to an alternative variant
      assert(coosy_);
      if (coefOrig->GetDependency() == CoefFunction::CONSTANT /*&&
          coosy_->HasConstantRotMatrix()*/) { // TODO
        PerformRotationConst( rotMatrix, coef, coefOrig );
      } else {
        PerformRotation( rotMatrix, coef, coefOrig );
      }

      // store back rotated material
      tensorCoef_[matType] = coef;
      
      // In the case of a persistent rotation, we override the
      // original tensor as well
      if ( persistent ) {
        tensorOrigCoef_[matType] = coef;
      }
    } else {
      string dim = "tensor";
      matTypeNotInDataBase( matType, dim );
    }
  }

  void BaseMaterial::RotateAllTensorsByRotationMat(  const Matrix<Double>& rotMatrix,
                                                     bool persistent ) {
    // Rotate all coefficient tensors
    BaseMaterial::CoefMap::iterator cIt = tensorCoef_.begin();
    for( ; cIt != tensorCoef_.end(); cIt++ ) {
      RotateTensorByRotationMatrix( rotMatrix, cIt->first, persistent );
    }
  }

  // Rotate all tensor material parameters by given coordinate system
  // and axis mapping.
  void BaseMaterial::RotateAllTensorsByCoordSys(CoordSystem *coordSys,
                                          const StdVector<std::string> &axisMap,
                                          const StdVector<Double> &axisFactors,
                                          bool persistent)
  {
    // TODO
    /*BaseMaterial::CoefMap::iterator coefIt = tensorCoef_.begin(),
                                    endIt = tensorCoef_.end();
    
    if( coordSys->HasConstantRotMatrix()) {
      Matrix<Double> rotMatrix, invRotMatrix;
      rotMatrix = coordSys->GetRotMatrixFromDirVec(axisMap, axisFactors);
      rotMatrix.Transpose(invRotMatrix);
      for ( ; coefIt != endIt; ++coefIt) {
        RotateTensorByRotationMatrix(invRotMatrix, coefIt->first, persistent);
      }
    }
    else {
      for ( ; coefIt != endIt; ++coefIt) {
        coefIt->second->SetCoordinateSystem(coordSys);
        coefIt->second->SetAxisMapping(axisMap, axisFactors);
      }
    }*/
  }

  void BaseMaterial::RotateTensorByPointCoord( const Vector<Double>&  coord,
                                               MaterialType matType ) {


    // Determine rotation angles from attached coordinate system
    Matrix<Double> rotMatrix;
    coosy_->GetFullGlobRotationMatrix( rotMatrix, coord );

    // Calculate rotation. In this case this is always
    // non-persistent, as the orientation may change for each point in the
    // related coordinate system
    RotateTensorByRotationMatrix( rotMatrix, matType, false );

  }

  // Pass coordinate system to material
  void BaseMaterial::SetCoordSys( CoordSystem* system ) {
    coosy_ = system;
  }

  // Get coordinate system from material
  CoordSystem* BaseMaterial::GetCoordSys() {
    return coosy_;
  }

  void BaseMaterial::PerformRotation( const Matrix<Double>& R, PtrCoefFct& rotatedCoef,
                                       PtrCoefFct origCoef) {
    
    // determine entry type (double / complex) of coefficient
    Global::ComplexPart part = 
        origCoef->IsComplex() ? Global::COMPLEX : Global::REAL;
    
    // get memory for transposed rotation matrix
    Matrix<Double> RT;
    RT.Resize(3,3);
    R.Transpose(RT);
    
    // convert both to coefficient functions
    shared_ptr<CoefFunctionConst<Double> > cR(new CoefFunctionConst<Double>());
    shared_ptr<CoefFunctionConst<Double> > cRT(new CoefFunctionConst<Double>());
    cR->SetTensor( R );
    cRT->SetTensor( RT );

    //get dimension of matrix
    UInt rowSize, colSize;
    origCoef->GetTensorSize(rowSize, colSize );

    Matrix<Complex> helpMat;

    // perform rotation of 3 x 3 matrix
    if ( rowSize == 3 && colSize == 3) {
      // tensor is a 3x3 matrix: sol = R * matrixOrig * RT
      CoefXprBinOp tmp( mp_, origCoef, cRT, CoefXpr::OP_MULT );
      CoefXprBinOp final( mp_, cR, tmp, CoefXpr::OP_MULT );
      rotatedCoef = CoefFunction::Generate( mp_, part, final );
    }
    else if ((rowSize == 3 && colSize == 6) ||
             (rowSize == 6 && colSize == 6))
    {
      // we also need Q;
      Matrix<Double> Q;

      // Composed Rotation Matrix
      // Ref.: M.Richter, "Entwicklung mechanischer Modelle zur analytischen
      // Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton",
      // Diss., Dresden, 2005, p. 27

      Q.Resize(6,6);  

      Q[0][0] = R[0][0]*R[0][0];
      Q[0][1] = R[0][1]*R[0][1];
      Q[0][2] = R[0][2]*R[0][2];
      Q[0][3] = 2.0*R[0][1]*R[0][2];
      Q[0][4] = 2.0*R[0][0]*R[0][2];
      Q[0][5] = 2.0*R[0][0]*R[0][1];

      Q[1][0] = R[1][0]*R[1][0];
      Q[1][1] = R[1][1]*R[1][1];
      Q[1][2] = R[1][2]*R[1][2];
      Q[1][3] = 2.0*R[1][1]*R[1][2];
      Q[1][4] = 2.0*R[1][0]*R[1][2];
      Q[1][5] = 2.0*R[1][0]*R[1][1];

      Q[2][0] = R[2][0]*R[2][0];
      Q[2][1] = R[2][1]*R[2][1];
      Q[2][2] = R[2][2]*R[2][2];
      Q[2][3] = 2.0*R[2][1]*R[2][2];
      Q[2][4] = 2.0*R[2][0]*R[2][2];
      Q[2][5] = 2.0*R[2][0]*R[2][1];

      Q[3][0] = R[1][0]*R[2][0];
      Q[3][1] = R[1][1]*R[2][1];
      Q[3][2] = R[1][2]*R[2][2];
      Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
      Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0];
      Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];

      Q[4][0] = R[0][0]*R[2][0];
      Q[4][1] = R[0][1]*R[2][1];
      Q[4][2] = R[0][2]*R[2][2];
      Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
      Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0];
      Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];

      Q[5][0] = R[0][0]*R[1][0];
      Q[5][1] = R[0][1]*R[1][1];
      Q[5][2] = R[0][2]*R[1][2];
      Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
      Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0];
      Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];


      //  std::cout << "R:\n" << R << std::endl;
      //  std::cout << "Q:\n" << Q << std::endl;
      //  std::cout << "Tensor orig:\n" << matTensor << std::endl;

      Matrix<Double> QT;
      QT.Resize(6,6);
      Q.Transpose(QT);
      
      // convert both to coefficient functions
      shared_ptr<CoefFunctionConst<Double> > cQ(new CoefFunctionConst<Double>());
      shared_ptr<CoefFunctionConst<Double> > cQT(new CoefFunctionConst<Double>());
      cQ->SetTensor( Q );
      cQT->SetTensor( QT );

      if ( rowSize == 3 && colSize == 6 ) {
        CoefXprBinOp tmp( mp_, origCoef, cQT, CoefXpr::OP_MULT );
        UInt nRows, nCols;
        StdVector<std::string>t1, t2;
        tmp.GetTensorXpr(nRows, nCols, t1, t2 );
        CoefXprBinOp final( mp_, cR, tmp, CoefXpr::OP_MULT );
        final.GetTensorXpr(nRows, nCols, t1, t2 );
        rotatedCoef = CoefFunction::Generate(mp_,  part, final );
      }
      else if (rowSize == 6 && colSize == 6 ) {
        CoefXprBinOp tmp( mp_, origCoef, cQT, CoefXpr::OP_MULT );
        PtrCoefFct tmp2 = CoefFunction::Generate(mp_,  part, tmp );
        CoefXprBinOp final( mp_, cQ, tmp, CoefXpr::OP_MULT );
        rotatedCoef = CoefFunction::Generate(mp_, part, final );
      
      }
    }
    else {
      EXCEPTION("Tensor rotation currently only works for 3D matrices!");
    }
  }
  
  void BaseMaterial::PerformRotationConst(const Matrix<Double>& R,
                                          PtrCoefFct& rotatedCoef,
                                          PtrCoefFct origCoef)
  {

    if (origCoef->GetDependency() != CoefFunction::CONSTANT) {
      EXCEPTION( "This type of rotation only works with constant coefficient matrices");
    }

    LocPointMapped lpm; // can be arbitrary

    if (origCoef->IsComplex()) {
      // Obtain constant matrix
      Matrix<Complex> origMat, rotMat;
      origCoef->GetTensor(origMat, lpm);

      // Do the rotation
      origMat.PerformRotation(R, rotMat);

      // Generate a constant coefficient matrix
      shared_ptr<CoefFunctionConst<Complex> > rotCoefConst(new CoefFunctionConst<Complex>());
      rotCoefConst->SetTensor(rotMat);
      rotatedCoef = rotCoefConst;
    }
    else {
      // Obtain constant matrix
      Matrix<Double> origMat, rotMat;
      origCoef->GetTensor(origMat, lpm);

      // Do the rotation
      origMat.PerformRotation(R, rotMat);

      // Generate a constant coefficient matrix
      shared_ptr<CoefFunctionConst<Double> > rotCoefConst(new CoefFunctionConst<Double>());
      rotCoefConst->SetTensor(rotMat);
      rotatedCoef = rotCoefConst;
    }
  }

  template<typename T>
  void BaseMaterial::PerformRotationVoigt( const Matrix<Double>& rotMatrix,
                                           Vector<T>& rotMatTensor,
                                           const Vector<T>& origMatTensor ) {

    // ensure that only 6x6 vector gets entered
    if( origMatTensor.GetSize() != 6 ) {
      EXCEPTION( "Only 6-element tensor can be rotated" );
    }

    Vector<T>& ret = rotMatTensor;
    const Vector<T>& s = origMatTensor;


    // Note: We have to transpose the rotation matrix, as we want to calculate
    //       R^T * [sigma] * R
    // but the formula below computes R * [sigma] * R^T, as given in
    // http://en.wikipedia.org/wiki/Cauchy_stress_tensor#Transformation_rule_of_the_stress_tensor

    Matrix<Double> r;
    rotMatrix.Transpose(r);
    ret.Resize(6);

    // sigma_11
    ret[0] = (r[0][0]*r[0][0])*s[0] + (r[0][1]*r[0][1])*s[1] + (r[0][2]*r[0][2])*s[2]
             + 2.*( r[0][0]*r[0][1]*s[5]
             + r[0][0]*r[0][2]*s[4]
             + r[0][1]*r[0][2]*s[3]);

    // sigma_22
    ret[1] = (r[1][0]*r[1][0])*s[0] + (r[1][1]*r[1][1])*s[1] + (r[1][2]*r[1][2])*s[2]
             + 2.*(r[1][0]*r[1][1]*s[5]
             + r[1][0]*r[1][2]*s[4]
             + r[1][1]*r[1][2]*s[3]);

    // sigma_33
    ret[2] = (r[2][0]*r[2][0])*s[0] + (r[2][1]*r[2][1])*s[1] + (r[2][2]*r[2][2])*s[2]
             + 2.*(r[2][0]*r[2][1]*s[5]
             + r[2][0]*r[2][2]*s[4]
             + r[2][1]*r[2][2]*s[3]);

    // sigma_23
    ret[3] = r[1][0]*r[2][0]*s[0] + r[1][1]*r[2][1]*s[1] + r[1][2]*r[2][2]*s[2]
             + (r[1][0]*r[2][1] + r[1][1]*r[2][0]) * s[5]
             + (r[1][1]*r[2][2] + r[1][2]*r[2][1]) * s[3]
             + (r[1][0]*r[2][2] + r[1][2]*r[2][0]) * s[4];

    // sigma_13
    ret[4] = r[0][0]*r[2][0]*s[0] + r[0][1]*r[2][1]*s[1] + r[0][2]*r[2][2]*s[2]
             + (r[0][0]*r[2][1] + r[0][1]*r[2][0]) * s[5]
             + (r[0][1]*r[2][2] + r[0][2]*r[2][1]) * s[3]
             + (r[0][0]*r[2][2] + r[0][2]*r[2][0]) * s[4];

    // sigma_12
    ret[5] = r[0][0]*r[1][0]*s[0] + r[0][1]*r[1][1]*s[1] + r[0][2]*r[1][2]*s[2]
             + (r[0][0]*r[1][1] + r[0][1]*r[1][0]) * s[5]
             + (r[0][1]*r[1][2] + r[0][2]*r[1][1]) * s[3]
             + (r[0][0]*r[1][2] + r[0][2]*r[1][0]) * s[4];

    //    std::cerr << "origTensor: " << s.Serialize() << std::endl;
    //    std::cerr << "rotatedTensor: " << ret.Serialize() << std::endl;

  }
  template void BaseMaterial::
  PerformRotationVoigt<Complex>( const Matrix<Double>& rotMatrix,
                            Vector<Complex>& rotMatTensor,
                            const Vector<Complex>& origMatTensor );

  template void BaseMaterial::
  PerformRotationVoigt<Double>( const Matrix<Double>& rotMatrix,
                           Vector<Double>& rotMatTensor,
                           const Vector<Double>& origMatTensor );

  void BaseMaterial::CalcFull3x3Tensor( MaterialType isoProp,
                                        MaterialType* orthoProp,
                                        MaterialType tensorProp ) {

    PtrCoefFct val1, val2, val3, isoVal, valTensor;
    StdVector<PtrCoefFct> tensorComp(9);

    // depending on symmetry, calculate full 3x3 permeability tensor
    switch(GetSymmetryType(tensorProp)) {

      case GENERAL:
        // in this case we have already the full material tensor
        break;

      case ISOTROPIC:
        isoVal = GetScalCoefFnc( isoProp, Global::COMPLEX );
        // set diagonal entries
        tensorComp[0] = isoVal;
        tensorComp[4] = isoVal;
        tensorComp[8] = isoVal;
        valTensor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                           3, 3, tensorComp );
        SetCoefFct( tensorProp, valTensor );
        break;

      case TRANS_ISOTROPIC:
        val1 = GetScalCoefFnc( isoProp, Global::COMPLEX );
        val3 = GetScalCoefFnc( orthoProp[2], Global::COMPLEX );
        // set diagonal entries
        tensorComp[0] = val1;
        tensorComp[4] = val1;
        tensorComp[8] = val3;
        valTensor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                           3, 3, tensorComp );
        SetCoefFct( tensorProp, valTensor );
        break;

      case ORTHOTROPIC:

        val1 = GetScalCoefFnc( orthoProp[0], Global::COMPLEX );
        val2 = GetScalCoefFnc( orthoProp[1], Global::COMPLEX );
        val3 = GetScalCoefFnc( orthoProp[2], Global::COMPLEX );
        tensorComp[0] = val1;
        tensorComp[4] = val2;
        tensorComp[8] = val3;
        valTensor = CoefFunction::Generate( mp_, Global::COMPLEX,
                                           3, 3, tensorComp );
        SetCoefFct( tensorProp, valTensor );
        break;

      default:
        EXCEPTION( "Calculation of full matrix for property '"
            << MaterialTypeEnum.ToString( tensorProp )
            << "' with symmetry type '"
            << SymmetryTypeEnum.ToString(GetSymmetryType(tensorProp))
            << "' not implemented!" );
    }
  }


  void BaseMaterial::InitHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                               bool isInverse, bool computeHystInverse, UInt dim ) {

    /*
     * is this function ever called?
     * -> grep shows NO call to InitHyst;
     *    instead everything is handled via CoefFunctionHyst
     * -> BUT linker needs some reference to Preisach and VectorPreisach
     *    if we comment out the calls below (even though they are never used)
     *    the classes are loaded in correct order for linker
     */
    // just make dummy calls for linker
    ParameterPreisachOperators operatorParams = ParameterPreisachOperators();
    ParameterPreisachWeights weightParams = ParameterPreisachWeights();
    bool isVirgin = false;
    bool ignoreAnhystPart = false;
    Integer numElem = 1;

    hyst_ = new Preisach(numElem,operatorParams,weightParams,isVirgin,ignoreAnhystPart);
    hyst_ = new VectorPreisachSutor(numElem,operatorParams,weightParams,dim,isVirgin);
    hyst_ = new VectorPreisachMayergoyz(numElem,operatorParams,weightParams,dim,isVirgin);
    EXCEPTION( "BaseMaterial::InitHyst should not be used anymore" );
//
//    isHystInverse_      = isInverse;
//    computeHystInverse_ = computeHystInverse;
//
//    string val = stringParams_[HYST_MODEL];
//    if ( val != "preisach" ) {
//      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
//    }
//    else {
//     // EXCEPTION( "BaseMaterial::InitHyst should not be used anymore" );
////
////      isHysteresis_ = true;
////
////      Double Xsat, Ysat;
////      GetScalar(Xsat, X_SATURATION, Global::REAL);
////      GetScalar(Ysat, Y_SATURATION, Global::REAL);
////      Matrix<Double> weights;
////      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
////      bool isVirgin = true;
////
////      if(dim == 1){
////
////        hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);
////
////      } else if(dim > 1 && dim <= 3){
////
////        Double rotationalResistance = 1.0;
////        GetScalar(rotationalResistance, ROT_RESISTANCE, Global::REAL);
////
////        int evalVersion;
////        GetScalar(evalVersion, EVAL_VERSION);
////
////        int isTesting;
////        GetScalar(isTesting, IS_TESTING);
////
////        Double angDistance;
////        Matrix<Double> easyAxis_Matrix;
////        Vector<Double> easyAxis = Vector<Double>(dim);
////        GetScalar(angDistance, ANG_DISTANCE, Global::REAL);
////
////      /*
////       * should be obsolete as hyst_ is initialized in coefFctHyst
////       */
////
////      bool classical;
////
////      if(evalVersion == 1){
////        classical = true; // original vector preisach model -> sutor2012
////
////        hyst_ = new VectorPreisachSutor_ListApproach(numElemSD, Xsat, Ysat,
////                                                   weights, rotationalResistance, dim_, isVirgin,
////                                                   classical, angDistance, 0, 0, 0, 0, false);
////      } else if(evalVersion == 2){
////        classical = false; // revised vector preisach model -> sutor2015
////
////        hyst_ = new VectorPreisachSutor_ListApproach(numElemSD, Xsat, Ysat,
////                                                   weights, rotationalResistance, dim_, isVirgin,
////                                                   classical, angDistance, 0, 0, 0, 0, false);
////      } else if(evalVersion == 10){
////        classical = true; // original vector preisach model -> sutor2015; matrix based implementation
////
////        hyst_ = new VectorPreisachSutor_MatrixApproach(numElemSD, Xsat, Ysat,
////                                                   weights, rotationalResistance, dim_, isVirgin,
////                                                   classical, angDistance, 0, 0, 0, 0, false);
////      } else if(evalVersion == 20){
////        classical = false; // revised vector preisach model -> sutor2015; matrix based implementation
////
////        hyst_ = new VectorPreisachSutor_MatrixApproach(numElemSD, Xsat, Ysat,
////                                                   weights, rotationalResistance, dim_, isVirgin,
////                                                   classical, angDistance, 0, 0, 0, 0, false);
////      } else {
////        EXCEPTION("evalVersion has to be one of the following: \n "
////            "1: classical vector model (sutor2012) \n"
////            "2: revised vector model (sutor2015) [DEFAULT] \n"
////            "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
////            "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
////      }
////
//////        if((evalVersion == 7)||(evalVersion == 8)){
//////          hyst_ = new VectorPreisachv7(numElemSD, Xsat, Ysat, weights,rotationalResistance,dim, isVirgin, isTesting!=0, (UInt) evalVersion);
//////        } else if((evalVersion == 9)||(evalVersion == 10)){
//////		  Vector<Double> easyAxis = Vector<Double>(dim_);
//////		  Double phaseLag = 0.0;
//////		  hyst_ = new VectorPreisachSutor(numElemSD, Xsat, Ysat, weights,rotationalResistance,dim, isVirgin, isTesting!=0, (UInt) evalVersion,phaseLag,easyAxis);
//////        } else {
//////          hyst_ = new VectorPreisach(numElemSD, Xsat, Ysat, weights,rotationalResistance,dim, isVirgin, isTesting!=0, (UInt) evalVersion);
//////        }
////
////      }
////
////
////      // set map: global to local element number
////      EntityIterator it = actSDList->GetIterator();
////      UInt iel = 0;
////      UInt globalElNr;
////      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
////
////      globalElNr = it.GetElem()->elemNum;
////      globalElem2Local_[globalElNr] = iel;
////      }
////    }
////
////    //allocate memory for previous results, needed for the
////    //effective material parameter formulation
////    if(dim == 1){
////      Xprevious_.Resize(numElemSD);
////      Yprevious_.Resize(numElemSD);
////      Xprevious_.Init();
////      Yprevious_.Init();
////    }  else if(dim > 1 && dim <= 3){
////      XpreviousVEC_ = new Vector<Double>[numElemSD];
////      YpreviousVEC_ = new Vector<Double>[numElemSD];
////
////      for(UInt i = 0; i < numElemSD; i++){
////        XpreviousVEC_[i].Resize(dim_);
////        XpreviousVEC_[i].Init();
////
////        YpreviousVEC_[i].Resize(dim_);
////        YpreviousVEC_[i].Init();
////       }
////
//    }
  }

  void BaseMaterial::InitVecHyst( UInt numElemSD, shared_ptr<ElemList> actSDList,
                                  UInt dim ) {


    string val = stringParams_[HYST_MODEL];
    if ( val != "preisach" ) {
      EXCEPTION( "Currently we just support Preisach Hysteresis Model" );
    }
    else {
      dim_ = dim;
      isHysteresis_ = true;

      Double Xsat, Ysat;
      GetScalar(Xsat, X_SATURATION, Global::REAL);
      GetScalar(Ysat, Y_SATURATION, Global::REAL);
      Matrix<Double> weights;
      GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);
      bool isVirgin = true; 

      dimVecHyst_ = dim;
      hyst_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      if ( dim == 2 ) 
        hystY_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      else if ( dim == 3) {
        hystY_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
        hystZ_ = new SimplePreisachInv(numElemSD, Xsat, Ysat, weights, isVirgin);
      }

      // set map: global to local element number
      EntityIterator it = actSDList->GetIterator();
      UInt iel = 0;
      UInt globalElNr;
      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
        globalElNr = it.GetElem()->elemNum;
        globalElem2Local_[globalElNr] = iel;
      }
    }

    //allocate memory for previous results, needed for the
    //effective material parameter formulation
    vecXprevious_.Resize(dim,numElemSD);
    vecYprevious_.Resize(dim,numElemSD);
    vecXprevious_.Init();
    vecYprevious_.Init();

    actDiffVal4VecHyst_.Resize(dim,numElemSD);
    previousDiffVal4VecHyst_.Resize(dim,numElemSD);
    Double nu;
    GetScalar(nu,MAG_RELUCTIVITY_SCALAR,Global::REAL);
    for (UInt i=0; i<dim; i++)
      for (UInt j=0; j<numElemSD; j++)
        previousDiffVal4VecHyst_[i][j] = nu;

  }


  Double BaseMaterial::GetScalarHystVal( UInt nrElem ) {

    if(dim_ != 1){
       EXCEPTION("Only implemented for scalar model");
     }
    UInt idx = globalElem2Local_[nrElem];
    return hyst_->getValue( idx );
  }


  Double BaseMaterial::GetScalarHystPrevVal( UInt nrElem ) {
    if(dim_ != 1){
       EXCEPTION("Only implemented for scalar model");
     }
    UInt idx = globalElem2Local_[nrElem];
    return Yprevious_[idx];
  }

  void BaseMaterial::GetRayleighCoeffStrings(std::string &alpha, std::string &beta)
  {
    if (this->raylDampType_ == ALPHA_BETA) {
      double a;
      double b;
      GetScalar(a, RAYLEIGH_ALPHA, Global::REAL);
      GetScalar(b, RAYLEIGH_BETA, Global::REAL);
      alpha = lexical_cast<std::string>(a);
      beta = lexical_cast<std::string>(b);
      return;
    }
    else {
      EXCEPTION("Rayleigh Alpha and Beta have not been set!")
    }
  }

  void BaseMaterial::GetFreqAdaptedRayleighCoeffStrings(std::string &alpha, std::string &beta)
  {
    if (this->raylDampType_ == ADAPTED_LOSS_TANGENS) {
      // First, check if the damping information is provided as a constant or frequency dependent
      std::string tanDeltaStr;
      GetString(tanDeltaStr, LOSS_TANGENS_DELTA);
      double deltaF = 0.01;
      // make sure to enclose the expression by brackets!
      std::string deltaFStr = lexical_cast<std::string>(deltaF);
      // simplified formula to keep mathParser string short
      alpha = tanDeltaStr + "* pi * f * (1 - (" + deltaFStr + ")^2)";
      beta = tanDeltaStr + "* 1 / (4 * pi * f)";
      return;
    }
    else {
      EXCEPTION("Adapted Loss Tangens has not been set!")
    }
  }

  MaterialType BaseMaterial::ConvertMaterialClass(MaterialClass mc)
  {
    switch(mc)
    {
    case MECHANIC:
      return MECH_STIFFNESS_TENSOR;
    case PIEZO:
      return PIEZO_TENSOR;
    case ELECTROSTATIC:
      return MECH_STIFFNESS_TENSOR;
    case SMOOTH:
      return SMOOTH_STIFFNESS_TENSOR;
    default:
      assert(false); // implement for your needs!
      break;
    }

    assert(false);
    return NO_MATERIAL;
  }


  void BaseMaterial::InitPiezoMicro( UInt numElemSD, shared_ptr<ElemList> actSDList, 
                                     BaseMaterial* mechMat, BaseMaterial* elecMat,
                                     SubTensorType tensorType, Double dt) {

      isPiezoMicroModel_ = true;

      piezoMicroModel_ = new PiezoMicroModelBK( numElemSD, this, 
                                                mechMat, elecMat, 
                                                tensorType, dt );

      // set map: global to local element number
      EntityIterator it = actSDList->GetIterator();
      UInt iel = 0;
      UInt globalElNr;
      for ( it.Begin(); !it.IsEnd(); it++, iel++) {
        globalElNr = it.GetElem()->elemNum;
        globalElem2Local_[globalElNr] = iel;
      }
  }
  

  void BaseMaterial::GetEffectiveTensors( Matrix<Double>& matMechC,
                                          Matrix<Double>& matMechS,
                                          Matrix<Double>& matElec,
                                          Matrix<Double>& matPiezo,
                                          Vector<Double>& stress, 
                                          Vector<Double>& elecField,
                                          UInt elemIdx, 
                                          bool recompute,
                                          bool previous ) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->GetEffectiveTensors( matMechC, matMechS, 
                                           matElec,  matPiezo, 
                                           stress, elecField, 
                                           idx, recompute, previous );
  }

  void BaseMaterial::GetEffectiveIrreversibleValues( Vector<Double>& Pirr,
                                                     Vector<Double>& Sirr,
                                                     UInt elemIdx,
                                                     bool recompute,
                                                     bool previous ) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->GetEffectiveIrreversibleValues( Pirr, Sirr, idx, recompute, previous );
  }

   void BaseMaterial::ComputeEffectiveCouplingTensor(Matrix<Double>& dMatEff, 
                                                     Vector<Double>& elecFieldAct,
                                                     Vector<Double>& elecFieldPrev,
                                                     UInt elemIdx) {

    UInt idx = globalElem2Local_[elemIdx];

    piezoMicroModel_->ComputeEffectiveCouplingTensor(dMatEff, elecFieldAct,
                                                     elecFieldPrev, idx);
   }


   PtrCoefFct BaseMaterial::GetTensorCoefFnc(MaterialType matType, SubTensorType type,
                                             Global::ComplexPart matDataType,
                                             bool transpose ) const
   {
     if ( tensorCoef_.find(matType) == tensorCoef_.end() ) {
       matTypeNotInDataBase(matType, "tensor");
     } 

     PtrCoefFct mFunct = GetSubTensorCoefFnc( matType, type, matDataType, transpose );
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }
   
   PtrCoefFct BaseMaterial::GetVectorCoefFnc(MaterialType matType,
                                             Global::ComplexPart matDataType) const
   {
     CoefMap::const_iterator it = vectorCoef_.find(matType);
     if ( it == vectorCoef_.end() ) {
       matTypeNotInDataBase(matType, "vector");
     }

     PtrCoefFct mFunct = it->second->GetComplexPart( matDataType );
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }

   PtrCoefFct  BaseMaterial::GetScalCoefFnc(MaterialType matType,
                                            Global::ComplexPart matDataType) const
   {
     CoefMap::const_iterator it = scalarCoef_.find(matType);
     if ( it == scalarCoef_.end() ) {
       matTypeNotInDataBase(matType, "scalar");
     }

     PtrCoefFct mFunct = it->second->GetComplexPart( matDataType );
     mFunct->SetCoordinateSystem(this->coosy_);
     return mFunct;
   }

   PtrCoefFct BaseMaterial::GetSubTensorCoefFnc( MaterialType matType, 
                                                 SubTensorType tensorType,
                                                 Global::ComplexPart matDataType,
                                                 bool transposed ) const
   {
     CoefMap::const_iterator it = tensorCoef_.find(matType);
     if( it ==  tensorCoef_.end() ) {
       matTypeNotInDataBase(matType, "tensor");
     }
     
     if (tensorType == FULL) {
       if( (matType == PREISACH_WEIGHTS) || (matType == PREISACH_WEIGHTS_STRAIN) || (matType == HYST_IRRSTRAIN_CI)){
         // special cases; these tensors are Nx1; need no transpose and can be taken as set in param handler
         // just return tensorCoef_
         return it->second;
       }
     }
 
     CoefXprSubTensor subTensorXpr(mp_, it->second );
     subTensorXpr.SetSubTensorType( tensorType, transposed );
     return CoefFunction::Generate( mp_, matDataType, subTensorXpr );
   }

   PtrCoefFct BaseMaterial::GetSubVectorCoefFnc( MaterialType matType,
                                                 SubTensorType tensorType,
                                                 Global::ComplexPart matDataType) const
   {
       EXCEPTION("Not implemented in base class");
   }

   PtrCoefFct BaseMaterial::GetTensorCoefFncNonLin( MaterialType matType,
                                                    SubTensorType type,
                                                    Global::ComplexPart matDataType,
                                                    PtrCoefFct dependency )
   {
     EXCEPTION("Currently only implemented for ElectroMagnetic material")
   }
   
   PtrCoefFct BaseMaterial::GetScalCoefFncNonLin_MagStrict(MaterialType matType,
                                                           Global::ComplexPart matDataType,
                                                           PtrCoefFct mechStrain )
   {
     EXCEPTION("Only implemented for ElectroMagnetic material")
   }

   PtrCoefFct BaseMaterial::GetScalCoefFncMultivariateNonLin(
       MaterialType matType,
       NonLinType nlType,
       Global::ComplexPart matDataType,
       StdVector<PtrCoefFct> dependencies,
       StdVector<RegionIdType> & regs)
   {
     EXCEPTION("not implemented in base class");
   }


  //  PtrCoefFct BaseMaterial::GetScalCoefFncModel(shared_ptr<CoefFunction> coefObj){
  //   shared_ptr<CoefFunction> coef;
  //   coef = coefObj;
  //   return coef;
  //  }

  //  PtrCoefFct BaseMaterial::GetTensorCoefFncModel(shared_ptr<CoefFunction> coefObj){
  //   shared_ptr<CoefFunction> coef;
  //   coef = coefObj;
  //   return coef;
  //  }

  //  PtrCoefFct BaseMaterial::GetVectorCoefFncModel(shared_ptr<CoefFunction> coefObj){
  //   shared_ptr<CoefFunction> coef;
  //   coef = coefObj;
  //   return coef;
  //  }

   PtrCoefFct BaseMaterial::GetScalCoefFncNonLin(MaterialType matType,
                                                 Global::ComplexPart matDataType,
                                                 PtrCoefFct dependency,
                                                 PtrCoefFct temperature_dependency)
   {
     shared_ptr<CoefFunctionApprox> coef;
     
     // Ensure that only real-valued parameters are used
     if( matDataType != Global::REAL ) {
       EXCEPTION( "Only real-valued nonlinear parameters are supported");
     }
     
     // check if isotropic material type is defined
     if( nonlinIsoParams_.find(matType) == nonlinIsoParams_.end() ) {
       EXCEPTION( "No nonlinear definition found for material type '"
           << MaterialTypeEnum.ToString(matType) << "'");
     }
     
     // check, if nonlinear curve was already calculated
     NonLinIsoMap::iterator it = nonlinIsoParams_.find(matType);
     if(it == nonlinIsoParams_.end() ) {
       EXCEPTION( "Could not find material type '"
           << MaterialTypeEnum.ToString( matType) << "' in material '"
           << name_ << "'");
     }
     MatDescriptorNl & matNl = it->second;
     
     // Check, if smooth spline approximation was already created 
     // and initialized
     if( !matNl.approxData ) {
       if ( matNl.approxType == SMOOTH_SPLINES ) {
         SmoothSpline * sp = new SmoothSpline( matNl.fileName, matType );
         sp->SetAccuracy( matNl.measAccuracy );
         sp->SetMaxY( matNl.maxVal );
         sp->CalcBestParameter();
         sp->CalcApproximation();
         sp->Print();
         matNl.approxData = sp;
       }
       else if ( matNl.approxType == LIN_INTERPOLATE ) {
         LinInterpolate * sp = new LinInterpolate( matNl.fileName, matType );
         //sp->SetAccuracy( matNl.measAccuracy );
         //sp->SetMaxY( matNl.maxVal );
         sp->Print();
         matNl.approxData = sp;
        }
       else {
         EXCEPTION( "No nonlinear approx/interpolate type not available '"
             << ApproxCurveTypeEnum.ToString(matNl.approxType) << "'");
       }
     }

     ApproxData * sp = matNl.approxData;
     
     // get linear starting value
     Double startVal = 0.0;
     this->GetScalar( startVal, matType, Global::REAL );
     coef.reset(new CoefFunctionApprox());
     coef->Init( startVal, sp, dependency );

     return coef;
   }

   // Definition of available tensor symmetry types
   static EnumTuple symTypeTuples[] =
   {
    EnumTuple(BaseMaterial::NOSYMMETRY, "noSymmetry"),
    EnumTuple(BaseMaterial::GENERAL, "general"),
    EnumTuple(BaseMaterial::ISOTROPIC, "isotropic"),
    EnumTuple(BaseMaterial::ORTHOTROPIC, "orthotropic"),
    EnumTuple(BaseMaterial::TRANS_ISOTROPIC, "transversalIsotropic")
   };
   Enum<BaseMaterial::SymmetryType> BaseMaterial::SymmetryTypeEnum = \
       Enum<BaseMaterial::SymmetryType>("Tensor symmetry types",
                                        sizeof(symTypeTuples) / sizeof(EnumTuple),
                                        symTypeTuples);

   // Definition of available material models
   static EnumTuple matModelTuples[] = {
       EnumTuple(BaseMaterial::MAT_MODEL_LINEAR, "linear"),
       EnumTuple(BaseMaterial::MAT_MODEL_LINEARVISCOELASTIC, "linearViscoElastic")
   };
   Enum<BaseMaterial::MaterialModel> BaseMaterial::MaterialModelEnum =
       Enum<BaseMaterial::MaterialModel>("Material model",
                                         sizeof(matModelTuples) / sizeof(EnumTuple),
                                         matModelTuples);

} // end of namespace
