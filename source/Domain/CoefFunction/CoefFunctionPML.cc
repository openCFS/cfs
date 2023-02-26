// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionPML.cc
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2012
 *       \author   ahueppe
 */
//================================================================================================



#include "CoefFunctionPML.hh"
#include "Utils/mathParser/mathParser.hh"
#include "boost/bind/bind.hpp"
#include "boost/lexical_cast.hpp"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

namespace CoupledField{
  //================================================================================================
  // Base Class
  //================================================================================================
  template<typename T>
  CoefFunctionPMLBase<T>::CoefFunctionPMLBase(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
                                    shared_ptr<EntityList> EntList, StdVector<RegionIdType> pdeDomains): 
                                    CoefFunction() {
    isAnalytic_ = false;
    isComplex_ =  std::is_same<T,Complex>::value;
    dependType_ = GENERAL;
    this->entities_.Push_back(EntList);
    dim_ = entities_[0]->GetGrid()->GetDim();
    //orderCoefFct_ = orderCoefFct;
    name_ = "CoefFunctionPMLBase";
    pmlType_ = DampFunction::NO_TYPE;

    //check if the coefficient function of the speed of sound is real valued and scalar
    if(speedOfSound->GetDimType() != CoefFunction::SCALAR || speedOfSound->IsComplex()){
      EXCEPTION("The PML coefficient function needs a real valued scalar coeffunction as speed of sound");
    }
    speedOfSound_ = speedOfSound;
  
    //this is just to be up to date with the desired frequency!
    // obtain handle from internal variable coefficient function
    mp_ = domain->GetMathParser();
    mHandle_ = mp_->GetNewHandle(true);

    mp_->SetExpr(mHandle_,"f");

    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack(
        boost::bind(&CoefFunctionPMLBase<T>::UpdateOmega, this ),
        mHandle_ );
    // important: Trigger first-time calculation
    UpdateOmega();
  }

  template<typename T>
  CoefFunctionPMLBase<T>::~CoefFunctionPMLBase() {
    //std::cout<<"DESRUC"<<std::endl;
    //It might be necessary to just disconnect the callback instead of releasing the handle
    mp_->ReleaseHandle( mHandle_ );
  }

  template<typename T>
  void CoefFunctionPMLBase<T>::UpdateOmega() {
    omega_ = this->mp_->Eval(mHandle_) * 2 * M_PI;
  }

  template<typename T>
  void CoefFunctionPMLBase<T>::CreateDampFunction(){

    switch(pmlType_){
      case DampFunction::CONSTANT:
        dampFunction_.reset(new DampFunctionConst(1.0));
        break;
      case DampFunction::INVERSE_DIST:
        dampFunction_.reset(new DampFunctionInvDist(1.0));
        break;
      case DampFunction::QUADRATIC:
        dampFunction_.reset(new DampFunctionQuad(1.0));
        break;
      case DampFunction::SMOOTH:
        dampFunction_.reset(new DampFunctionSmooth(1.0));
        break;
      case DampFunction::TANGENS:
        dampFunction_.reset(new DampFunctionTangens());
        break;
      case DampFunction::RATIONAL:
        dampFunction_.reset(new DampFunctionRational());
        break;
      case DampFunction::EXPONENTIAL:
        dampFunction_.reset(new DampFunctionExponential());
        break;
      default:
        EXCEPTION("PML type "<< DampFunction::DampingTypeEnum.ToString(pmlType_) << " not supported");
        break;
    }
  }

  //================================================================================================
  // Classic (=Cartesian) PML
  //================================================================================================
  template<typename T>
  CoefFunctionPML<T>::CoefFunctionPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound,
              shared_ptr<EntityList> EntList, StdVector<RegionIdType> pdeDomains, bool isVector) : 
                  CoefFunctionPMLBase<T>(pmlDef, speedOfSound, EntList, pdeDomains) {
    this->name_ = "CoefFunctionPML";
    //prepare dimenstions of propagation region
    innerMinMaxComp_.Resize(3,2);
    innerMinMaxComp_.Init();

    outerMinMaxComp_.Resize(3,2);
    outerMinMaxComp_.Init();

    isVector_ = isVector;
    this->formulationType_ = CLASSIC;
    ReadDataPML(pmlDef,pdeDomains);

    if (isVector_ ) {
      this->dimType_ = CoefFunction::VECTOR;
    } else {
      this->dimType_ = CoefFunction::SCALAR;
    }
  }

  template<typename T>
  CoefFunctionPML<T>::~CoefFunctionPML(){}


  template<typename T>
  void CoefFunctionPML<T>::GetTensor(Matrix<Complex>& tensor,
                                const LocPointMapped& lpm ){
    //this is diagonal tensor with the coefficients

    tensor.Resize(this->dim_, this->dim_);
    tensor.Init();
    Double locThick=0.0;
    Double position=0.0;
    Complex one(1.0,0.0);
    Double sos;
    this->speedOfSound_->GetScalar(sos,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        Complex fac(0.0,-1.0 * sos* this->dampFunction_->ComputeFactor(position,locThick));
        tensor[i][i] = this->omega_ / (this->omega_ + fac);
      }else{
        tensor[i][i] = one;
      }
    }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetTensor(Matrix<Double>& tensor,
                                const LocPointMapped& lpm ){
    //this is diagonal tensor with the coefficients
    tensor.Resize(this->dim_,this->dim_);
    tensor.Init();
    Double locThick=0.0;
    Double position=0.0;
    Double sos;
    this->speedOfSound_->GetScalar(sos,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        tensor[i][i] = this->dampFunction_->ComputeFactor(position,locThick) * sos;
      }else{
        tensor[i][i] = 0.0;
      }
    }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetVector(Vector<Complex>& vec,
                                const LocPointMapped& lpm ){
    //this is diagonal tensor with the coefficients
    vec.Resize(this->dim_);
    vec.Init();
    Double locThick=0.0;
    Double position=0.0;
    Complex dummy(1.0,0.0);
    Double sos;
    this->speedOfSound_->GetScalar(sos,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        Complex fac(1.0,-1.0 * sos* this->dampFunction_->ComputeFactor(position,locThick)/this->omega_);
        vec[i] = dummy / fac;
      }else{
        vec[i] = dummy;
      }
    }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetVector(Vector<Double>& vec,
                                const LocPointMapped& lpm ){

    //first loop over every entry and determine the factor
    vec.Resize(this->dim_,0.0);
    Double locThick=0.0;
    Double position=0.0;
    Double sos;
    this->speedOfSound_->GetScalar(sos,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        vec[i] = this->dampFunction_->ComputeFactor(position,locThick) * sos;
      }else{
        vec[i] = 0.0;
      }
    }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetScalar(Complex& val,
                                const LocPointMapped& lpm ){
    //computes (1-j pml_x/omega) *(1-j pml_y/omega)*(1-j pml_z/omega))
    Double locThick=0.0;
    Double position=0.0;
    Complex dummy(1.0,0.0);
    val = dummy;
    Double sos;
    this->speedOfSound_->GetScalar(sos,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        Complex fac(1.0,-1.0  *  sos * this->dampFunction_->ComputeFactor(position,locThick)/this->omega_);
        val *=  fac;
      }else{
        val *= dummy;
      }
    }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetScalar(Double& val,
                                const LocPointMapped& lpm ){

    //computes 1/(pml_x*pml_y*pml_z)
    //right now we ust return 1...
    EXCEPTION("GETSCALAR IS INVALID")
    val = 1.0;
    return;

    //UInt gridDim = entities_[0]->GetGrid()->GetDim();
    //Double locThick=0.0;
    //Double position=0.0;
    //val = 1.0;
    //for(UInt i=0;i<gridDim;++i){
    //  GetThicknessAtPoint(locThick,position,lpm,i);
    //  if(abs(locThick)>0.0){
    //    val *= this->dampFunction_->ComputeFactor(position,locThick);
    //  }else{
    //    val *= 1.0;
    //  }
    //}
  }


  template<typename T>
  void CoefFunctionPML<T>::ReadDataPML(PtrParamNode pmlDef,StdVector<RegionIdType> pdeDomains){

    std::string cSysId;
    pmlDef->GetValue("coordSysId",cSysId,ParamNode::PASS);

    if(cSysId == ""){
      cSysId = "default";
    }

    this->SetCoordinateSystem(domain->GetCoordSystem(cSysId));

    //in some cases we can have a propagation node
    PtrParamNode propNode = pmlDef->Get("propRegion",ParamNode::PASS);
    if(propNode){
      //read data from prop node
      ParamNodeList dirNodes = propNode->GetChildren();
      std::string dirName;
      UInt dirIndex = 0;
      for( UInt iDir = 0; iDir < dirNodes.GetSize(); ++iDir ) {
        dirNodes[iDir]->GetValue( "comp", dirName );
        dirIndex = this->coordSys_->GetVecComponent( dirName );
        PtrParamNode ppNode;
        ppNode = dirNodes[iDir]->Get("min");       
        innerMinMaxComp_[dirIndex-1][0] = ppNode->MathParse<Double>();
        ppNode = dirNodes[iDir]->Get("max");       
        innerMinMaxComp_[dirIndex-1][1] = ppNode->MathParse<Double>();
      }
      if(dirNodes.GetSize() != this->dim_)
        GuessLayerData(pdeDomains);
    }else{
      //try to guess the data from grid
      GuessLayerData(pdeDomains);
    }

    //in any case we determine the layer data for the outer bounds
    this->entities_[0]->GetGrid()->CalcBoundingBoxOfRegion(this->entities_[0]->GetRegion(),outerMinMaxComp_,this->coordSys_);

    // => Put the following information to the info xml
    //std::cout << "computed boundingbox of pml region:" << std::endl;
    //for(UInt i=0;i<this->dim_;++i){
    //  std::string vComp = coordSys_->GetDofName(i+1);
    //  std::cout << vComp << "_min: " << outerMinMaxComp_[i][0] << std::endl;
    //  std::cout << vComp << "_max: " << outerMinMaxComp_[i][1] << std::endl;
    //}

    // type of PML damping function
    std::string typeOfPml;
    pmlDef->GetValue("type",typeOfPml);

    this->pmlType_ = DampFunction::DampingTypeEnum.Parse(typeOfPml);

    CreateDampFunction();

    Double dampFactor;
    pmlDef->GetValue("dampFactor",dampFactor);

    this->dampFunction_->DampFactor = dampFactor;
  }

  template<typename T>
  void CoefFunctionPML<T>::GuessLayerData(StdVector<RegionIdType> pdeDomains){
    StdVector<RegionIdType> propReg;
    for(UInt aDom = 0; aDom < pdeDomains.GetSize();++aDom){
      if(pdeDomains[aDom] != this->entities_[0]->GetRegion()){
        propReg.Push_back(pdeDomains[aDom]);
      }
    }

    if( propReg.GetSize() == 0 ){
      EXCEPTION("The PDE  has only PML damped regions. \nThus the "
          << "automatic calculation of damping parameters does not work!\n"
          << "Use the <propRegion> element instead do define the propagation "
          << "region explicitly!");
    }

    //now we check which dimensions are missing and determine from grid
    UInt gridDim = this->entities_[0]->GetGrid()->GetDim();
    Matrix<Double> regBound;
    Double bigVal = 1e40;
    for(UInt i=0;i<gridDim;++i){
      innerMinMaxComp_[i][0] = bigVal;
      innerMinMaxComp_[i][1] = -bigVal;
    }
    for(UInt aDom = 0; aDom < propReg.GetSize();++aDom){
      // ok so we compute the boundingbox of each region
      this->entities_[0]->GetGrid()->CalcBoundingBoxOfRegion(propReg[aDom],regBound,this->coordSys_);
      for(UInt i=0;i<gridDim;++i){
        if(innerMinMaxComp_[i][0]>regBound[i][0])
          innerMinMaxComp_[i][0] = regBound[i][0];

        if(innerMinMaxComp_[i][1]<regBound[i][1])
          innerMinMaxComp_[i][1] = regBound[i][1];

      }
    }
    // Put the following definition to the info xml.
  // std::cout << "computed boundingbox of propagation region:" << std::endl;
  // for(UInt i=0;i<gridDim;++i){
  //   std::string vComp = coordSys_->GetDofName(i+1);
  //   std::cout << vComp << "_min: " << innerMinMaxComp_[i][0] << std::endl;
  //   std::cout << vComp << "_max: " << innerMinMaxComp_[i][1] << std::endl;
  // }
  }

  template<typename T>
  void CoefFunctionPML<T>::GetThicknessAtPoint(Double& thickness, Double & position, LocPointMapped lpm,UInt dir){
    //try to guess where we are at
    thickness = 0.0;
    Vector<Double> globPoint;
    Vector<Double> coordLocalPoint;
    position = 0.0;
    lpm.shapeMap->Local2Global(globPoint,lpm.lp.coord);
    this->coordSys_->Global2LocalCoord(coordLocalPoint,globPoint);

    if(coordLocalPoint[dir] > innerMinMaxComp_[dir][1] &&
      coordLocalPoint[dir] < outerMinMaxComp_[dir][1]){
      thickness = abs(outerMinMaxComp_[dir][1] - innerMinMaxComp_[dir][1]);
      position = abs( coordLocalPoint[dir] - innerMinMaxComp_[dir][1]);
    }
    if(coordLocalPoint[dir] < innerMinMaxComp_[dir][0] &&
      coordLocalPoint[dir] > outerMinMaxComp_[dir][0]){
      thickness = abs(outerMinMaxComp_[dir][0] - innerMinMaxComp_[dir][0]);
      position = abs(coordLocalPoint[dir] - innerMinMaxComp_[dir][0] );
    }

  }

  //================================================================================================
  // DampFunction
  //================================================================================================

  // Definition of types of pml function
  static EnumTuple dampTypeTubles[] = {
    EnumTuple(DampFunction::CONSTANT,       "constant"),
    EnumTuple(DampFunction::INVERSE_DIST,   "inverseDist"),
    EnumTuple(DampFunction::QUADRATIC,      "quadDist"),
    EnumTuple(DampFunction::SMOOTH,         "smoothDamp"),
    EnumTuple(DampFunction::TANGENS,        "tangens"),
    EnumTuple(DampFunction::RATIONAL,       "rational"),
    EnumTuple(DampFunction::EXPONENTIAL,    "exponential")
  };


  Enum<DampFunction::DampingType> DampFunction::DampingTypeEnum = \
    Enum<DampFunction::DampingType>("Types of damping",
        sizeof(dampTypeTubles) / sizeof(EnumTuple),
        dampTypeTubles);

  //================================================================================================
  // Shifted PML
  //================================================================================================

  template<typename T>
  CoefFunctionShiftedPML<T>::CoefFunctionShiftedPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound, shared_ptr<EntityList> EntList,
                                                    StdVector<RegionIdType> pdeDomains, bool isVector)
    : CoefFunctionPML<T>(pmlDef, speedOfSound, EntList, pdeDomains, isVector)
  {
    this->name_ = "CoefFunctionShiftedPML";
    std::string scalingCoefStr, shiftCoefStr;
    Double scalingPow = 0.0, shiftPow = 0.0;

    PtrParamNode scalingNode = pmlDef->Get("scalingCoef", ParamNode::PASS);
    if (scalingNode)
    {
      scalingNode->GetValue("value", scalingCoefStr);
      scalingNode->GetValue("power", scalingPow);
      scalingFunc_.reset(new DampFunctionPolyDirect(scalingPow));
      scalingCoef_ = CoefFunction::Generate(this->mp_, Global::REAL, scalingCoefStr);
    }

    PtrParamNode shiftNode = pmlDef->Get("frqShiftCoef", ParamNode::PASS);
    if (shiftNode)
    {
      shiftNode->GetValue("value", shiftCoefStr);
      shiftNode->GetValue("power", shiftPow);
      shiftFunc_.reset(new DampFunctionPolyInverse(shiftPow));
      shiftCoef_ = CoefFunction::Generate(this->mp_, Global::REAL, shiftCoefStr);
    }

    this->formulationType_ = SHIFTED;
  }

  template<typename T>
  CoefFunctionShiftedPML<T>::~CoefFunctionShiftedPML() { }

  template<typename T>
  void CoefFunctionShiftedPML<T>::GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm)
  {
    // this is diagonal tensor with the coefficients 1/s

    tensor.Resize(this->dim_, this->dim_);
    tensor.Init();
    Vector<Complex> vector;
    GetVector(vector, lpm);

    for (UInt i = 0; i < this->dim_; ++i)
    {
      tensor[i][i] = vector[i];
    }
  }

  template<typename T>
  void CoefFunctionShiftedPML<T>::GetVector(Vector<Complex>& vector, const LocPointMapped& lpm)
  {
    // vector of the form { 1/sx, 1/sy, 1/sz }

    vector.Resize(this->dim_);
    vector.Init();
    Double locThick = 0.0;
    Double position = 0.0;
    Complex dummy(1.0, 0.0);

    // s = kappa + sigma/(alpha + i*w) =
    //   = kappa + alpha*sigma/(alpha^2 + w^2) - i*w*sigma/(alpha^2 + w^2) = kappa + alpha*frac - i*w*frac
    Double sos, alpha, kappa, sigma, frac;
    this->speedOfSound_->GetScalar(sos, lpm);

    Double kappa0 = 0.0, alpha0 = 0.0;
    scalingCoef_->GetScalar(kappa0, lpm);
    shiftCoef_->GetScalar(alpha0, lpm);
    // These two values are to be used as coefficients in scaling and frequency shifting functions:
    // kappa(x) = 1 + kappa0*Fk(x); alpha(x) = alpha0*Fa(x)
    scalingFunc_->DampFactor = kappa0;
    shiftFunc_->DampFactor = alpha0;;

    for (UInt i = 0; i < this->dim_; ++i)
    {
      this->GetThicknessAtPoint(locThick, position, lpm, i);
      if (abs(locThick) > 0.0)
      {
        sigma = sos*this->dampFunction_->ComputeFactor(position, locThick);
        kappa = 1.0 + scalingFunc_->ComputeFactor(position, locThick);
        alpha = shiftFunc_->ComputeFactor(position, locThick);
        frac = sigma/(alpha*alpha + this->omega_*this->omega_);
        Complex fac(kappa + alpha*frac, -1.0*this->omega_*frac);
        vector[i] = dummy/fac;
      }
      else
      {
        vector[i] = dummy;
      }
    }
  }

  template<typename T>
  void CoefFunctionShiftedPML<T>::GetScalar(Complex& scalar, const LocPointMapped& lpm)
  {
    // computes sx*sy*sz

    Double locThick = 0.0;
    Double position = 0.0;
    Complex dummy(1.0, 0.0);
    scalar = dummy;

    // s = kappa + sigma/(alpha + i*w) =
    //   = kappa + alpha*sigma/(alpha^2 + w^2) - i*w*sigma/(alpha^2 + w^2) = kappa + alpha*frac - i*w*frac
    Double sos, alpha, kappa, sigma, frac;
    this->speedOfSound_->GetScalar(sos, lpm);

    Double kappa0 = 0.0, alpha0 = 0.0;
    scalingCoef_->GetScalar(kappa0, lpm);
    shiftCoef_->GetScalar(alpha0, lpm);
    // These two values are to be used as coefficients in scaling and frequency shifting functions:
    // kappa(x) = 1 + kappa0*Fk(x); alpha(x) = alpha0*Fa(x)
    scalingFunc_->DampFactor = kappa0;
    shiftFunc_->DampFactor = alpha0;

    for (UInt i = 0; i < this->dim_; ++i)
    {
      this->GetThicknessAtPoint(locThick, position, lpm, i);
      if (abs(locThick) > 0.0)
      {
        sigma = sos*this->dampFunction_->ComputeFactor(position, locThick);
        kappa = 1.0 + scalingFunc_->ComputeFactor(position, locThick);
        alpha = shiftFunc_->ComputeFactor(position, locThick);
        frac = sigma/(alpha*alpha + this->omega_*this->omega_);
        Complex fac(kappa + alpha*frac, -1.0*this->omega_*frac);
        scalar *= fac;
      }
      else
      {
        scalar *= dummy;
      }
    }
  }

  //================================================================================================
  // Curvilinear PML
  //================================================================================================

  template<typename T>
  CoefFunctionCurvilinearPML<T>::CoefFunctionCurvilinearPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound, shared_ptr<EntityList> EntList,
                        StdVector<RegionIdType> pdeDomains, OutputType outputType) : CoefFunctionPMLBase<T>(pmlDef, speedOfSound, EntList, pdeDomains) {
    // set name and type
    this->name_ = "CoefFunctionCurvilinearPML";
    this->formulationType_ = CURVILINEAR;
    // get the grid pointer
    grid_ =  this->entities_[0]->GetGrid();
    // set the type of output
    outputType_ = outputType;
    // assign the dimType_ of the coefFunction according to the output type
    SetDimType();
    // read from PML node
    ReadDataPML(pmlDef);
    // check for valid declaration in XML
    CheckForInvalidParams(pmlDef);
    // check if PML layer was generated automatically
    CheckForLayerGenerationNode(pmlDef);
    // trigger computation of geometry data and get pointer to it (is only computed once)
    grid_->GetGeometryOnRegionNodes(ptrNodeGeom_, volRegion_, true);
    // get nodal distances to the PML interface and store in a StdVector
    GetThicknessOnNodes();
    // set identity operators for mapping nodal values to lpm
    CreateMappingOperators();
  }

  template<typename T>
  CoefFunctionCurvilinearPML<T>::~CoefFunctionCurvilinearPML() { };

    template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::TENSOR:
        // compute parameters at lpm
        GetTensorParams(lpm);
        if (this->dim_ == 3) {
          // compute the current wave number
          Double K = this->omega_ / sos_;

          // helper functions
          Vector<Double> h = Vector<Double>(3);
          h[0] = dampFunc_;
          h[1] = intDampFunc_ * kmin_ / (1.0 + dist_ * kmin_);
          h[2] = intDampFunc_ * kmax_ / (1.0 + dist_ * kmax_);
          Double denominator;
          Double K2 = pow(K,2);  //square of wave number

          // get the entries of (I+jD)^-1 with separated real- and imaginary part (is named s but resembles 1/s_i)
          Vector<Complex> s = Vector<Complex>(3);
          for (UInt iDim = 0; iDim < 3; ++iDim) {
            denominator = 1.0 + h[iDim] / K2;
            s[iDim] = Complex(1.0 / denominator, (-h[iDim] / K) / denominator);
          }

          // compute the tensor as the inverse of the Jakobi matrix. The matrix can be factorized in three parts:
          // J^-1 = A^T  * (I + jD)^-1 * A.
          // A is an orthogonal (rotation) matrix containing the curvilinear base vectors (n, t1, t2)
          // I is the identity matrixm, j the imaginary number
          // D is a diagonal matrix containing i.a. the damping functions
          // Inverting (I+jD) is hence simply inverting its entries. Inverting A results in its transpose.
          // Here the matrix is assembled directly in the multiplied version. The result is a symmetric matrix of the form:
          // 
          // / n1^2/s1  + t21^2/s2   + t31^2/s3       |    ...                                    |    ...                           \
          // |
          // | n1*n2/s1 + t21*t22/s2 + t31*t32/s3     |    n2^2/s1  + t22^2/s2   + t32^2/s3       |    ...                           |
          // |
          // \ n1*n3/s1 + t21*t23/s2 + t31*t33/s3     |    n2*n3/s1 + t22*t23/s2 + t32*t33/s3     |    n3^2/s1 + t23^2/s2 + t33^2/s3 /
          tensor.Resize(3, 3);
          tensor[0][0] = pow(n_[0],2)*s[0] + pow(tmin_[0],2)*s[1] + pow(tmax_[0],2)*s[2];
          tensor[1][1] = pow(n_[1],2)*s[0] + pow(tmin_[1],2)*s[1] + pow(tmax_[1],2)*s[2];
          tensor[2][2] = pow(n_[2],2)*s[0] + pow(tmin_[2],2)*s[1] + pow(tmax_[2],2)*s[2];
          tensor[1][0] = n_[0]*n_[1]*s[0] + tmin_[0]*tmin_[1]*s[1] + tmax_[0]*tmax_[1]*s[2];
          tensor[2][0] = n_[0]*n_[2]*s[0] + tmin_[0]*tmin_[2]*s[1] + tmax_[0]*tmax_[2]*s[2];
          tensor[2][1] = n_[1]*n_[2]*s[0] + tmin_[1]*tmin_[2]*s[1] + tmax_[1]*tmax_[2]*s[2];
          tensor[0][1] = tensor[1][0];
          tensor[0][2] = tensor[2][0];
          tensor[1][2] = tensor[2][1];

        } else {
          EXCEPTION("CoefFunctionCurvilinearPML::GetTensor in 2D not implemented yet");
        }
        break;
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetTensor(Complex...) is used for OutputType: " <<
                  "TENSOR only.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalar(Complex& val, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::DETERMINANT:
        // compute parameters at lpm
        GetDeterminantParams(lpm);
        if (this->dim_ == 3)
        {
          // compute the current wave number
          Double K = this->omega_ / sos_;

          // helper functions
          Vector<Double> h = Vector<Double>(3);
          h[0] = dampFunc_;
          h[1] = intDampFunc_ * kmin_ / (1.0 + dist_ * kmin_);
          h[2] = intDampFunc_ * kmax_ / (1.0 + dist_ * kmax_);
          Double K2 = pow(K,2);  //square of wave number

          // compute the determinant...
          // the eigenvalues are of the form s_i = (1 + 1i/K * h[i]). 
          // Hence, the determinant J = s_0*s_1*s_2 computes after multiplication and separation of real/imaginary part:
          val = Complex(1.0 - (h[0]*h[2] + h[1]*h[2] - h[0]*h[1]) / K2, 
                              (h[0]+h[1]+h[2]) / K -  (h[0]*h[1]*h[2])) / pow(K, 3);
        } else {
          EXCEPTION("CoefFunctionCurvilinearPML::GetTensor in 2D not implemented yet");
        }
        break;
      case OutputType::MIN_PRINC_CURV:
        GetMinPrincCurvParams(lpm);
        // convert to complex format
        val = Complex(kmin_, 0);
        break;
      case OutputType::MAX_PRINC_CURV:
        assert(this->dim_ == 3);
        GetMaxPrincCurvParams(lpm);
        // convert to complex format
        val = Complex(kmax_, 0);
        break;
      case OutputType::DISTANCE:
        GetDistParams(lpm);
        // convert to complex format
        val = Complex(dist_, 0);
        break;
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetScalar(Complex...) is used for OutputType: " <<
                  "DETERMINANT, MIN_PRINC_CURV, MAX_PRINC_CURV, and DISTANCE only.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalar(Double& val, const LocPointMapped& lpm ) {
    // distinguish which output type is set and compute / interpolate only required parameters
    switch(outputType_) {
      case OutputType::MIN_PRINC_CURV:
        GetMinPrincCurvParams(lpm);
        val = kmin_;
        break;
      case OutputType::MAX_PRINC_CURV:
        assert(this->dim_ == 3);
        GetMaxPrincCurvParams(lpm);
        val = kmax_;
        break;
      case OutputType::DISTANCE:
        GetDistParams(lpm);
        val = dist_;
        break;
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetScalar(Double...) is used for OutputType: " <<
                  "MIN_PRINC_CURV, MAX_PRINC_CURV, and DISTANCE only.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetVector(Vector<Double>& val, const LocPointMapped& lpm ) {
    // distinguish which output type is set and compute / interpolate only required parameters
    switch(outputType_) {
      case OutputType::NORM_VEC:
        GetNormVecParams(lpm);
        val = n_;
        break;
      case OutputType::MIN_PRINC_VEC:
        GetMinPrincVecParams(lpm);
        val = tmin_;
        break;
      case OutputType::MAX_PRINC_VEC:
        assert(this->dim_ == 3);
        GetMaxPrincVecParams(lpm);
        val = tmax_;
        break;
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetVector(Vector<Double>...) is used for OutputType: " <<
                  "NORM_VEC, MIN_PRINC_VEC, and MAX_PRINC_VEC only.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetVector(Vector<Complex>& val, const LocPointMapped& lpm ) {
    // distinguish which output type is set and compute / interpolate only required parameters
    // resize output vector
    val.Resize(this->dim_);
    switch(outputType_) {
      case OutputType::DAMP_FACTOR:
        // get the required parameters
        GetDampFactorParams(lpm);
        // compute eigenvalues of the PML matrix
        if (this->dim_ == 3) {
          // compute the current wave number
          Double K = this->omega_ / sos_;
          // helper functions
          Vector<Double> h = Vector<Double>(3);
          h[0] = dampFunc_;
          h[1] = intDampFunc_ * kmin_ / (1.0 + dist_ * kmin_);
          h[2] = intDampFunc_ * kmax_ / (1.0 + dist_ * kmax_);

          // get the diagonal entries of (I+jD) with separated real- and imaginary part
          for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
            val[iDim] = Complex(1, h[iDim] / K);
          }
        } else {
          EXCEPTION("CoefFunctionCurvilinearPML::GetVector(Complex...) in 2D not implemented yet");
        }
        break;
      case OutputType::NORM_VEC:
        GetNormVecParams(lpm);
        // convert to complex format
        for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
          val[iDim] = Complex(n_[iDim], 0);
        }
        break;
      case OutputType::MIN_PRINC_VEC:
        GetMinPrincVecParams(lpm);
        // convert to complex format
        for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
          val[iDim] = Complex(tmin_[iDim], 0);
        }
        break;
      case OutputType::MAX_PRINC_VEC:
        assert(this->dim_ == 3);
        GetMaxPrincVecParams(lpm);
        // convert to complex format
        for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
          val[iDim] = Complex(tmax_[iDim], 0);
        }
        break;
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetVector(Vector<Complex>...) is used for OutputType: " <<
                  "DAMP_FACTOR, NORM_VEC, MIN_PRINC_VEC, and MAX_PRINC_VEC only.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetTensorParams(const LocPointMapped& lpm) {
    // here we need all quantites...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // index for vectorial quantities on nodes stored in Vector
    UInt vecIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    // in 2D we only use the tmin_ as tangential vector and kmin_ as curvature
    if (this->dim_ == 2) {
      Vector<Double> normVec(numElemNodes * this->dim_);
      Vector<Double> tangVec(numElemNodes * this->dim_);
      Vector<Double> curv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        // loop over vector entries and store in stacked vector
        for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
          vecIdx = iDim + this->dim_ * iNodes;
          normVec[vecIdx] = ptrNodeGeom_->normalVectors_[nodeIdx][iDim];
          tangVec[vecIdx] = ptrNodeGeom_->minPrincipalVectors_[nodeIdx][iDim];
        }
        curv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> k(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->vectorMappingOperator_->ApplyOp(n_,lpm,ptrFe,normVec);
      this->vectorMappingOperator_->ApplyOp(tmin_,lpm,ptrFe,tangVec);
      this->scalarMappingOperator_->ApplyOp(k,lpm,ptrFe,curv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = k[0];
      dist_ = d[0];
    } // dim_ == 2

    else if (this->dim_ == 3) {
      // vectors with extracted geometry data
      Vector<Double> normVec(numElemNodes * this->dim_);
      Vector<Double> minPrincVec(numElemNodes * this->dim_);
      Vector<Double> maxPrincVec(numElemNodes * this->dim_);
      Vector<Double> minPrincCurv(numElemNodes);
      Vector<Double> maxPrincCurv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        // loop over vector entries and store in stacked vector
        for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
          vecIdx = iDim + this->dim_ * iNodes;
          normVec[vecIdx] = ptrNodeGeom_->normalVectors_[nodeIdx][iDim];
          minPrincVec[vecIdx] = ptrNodeGeom_->minPrincipalVectors_[nodeIdx][iDim];
          maxPrincVec[vecIdx] = ptrNodeGeom_->maxPrincipalVectors_[nodeIdx][iDim];
        }
        minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> kmin(1);
      Vector<Double> kmax(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->vectorMappingOperator_->ApplyOp(n_,lpm,ptrFe,normVec);
      this->vectorMappingOperator_->ApplyOp(tmin_,lpm,ptrFe,minPrincVec);
      this->vectorMappingOperator_->ApplyOp(tmax_,lpm,ptrFe,maxPrincVec);
      this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
      this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = kmin[0];
      kmax_ = kmax[0];
      dist_ = dist[0];
    } // dim_ == 3
    // get the damping function and its integral at the mapped distance
    dampFunc_ = this->dampFunction_->ComputeFactor(dist_, layerThickness_);
    intDampFunc_ = this->dampFunction_->ComputeIntegralFactor(dist_, layerThickness_);

    // finally, get the speed of sound at current lpm
    this->speedOfSound_->GetScalar(sos_,lpm);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetDeterminantParams(const LocPointMapped& lpm) {
    // here we only need curvarute and distances...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    // in 2D we only use the tmin_ as tangential vector and kmin_ as curvature
    if (this->dim_ == 2) {
      // vectors with extracted geometry data
      Vector<Double> curv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        curv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> k(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->scalarMappingOperator_->ApplyOp(k,lpm,ptrFe,curv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = k[0];
      dist_ = d[0];
    } // dim_ == 2
    else if (this->dim_ == 3) {
      // vectors with extracted geometry data
      Vector<Double> minPrincCurv(numElemNodes);
      Vector<Double> maxPrincCurv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> kmin(1);
      Vector<Double> kmax(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
      this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = kmin[0];
      kmax_ = kmax[0];
      dist_ = dist[0];
    } // dim_ == 3
    // get the damping function and its integral at the mapped distance
    dampFunc_ = this->dampFunction_->ComputeFactor(dist_, layerThickness_);
    intDampFunc_ = this->dampFunction_->ComputeIntegralFactor(dist_, layerThickness_);

    // finally, get the speed of sound at current lpm
    this->speedOfSound_->GetScalar(sos_,lpm);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetDampFactorParams(const LocPointMapped& lpm) {
    // here we don't need to interpolate the geometry vectors...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();
    
    // in 2D we only use the kmin_ as curvature
    if (this->dim_ == 2) {
      Vector<Double> curv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        curv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> k(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->scalarMappingOperator_->ApplyOp(k,lpm,ptrFe,curv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = k[0];
      dist_ = d[0];
    } // dim_ == 2

    else if (this->dim_ == 3) {
      // vectors with extracted geometry data
      Vector<Double> minPrincCurv(numElemNodes);
      Vector<Double> maxPrincCurv(numElemNodes);
      // and on node-to-interface distances
      Vector<Double> dist(numElemNodes);

      // loop over element nodes and get quantities
      for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
        nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
        minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
        maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx];
        dist[iNodes] = thicknessOnNodes_[nodeIdx];
      }

      // create helper variables for scalars because ApplyOp only takes Vectors
      Vector<Double> kmin(1);
      Vector<Double> kmax(1);
      Vector<Double> d(1);

      // interpolate from nodes to the local point mapped
      this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
      this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
      this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
      // extract from vector
      kmin_ = kmin[0];
      kmax_ = kmax[0];
      dist_ = dist[0];
    } // dim_ == 3
    // get the damping function and its integral at the mapped distance
    dampFunc_ = this->dampFunction_->ComputeFactor(dist_, layerThickness_);
    intDampFunc_ = this->dampFunction_->ComputeIntegralFactor(dist_, layerThickness_);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetNormVecParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the normal vectors...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // index for vectorial quantities on nodes stored in Vector
    UInt vecIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();
    Vector<Double> normVec(numElemNodes * this->dim_);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      // loop over vector entries and store in stacked vector
      for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
        vecIdx = iDim + this->dim_ * iNodes;
        normVec[vecIdx] = ptrNodeGeom_->normalVectors_[nodeIdx][iDim];
      }
    }
    // interpolate from nodes to the local point mapped
    this->vectorMappingOperator_->ApplyOp(n_,lpm,ptrFe,normVec);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetMinPrincVecParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the minimum principal vectors...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // index for vectorial quantities on nodes stored in Vector
    UInt vecIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    Vector<Double> minPrincVec(numElemNodes * this->dim_);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      // loop over vector entries and store in stacked vector
      for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
        vecIdx = iDim + this->dim_ * iNodes;
        minPrincVec[vecIdx] = ptrNodeGeom_->minPrincipalVectors_[nodeIdx][iDim];
      }
    }
    // interpolate from nodes to the local point mapped
    this->vectorMappingOperator_->ApplyOp(tmin_,lpm,ptrFe,minPrincVec);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetMinPrincCurvParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the minimum principal curvatures...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    Vector<Double> curv(numElemNodes);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      curv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
    }
    // create helper variables for scalars because ApplyOp only takes Vectors
    Vector<Double> k(1);
    // interpolate from nodes to the local point mapped
    this->scalarMappingOperator_->ApplyOp(k,lpm,ptrFe,curv);
    // extract from vector
    kmin_ = k[0];
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetMaxPrincVecParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the maximum principal vectors...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // index for vectorial quantities on nodes stored in Vector
    UInt vecIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    Vector<Double> maxPrincVec(numElemNodes * this->dim_);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      // loop over vector entries and store in stacked vector
      for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
        vecIdx = iDim + this->dim_ * iNodes;
        maxPrincVec[vecIdx] = ptrNodeGeom_->maxPrincipalVectors_[nodeIdx][iDim];
      }
    }
    // interpolate from nodes to the local point mapped
    this->vectorMappingOperator_->ApplyOp(tmax_,lpm,ptrFe,maxPrincVec);
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetMaxPrincCurvParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the maximum principal curvatures...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    Vector<Double> curv(numElemNodes);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      curv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx];
    }
    // create helper variables for scalars because ApplyOp only takes Vectors
    Vector<Double> k(1);
    // interpolate from nodes to the local point mapped
    this->scalarMappingOperator_->ApplyOp(k,lpm,ptrFe,curv);
    // extract from vector
    kmax_ = k[0];
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetDistParams(const LocPointMapped& lpm) {
    // here we need to interpolate only the nodal distances to the interface...
    // get a pointer to the considered element from the lpm
    const Elem* ptrElem = NULL;
    ptrElem = lpm.ptEl;
    // the global element ID
    UInt elemId = ptrElem->GetElemNum();
    // get the nodes connected to the element
    StdVector<UInt> nodeIds;
    grid_->GetElemNodes(nodeIds, elemId);
    // index of a node id
    UInt nodeIdx;
    // number of nodes in element
    UInt numElemNodes = nodeIds.GetSize();

    // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
    shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

    // get Base Fe, which provides the shape functions for interpolating with the identity operator 
    BaseFE* ptrFe = ptrEsm->GetBaseFE();

    Vector<Double> dist(numElemNodes);
    // loop over element nodes and get quantities
    for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
      nodeIdx = GetIdxByNodeId(nodeIds[iNodes]);
      dist[iNodes] = thicknessOnNodes_[nodeIdx];
    }
    // create helper variables for scalars because ApplyOp only takes Vectors
    Vector<Double> d(1);
    // interpolate from nodes to the local point mapped
    this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
    // extract from vector
    dist_ = d[0];
  };

  template<typename T>
  UInt CoefFunctionCurvilinearPML<T>::GetIdxByNodeId(const UInt& nodeId) const {
    UInt idx = 0;
    if (layerGenNode_) {
      UInt tmpIdx, tmpId;
      // assume ordered nodes and use the knowledge of the layer generation params to obtain the index
      // loop over surface regions in layer and check if id lies on current surface
      for (UInt iLayers = 1; iLayers <= numLayers_; ++iLayers) {
        tmpIdx = iLayers * numSurfNodes_ -1; //last index on one surface
        tmpId = ptrNodeGeom_->nodeIds_[tmpIdx];
        // if the id lies on the current surface, set ID
        if (nodeId <= tmpId) {
          idx = tmpIdx - tmpId + nodeId; // set index
          // check if search found the correct value 
          // (only the case if nodes are sorted and nodeIDs of a surface have no 'holes')
          // use conventional search otherwise
          idx = (idx >= numLayers_ * numSurfNodes_) ? 0 : idx; // if idx is too high set to 0
          if (ptrNodeGeom_->nodeIds_[idx] != nodeId) { // check if idx is correct
            unsigned int guess = tmpIdx - numSurfNodes_ + 1;
            int tmpIdx2 = ptrNodeGeom_->nodeIds_.Find(nodeId, guess, false);
            idx = tmpIdx2;
          }
          iLayers = numLayers_+1; // break the loop because idx was found
        }
      }
    } else {
      // should be made more efficient in future!!
      unsigned int guess = 0;
      int tmpIdx2 = ptrNodeGeom_->nodeIds_.Find(nodeId, guess, false);
      idx = tmpIdx2;
    }
    return idx;
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetThicknessOnNodes() {
    if (layerGenNode_) {
      // every iso surface in the automatically generated layer has the 
      // same number of nodes. Hence, finding the thickness at nodes is easy...
      // iterate over all surface layers (generated and the interface) and assign
      for (UInt iLayers = 0; iLayers <= numLayers_; iLayers++) {
        for (UInt iNodes = 0; iNodes < numSurfNodes_; iNodes++) {
          thicknessOnNodes_.Push_back(elemHeight_ * iLayers);
        }
      }
    } else {
      // room for future implementation that e.g. reads data from a file
      EXCEPTION("CoefFunctionCurvilinearPML::GetThicknessOnNodes currently only implemented with auto layer generation specified.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CreateMappingOperators(){
    // create identity operators for 2D or 3D problems
    if (this->dim_ == 2) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,2,1,Double>());
      // in the tensor or vector case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::TENSOR || this->dimType_ == CoefFunction::VECTOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,2,3,Double>());
    }
    else if (this->dim_ == 3) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,3,1,Double>());
      // in the tensor or vector case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::TENSOR || this->dimType_ == CoefFunction::VECTOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,3,3,Double>());
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::SetDimType() {
    switch (outputType_) {
      case OutputType::TENSOR:
        this->dimType_ = CoefFunction::TENSOR;
        break;
      case OutputType::DETERMINANT:
        this->dimType_ = CoefFunction::SCALAR;
        break;
      case OutputType::DAMP_FACTOR:
        this->dimType_ = CoefFunction::VECTOR;
        break;
      case OutputType::NORM_VEC:
        this->dimType_ = CoefFunction::VECTOR;
        break;
      case OutputType::MIN_PRINC_VEC:
        this->dimType_ = CoefFunction::VECTOR;
        break;
      case OutputType::MIN_PRINC_CURV:
        this->dimType_ = CoefFunction::SCALAR;
        break;
      case OutputType::MAX_PRINC_VEC:
        // not used in 2D
        if (this->dim_ != 3) {
          EXCEPTION("In CoefFunctionCurvilinearPML: outputType '" << outputType_ <<
                    "' is only valid for 3D problems.");
        } else {
          this->dimType_ = CoefFunction::VECTOR;
          break;
        }
      case OutputType::MAX_PRINC_CURV:
        // not used in 2D
        if (this->dim_ != 3) {
          EXCEPTION("In CoefFunctionCurvilinearPML: outputType '" << outputType_ <<
                    "' is only valid for 3D problems.");
        } else {
          this->dimType_ = CoefFunction::SCALAR;
          break;
        }
      case OutputType::DISTANCE:
        this->dimType_ = CoefFunction::SCALAR;
        break;
      default:
        EXCEPTION("In CoefFunctionCurvilinearPML: outputType '" << outputType_ <<
                        "' is invalid.");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::ReadDataPML(PtrParamNode pmlDef){
    // read and set coordinate system 
    std::string cSysId;
    pmlDef->GetValue("coordSysId",cSysId,ParamNode::PASS);
    // currently this type wil only be implemented to work in Cartesian coordinates
    // to work with other coordinate systems we need to compute different transformation matrices
    if(cSysId == "" || cSysId == "default")
      cSysId = "default";
    else
      EXCEPTION("Curvilinear PML currently only works in the default Cartesian coordinate system!");
    this->SetCoordinateSystem(domain->GetCoordSystem(cSysId));

    // type of PML damping function
    std::string typeOfPml;
    pmlDef->GetValue("type",typeOfPml);
    this->pmlType_ = DampFunction::DampingTypeEnum.Parse(typeOfPml);
    // set the DampingFunction object to the corresponding type 
    CreateDampFunction();

    // check for damping factor
    Double dampFactor;
    pmlDef->GetValue("dampFactor",dampFactor);
    this->dampFunction_->DampFactor = dampFactor;

    // extract the ID of the PML volume region
    volRegion_ = this->entities_[0]->GetRegion();
  }

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CheckForInvalidParams(PtrParamNode pmlDef) {
    // check for propRegion, scaling or frequency shift coeff in the xml
    // if these are set for curvilinear PML, ignore and warn
    PtrParamNode propRegionNode = pmlDef->Get("propRegion", ParamNode::PASS);
    if (propRegionNode) {
      std::string propRegionNodeFormul; 
      propRegionNode->GetParent()->GetValue("formulation", propRegionNodeFormul, ParamNode::PASS);
      if (propRegionNodeFormul == "curvilinear")
        WARN("propRegion is currently not implemented for curvilinear PML and will be ignored!");
    }
    PtrParamNode scalingNode = pmlDef->Get("scalingCoef", ParamNode::PASS);
    if (scalingNode) {
      std::string scalingNodeFormul; 
      scalingNode->GetParent()->GetValue("formulation", scalingNodeFormul, ParamNode::PASS);
      if (scalingNodeFormul == "curvilinear")
        WARN("scalingCoef is currently not implemented for curvilinear PML and will be ignored!");
    }
    PtrParamNode shiftNode = pmlDef->Get("frqShiftCoef", ParamNode::PASS);
    if (shiftNode) {
      std::string shiftNodeFormul; 
      shiftNode->GetParent()->GetValue("formulation", shiftNodeFormul, ParamNode::PASS);
      if (shiftNodeFormul == "curvilinear")
        WARN("frqShiftCoef is currently not implemented for curvilinear PML and will be ignored!");
    }
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CheckForLayerGenerationNode(PtrParamNode pmlDef) {
    // check wether the PML layer has been generated using autoLayerGenerationList and store the parameter node
    PtrParamNode layerGenNode = pmlDef->GetParent()->GetParent()->GetParent()->
                                        GetParent()->GetParent()->Get("domain")->
                                        Get("autoLayerGenerationList", ParamNode::PASS);
    // check if the autoLayerGenerationList has actually generated the PML region
    if (layerGenNode) {
      std::string generatedLayerName;
      RegionIdType generatedLayerId;
      PtrParamNode actLayerGenNode;
      // iterate over the autoLayerGenerationList...
      UInt numGenLayers = layerGenNode->GetChildren().size();
      bool nodeFound = false;
      for (UInt iLayers = 0; iLayers < numGenLayers; iLayers++) {
        PtrParamNode actLayerGenNode = layerGenNode->GetChildren()[iLayers];
        // get name of generated region
        actLayerGenNode->GetValue("name", generatedLayerName);
        // get its ID
        generatedLayerId = grid_->GetRegionId(generatedLayerName);
        // check if generated ID fits to the PML region
        if (volRegion_ == generatedLayerId) {
          layerGenNode_ = actLayerGenNode;

          // assign layer generation parameters for later use
          std::string surfRegionName;
          RegionIdType surfRegionId;
          layerGenNode_->GetValue("numLayers", numLayers_);
          layerGenNode_->GetValue("elemHeight", elemHeight_);
          layerGenNode_->GetValue("surfRegionToActOn", surfRegionName);
          surfRegionId = grid_->GetRegionId(surfRegionName);

          // set the total layer thickness
          layerThickness_ = numLayers_ * elemHeight_;
          // obtain the number of nodes on the interface region
          numSurfNodes_ = grid_->GetNumNodes(surfRegionId);
          // set bool to succeed search
          nodeFound = true;
        }
      }
      if (nodeFound == false) {
        EXCEPTION("Specified PML layer not found in autoLayerGenerationList");
      }
    }
  };

  // Explicit template instantiation
  template class CoefFunctionPMLBase<Double>;
  template class CoefFunctionPMLBase<Complex>;
  template class CoefFunctionPML<Double>;
  template class CoefFunctionPML<Complex>;
  template class CoefFunctionShiftedPML<Double>;
  template class CoefFunctionShiftedPML<Complex>;
  template class CoefFunctionCurvilinearPML<Double>;
  template class CoefFunctionCurvilinearPML<Complex>;
}
