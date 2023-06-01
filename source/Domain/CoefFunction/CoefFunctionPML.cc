// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionPML.hh
 *       \brief    This file contains the CoefFunctionPMLBase class, the classic/Cartesian 
 *                 CoefFunctionPML and the (also Cartesian) CoefFunctionShiftedPML
 *
 *       \date     Mar 12, 2012, edit by pheidegger in Mar 2023
 *       \author   ahueppe, pheidegger
 */
//================================================================================================

#include "CoefFunctionPML.hh"
namespace CoupledField{
  //================================================================================================
  // Base Class
  //================================================================================================
  template<typename T>
  CoefFunctionPMLBase<T>::CoefFunctionPMLBase(PtrParamNode pmlDef, PtrCoefFct matCoef,
                                              shared_ptr<EntityList> EntList, 
                                              StdVector<RegionIdType> pdeDomains,
                                              CoefFunction::CoefDimType dimType): 
                                    CoefFunction() {

    isAnalytic_ = false;
    isComplex_ =  std::is_same<T,Complex>::value;
    dependType_ = GENERAL;
    this->entities_.Push_back(EntList);
    dim_ = entities_[0]->GetGrid()->GetDim();

    //set type: SCALAR, VECTOR, TENSOR :
    this->dimType_ = dimType;

    //orderCoefFct_ = orderCoefFct;
    name_ = "CoefFunctionPMLBase";
    pmlType_ = DampFunction::NO_TYPE;

    //check if the coefficient function of the speed of sound is real valued and scalar
    if(matCoef->GetDimType() != CoefFunction::SCALAR || matCoef->IsComplex()){
      EXCEPTION("The PML coefficient function needs a real valued scalar coeffunction as speed of wave");
    }

    matCoef_ = matCoef;
  
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
  CoefFunctionPML<T>::CoefFunctionPML(PtrParamNode pmlDef, PtrCoefFct matCoef,
              shared_ptr<EntityList> EntList, StdVector<RegionIdType> pdeDomains, CoefFunction::coefDimType dimType) : 
                  CoefFunctionPMLBase<T>(pmlDef, matCoef, EntList, pdeDomains, dimType) {
                    
    this->name_ = "CoefFunctionPML";
    //prepare dimenstions of propagation region
    innerMinMaxComp_.Resize(3,2);
    innerMinMaxComp_.Init();

    outerMinMaxComp_.Resize(3,2);
    outerMinMaxComp_.Init();

    this->formulationType_ = CLASSIC;
    ReadDataPML(pmlDef,pdeDomains);

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
    Double waveSpeed;
    this->matCoef_->GetScalar(waveSpeed,lpm);
    
    if(this->used4PDE_ == "fullwave-E"){
      //special tensor for Electromagnetic Waves
      //see book FINITE ELEMENT ANALYSIS OF ANTENNAS AND ARRAYS, Jian-Ming Jin, Douglas J. Riley, pp. 64
      //however, we scale with wave speed!
      for(UInt i=0;i<this->dim_;++i) {
        GetThicknessAtPoint(locThick,position,lpm,i);
        if(abs(locThick)>0.0){
          Complex fac(0.0, this->dampFunction_->ComputeFactor(position,locThick));
          sFactors[i] =  1.0 - fac * waveSpeed / (this->omega_);
        }else{
          sFactors[i] = one;
        }
      } 
      
      // removed isActive_ as it was always false!
      tensor[0][0] = sFactors[1]*sFactors[2] / sFactors[0];
      tensor[1][1] = sFactors[0]*sFactors[2] / sFactors[1];
      tensor[2][2] = sFactors[0]*sFactors[1] / sFactors[2];
    }
    else{
      for(UInt i=0;i<this->dim_;++i){
        GetThicknessAtPoint(locThick,position,lpm,i);
        if(abs(locThick)>0.0){
          Complex fac(0.0, -1.0 * waveSpeed * this->dampFunction_->ComputeFactor(position,locThick));
          tensor[i][i] = this->omega_ / (this->omega_ + fac);
        }else{
          tensor[i][i] = one;
        }
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
    Double waveSpeed;
    this->matCoef_->GetScalar(waveSpeed,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        tensor[i][i] = this->dampFunction_->ComputeFactor(position,locThick) * waveSpeed;
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
    Double waveSpeed;
    this->matCoef_->GetScalar(waveSpeed,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        Complex fac(1.0, -1.0 * waveSpeed * this->dampFunction_->ComputeFactor(position,locThick)/this->omega_);
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
    Double waveSpeed;
    this->matCoef_->GetScalar(waveSpeed,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        vec[i] = this->dampFunction_->ComputeFactor(position,locThick) * waveSpeed;
      }else{
        vec[i] = 0.0;
      }
    }
  }

template<typename T>
void CoefFunctionPML<T>::UpdateOmega(){
  omega_ = this->mp_->Eval(mHandle_) * 2 * M_PI;
}

  template<typename T>
  void CoefFunctionPML<T>::GetScalar(Complex& val,
                                const LocPointMapped& lpm ){
    //computes (1-j pml_x/K) *(1-j pml_y/K)*(1-j pml_z/K))
    Double locThick=0.0;
    Double position=0.0;
    Complex dummy(1.0,0.0);
    val = dummy;
    Double waveSpeed;
    this->matCoef_->GetScalar(waveSpeed,lpm);
    for(UInt i=0;i<this->dim_;++i){
      GetThicknessAtPoint(locThick,position,lpm,i);
      if(abs(locThick)>0.0){
        Complex fac(1.0, -1.0 * waveSpeed * this->dampFunction_->ComputeFactor(position,locThick)/this->omega_);
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

template<typename T>
CoefFunctionShiftedPML<T>::CoefFunctionShiftedPML(PtrParamNode pmlDef, PtrCoefFct matCoef, shared_ptr<EntityList> EntList,
                                                  StdVector<RegionIdType> pdeDomains, CoefFunction::CoefDimType dimType)
  : CoefFunctionPML<T>(pmlDef, matCoef, EntList, pdeDomains, dimType)
{
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
    Double waveSpeed, alpha, kappa, sigma, frac;
    this->matCoef_->GetScalar(waveSpeed, lpm);

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
        sigma = waveSpeed*this->dampFunction_->ComputeFactor(position, locThick);
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
    Double waveSpeed, alpha, kappa, sigma, frac;
    this->matCoef_->GetScalar(waveSpeed, lpm);

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
        sigma = waveSpeed*this->dampFunction_->ComputeFactor(position, locThick);
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
  // Explicit template instantiation
  template class CoefFunctionPMLBase<Double>;
  template class CoefFunctionPMLBase<Complex>;
  template class CoefFunctionPML<Double>;
  template class CoefFunctionPML<Complex>;
  template class CoefFunctionShiftedPML<Double>;
  template class CoefFunctionShiftedPML<Complex>;
}
