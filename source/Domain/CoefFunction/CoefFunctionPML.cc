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
                        StdVector<RegionIdType> pdeDomains, bool isTensor) : CoefFunctionPMLBase<T>(pmlDef, speedOfSound, EntList, pdeDomains) {
    // set name and type
    this->name_ = "CoefFunctionCurvilinearPML";
    this->formulationType_ = CURVILINEAR;
    
    // get the grid pointer
    grid_ =  this->entities_[0]->GetGrid();
    
    // set if the CoefFunction is used as scalar or tensor
    if (isTensor == true)
      this->dimType_ = CoefFunction::TENSOR;
    else
      this->dimType_ = CoefFunction::SCALAR;

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
  CoefFunctionCurvilinearPML<T>::~CoefFunctionCurvilinearPML() { }

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetParamsAtLocalPoint(const LocPointMapped& lpm) {
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

    // distinguish how many params we require...
    // when we compute a tensor we require for the whole geometry
    if (this->dimType_ == CoefFunction::TENSOR) {
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
    } // dimType_ == CoefFunction::TENSOR

    // for the scalar (=Jakobi determinant), we only need of the some parameters
    else if (this->dimType_ == CoefFunction::SCALAR) {
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
    } // dimType_ == CoefFunction::SCALAR

    // get the damping function and its integral at the mapped distance
    dampFunc_ = this->dampFunction_->ComputeFactor(dist_, layerThickness_);
    intDampFunc_ = this->dampFunction_->ComputeIntegralFactor(dist_, layerThickness_);

    // finally, get the speed of sound at current lpm
    this->speedOfSound_->GetScalar(sos_,lpm);
  }

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) {
    if (this->dim_ == 3) {
      // interpolate nodal values to local point
      GetParamsAtLocalPoint(lpm);

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
      EXCEPTION("");
    }
    
  };

  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalar(Complex& val, const LocPointMapped& lpm ) {
    if (this->dim_ == 3)
    {
      // interpolate nodal values to local point
      GetParamsAtLocalPoint(lpm);

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
      EXCEPTION("");
    }
  };




  template<typename T>
  void CoefFunctionCurvilinearPML<T>::ComputeTensorsAndScalarsOnNodes() {
    Double sos = 343;

    // declare variables...
    // number of nodes in layer
    UInt numNodes = ptrNodeGeom_->numNodes_;
    // damp function evaluated at current node and integral of the damp function from the interface to the node
    Double dampFunc, intDampFunc;
    // the orthogonal part of the Jakobi matrix that holds the curvilinear base vectors
    Matrix<Complex> orthoMat = Matrix<Complex>(3, 3);
    Matrix<Complex> orthoMatT = Matrix<Complex>(3,3); // transposed
    // the diagonal part that does the (complex) metric stretching
    Matrix<Complex> stretchMat = Matrix<Complex>(3, 3);
    // set the helper function that does the damping
    StdVector<Double> dampingPart = StdVector<Double>(3);
    // another helper variable
    Double denominator;

    // resize the StdVectors
    scalarsOnNodes_.Resize(numNodes);
    tensorsOnNodes_.Resize(numNodes);

    // compute for every node in ptrNodeGeom_
    for (UInt iNodes = 0; iNodes < numNodes; iNodes++) {
      dampFunc = this->dampFunction_->ComputeFactor(thicknessOnNodes_[iNodes], layerThickness_);
      intDampFunc = this->dampFunction_->ComputeIntegralFactor(thicknessOnNodes_[iNodes], layerThickness_);

      dampingPart[0] = dampFunc;
      dampingPart[1] = intDampFunc * ptrNodeGeom_->minPrincipalCurvatures_[iNodes] / (1 + thicknessOnNodes_[iNodes] * ptrNodeGeom_->minPrincipalCurvatures_[iNodes]);
      dampingPart[2] = intDampFunc * ptrNodeGeom_->maxPrincipalCurvatures_[iNodes] / (1 + thicknessOnNodes_[iNodes] * ptrNodeGeom_->maxPrincipalCurvatures_[iNodes]);

      // compute and store the factorized Jakobi matrix...
      // fill the matrices with entries
      for (UInt iColumn = 0; iColumn < 3; iColumn++) {
        orthoMat[0][iColumn] = Complex(ptrNodeGeom_->normalVectors_[iNodes][iColumn], 0.0); // imaginary part is 0
        orthoMat[1][iColumn] = Complex(ptrNodeGeom_->minPrincipalVectors_[iNodes][iColumn], 0.0);
        orthoMat[2][iColumn] = Complex(ptrNodeGeom_->maxPrincipalVectors_[iNodes][iColumn], 0.0);

        denominator = 1 + pow(dampingPart[iColumn]*sos / this->omega_, 2);
        stretchMat[iColumn][iColumn] = Complex(1 / denominator, -1 * (dampingPart[iColumn]*sos / this->omega_) / denominator);
      }
      orthoMat.Transpose(orthoMatT);
      tensorsOnNodes_[iNodes] = orthoMatT * stretchMat * orthoMat;

      // compute and store the determinant...
      // the eigenvalues are of the form s_i = (1 + 1i/K * dampingPart[i]). 
      // Hence, the determinant J = s_0*s_1*s_2 computes after multiplication and separation of real/imaginary part:
      scalarsOnNodes_[iNodes] = Complex(1 - pow(sos/this->omega_, 2) * (dampingPart[0]*dampingPart[2] + dampingPart[1]*dampingPart[2] - dampingPart[0]*dampingPart[1]), 
                      sos/this->omega_ * (dampingPart[0]+dampingPart[1]+dampingPart[2]) - pow(sos/this->omega_, 3) * (dampingPart[0]*dampingPart[1]*dampingPart[2]));
    }
  }


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetTensorOnNode(Matrix<Complex>& tensor, const UInt& nodeId) {
        Double sos = 343;
    // get the index of the node
    UInt nodeIdx = GetIdxByNodeId(nodeId);

    // get the distance of the node to the interface
    Double dist = thicknessOnNodes_[nodeId];

    // get the damp function evaluated at current node
    Double dampFunc = this->dampFunction_->ComputeFactor(dist, layerThickness_);
    // and its integral from the interface to the node
    Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(dist, layerThickness_);

    // the orthogonal part of the Jakobi matrix that holds the curvilinear base vectors
    Matrix<Complex> orthoMat = Matrix<Complex>(3, 3);
    Matrix<Complex> orthoMatT = Matrix<Complex>(3,3); // transposed
    
    // the diagonal part that does the metric stretching
    Matrix<Complex> stretchMat = Matrix<Complex>(3, 3);

    // set the helper function that does the damping
    StdVector<Double> dampingPart = StdVector<Double>(3);
    dampingPart[0] = dampFunc;
    dampingPart[1] = intDampFunc * ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx] / (1 + dist * ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx]);
    dampingPart[2] = intDampFunc * ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx] / (1 + dist * ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx]);
    // another helper variable
    Double denominator;
    
    // fill the matrices with entries
    for (UInt iColumn = 0; iColumn < 3; iColumn++) {
      orthoMat[0][iColumn] = Complex(ptrNodeGeom_->normalVectors_[nodeIdx][iColumn], 0.0); // imaginary part is 0
      orthoMat[1][iColumn] = Complex(ptrNodeGeom_->minPrincipalVectors_[nodeIdx][iColumn], 0.0);
      orthoMat[2][iColumn] = Complex(ptrNodeGeom_->maxPrincipalVectors_[nodeIdx][iColumn], 0.0);

      denominator = 1 + pow(dampingPart[iColumn]*sos / this->omega_, 2);
        stretchMat[iColumn][iColumn] = Complex(1 / denominator, -1 * (dampingPart[iColumn]*sos / this->omega_) / denominator);
    }
    orthoMat.Transpose(orthoMatT);

    tensor = orthoMatT * stretchMat * orthoMat;
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalarOnNode(Complex& scalar, const UInt& nodeId) {
        Double sos = 343;
    // get the index of the node
    UInt nodeIdx = GetIdxByNodeId(nodeId);

    // get the distance of the node to the interface
    Double dist = thicknessOnNodes_[nodeId];

    // get the damp function evaluated at current node
    Double dampFunc = this->dampFunction_->ComputeFactor(dist, layerThickness_);
    // and its integral from the interface to the node
    Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(dist, layerThickness_);

    // set the helper function that does the damping
    StdVector<Double> dampingPart = StdVector<Double>(3);
    dampingPart[0] = dampFunc;
    dampingPart[1] = intDampFunc * ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx] / (1 + dist * ptrNodeGeom_->minPrincipalCurvatures_[nodeIdx]);
    dampingPart[2] = intDampFunc * ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx] / (1 + dist * ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdx]);
    
    // the eigenvalues are of the form s_i = (1 + 1i/K * dampingPart[i]). Hence, the determinant J = s_0*s_1*s_2 computes after multiplication and separation of real/imaginary part:
    scalar = Complex(1 - pow(sos/this->omega_, 2) * (dampingPart[0]*dampingPart[2] + dampingPart[1]*dampingPart[2] - dampingPart[0]*dampingPart[1]), 
                      sos/this->omega_ * (dampingPart[0]+dampingPart[1]+dampingPart[2]) - pow(sos/this->omega_, 3) * (dampingPart[0]*dampingPart[1]*dampingPart[2]));
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
      EXCEPTION("");
      // should be made more efficient in future!!
      //int idx2 = ptrNodeGeom_->nodeIds_.Find(nodeId);
      //idx = idx2;
    }
    return idx;
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
    } else {
      EXCEPTION("Curvilinear PML currently only works with an autonatically generated layer in autoLayerGenerationList");
    }
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
      EXCEPTION("Curvilinear PML currently only works with an automatically generated layer specified in autoLayerGenerationList");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CreateMappingOperators(){
    // create identity operators for 2D or 3D problems
    if (this->dim_ == 2) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,2,1,Double>());
      // in the tensor case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::TENSOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,2,3,Double>());
    }
    else if (this->dim_ == 3) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,3,1,Double>());
      // in the tensor case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::TENSOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,3,3,Double>());
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
