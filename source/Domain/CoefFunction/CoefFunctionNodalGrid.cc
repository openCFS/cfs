// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionNodalGrid.cc
 *       \brief    Implementation file for the nodal Grid interpolation CoefFunction
 *
 *       \date     Jun 27, 2012
 *       \author   ahueppe
 */
//================================================================================================

#include "CoefFunctionNodalGrid.hh"


#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "DataInOut/ResultHandler.hh"

#include <boost/lexical_cast.hpp>
#include <boost/tr1/type_traits.hpp>

namespace CoupledField{

  template<class DATA_TYPE>
  CoefFunctionNodalGrid<DATA_TYPE>::CoefFunctionNodalGrid(PtrParamNode configNode){

    this->ReadXmlNode(configNode);

    UInt dofDim = 0;
     if(resultInfo_->entryType == ResultInfo::SCALAR){
       dofDim = 1;
       dimType_ = SCALAR;
     }else if(resultInfo_->entryType == ResultInfo::VECTOR){
       dofDim = resultInfo_->dofNames.GetSize();
       dimType_ = VECTOR;
     }else if(resultInfo_->entryType == ResultInfo::TENSOR){
       dofDim = resultInfo_->dofNames.GetSize();
       dimType_ = TENSOR;
     }
     Grid* curGrid = domain->GetGrid(gridId_);
     UInt dim = curGrid->GetDim();
     //right now we hardcode this here as we have always nodal results
     //one can think about setting this operator externally but we come to this later
     if(dofDim == 1){
       if(dim == 2)
         this->myOperator_.reset(new IdentityOperator<FeH1,2,1,DATA_TYPE>());
       else if(dim == 3)
         this->myOperator_.reset(new IdentityOperator<FeH1,3,1,DATA_TYPE>());
     }else if(dofDim == 2){
       if(dim == 2)
         this->myOperator_.reset(new IdentityOperator<FeH1,2,2,DATA_TYPE>());
       else if(dim == 3)
         this->myOperator_.reset(new IdentityOperator<FeH1,3,2,DATA_TYPE>());
     }else if(dofDim == 3){
       if(dim == 2)
         this->myOperator_.reset(new IdentityOperator<FeH1,2,3,DATA_TYPE>());
       else if(dim == 3)
         this->myOperator_.reset(new IdentityOperator<FeH1,3,3,DATA_TYPE>());
     }else{
       EXCEPTION("Dimension of result not supported. Is it tensorial???")
     }

    // obtain handle from internal variable coefficient function
    this->mp_ = domain->GetMathParser();
    this->mHandleStep_ = this->mp_->GetNewHandle(true);

    // register callback mechanism if expression changes
    //mp_->AddExpChangeCallBack(
    //    boost::bind(&CoefFunctionGridBase<DATA_TYPE>::Recalculate, this ),
    //    mHandleStep_ );

    bool complex = std::tr1::is_same<DATA_TYPE,Complex>::value;
    std::string value;
    if(complex)
      value = "f";
    else
      value = "t";

    this->mp_->SetExpr(this->mHandleStep_, value);

    this->curStep_ = 1;
    this->curTStep_ = 0;
  }

  template<class DATA_TYPE>
  UInt CoefFunctionNodalGrid<DATA_TYPE>::GetVecSize() const {
    return domain->GetGrid(gridId_)->GetDim();
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("Not tensor valued data available yet");
  }

  template<class DATA_TYPE>
  bool CoefFunctionNodalGrid<DATA_TYPE>:: IsComplex(){
    return std::tr1::is_same<DATA_TYPE,Complex>::value;
  }


  template<class DATA_TYPE>
  CoefFunctionNodalGrid<DATA_TYPE>::~CoefFunctionNodalGrid(){

  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                              const LocPointMapped& lpm ){
    EXCEPTION("TENSORIAL RESULTS NOT SUPPORTED FOR INTERPOLATION FROM EXTERNAL GRIDS");

  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                              const LocPointMapped& lpm ){
    //1. call to grid to find the element
    //2. compute solution of this element
    //3. apply operator
    assert(this->dimType_ != CoefFunction::TENSOR);
    if(this->dimType_ == CoefFunction::SCALAR)
      CoefMat.Resize(1);
    else if (this->dimType_ == CoefFunction::VECTOR)
      CoefMat.Resize(domain->GetGrid(this->gridId_)->GetDim());

    CoefMat.Init();
#ifdef USE_CGAL
    //const Elem* targetElem = lpm.ptEl;
    Vector<Double> globCoord;
    Vector<Double> localCoord;
    Vector<DATA_TYPE> elemSol;

    const Elem* sourceElem = NULL;
    LocPoint lp;
    if(gridId_=="default"){
      //ok so computational grid is the same as source grid we obtain directly the element
      //from the grid
      UInt elemenum = lpm.ptEl->elemNum;
      sourceElem = domain->GetGrid(gridId_)->GetElem(elemenum);
      lp = lpm.lp;
    }else{

      lpm.shapeMap->Local2Global(globCoord,lpm.lp);
      sourceElem = domain->GetGrid(gridId_)->
                   GetElemAtGlobalCoord(globCoord,lp,srcRegions_);
    }

    if(!sourceElem){
      WARN("Could not find element for coord " << globCoord);
      return;
    }

    shared_ptr<ElemShapeMap> esm = domain->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );

    LocPointMapped lpmS;
    lpmS.Set(lp,esm);

    this->GetElemSolution( elemSol, sourceElem->elemNum);
    BaseFE * ptFe = esm->GetBaseFE();
    myOperator_->ApplyOp(CoefMat,lpmS,ptFe,elemSol);

    //now obtain factor
    CoefMat *= this->GetFactor(lpmS);

  #else
    EXCEPTION("This Executeable does not support interpolation and thus does not support the inclustion of external results");
  #endif
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                             const LocPointMapped& lpm ){

    assert(this->dimType_ == CoefFunction::SCALAR );

    //const Elem* targetElem = lpm.ptEl;
    Vector<Double> globCoord;
    Vector<Double> localCoord;
    Vector<DATA_TYPE> elemSol;
    Vector<DATA_TYPE> ptSol(1);
    ptSol.Init();

  #ifdef USE_CGAL
    const Elem* sourceElem = NULL;
    LocPoint lp;
    if(gridId_=="default"){
      //ok so computational grid is the same as source grid we obtain directly the element
      //from the grid
      UInt elemenum = lpm.ptEl->elemNum;
      sourceElem = domain->GetGrid(gridId_)->GetElem(elemenum);
      lp = lpm.lp;
    }else{

      lpm.shapeMap->Local2Global(globCoord,lpm.lp);
      sourceElem = domain->GetGrid(gridId_)->
                   GetElemAtGlobalCoord(globCoord,lp,srcRegions_);
    }

    if(!sourceElem){
      WARN("Could not find an element for the given local point! Setting it to zero...");
      return;
    }

    shared_ptr<ElemShapeMap> esm = domain->GetGrid(gridId_)->GetElemShapeMap( sourceElem, true );
    LocPointMapped lpmS;
    lpmS.Set(lp,esm);
    this->GetElemSolution( elemSol, sourceElem->elemNum);
    BaseFE * ptFe = esm->GetBaseFE();
    myOperator_->ApplyOp(ptSol,lpmS,ptFe,elemSol);
    CoefMat = ptSol[0];

    //now obtain factor
    CoefMat *= this->GetFactor(lpmS);
  #else
    EXCEPTION("GridCoefFunction just works with interpolation enabled!");
  #endif

  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum){
    StdVector<UInt> connect;
    Grid* cG = domain->GetGrid(this->gridId_);
    cG->GetElemNodes(connect,eNum);
    sol.Resize(connect.GetSize()*dofDim_);
    for(UInt aNode=0;aNode < connect.GetSize();aNode++){
      UInt nodeIdx = nodeIdxMap_[connect[aNode]];
      for(UInt d = 0; d<dofDim_;++d){
        sol[aNode*dofDim_+d] = solVec_[eqnNumbers_[nodeIdx][d]];
      }
    }
  }

  template<class DATA_TYPE>
  std::string CoefFunctionNodalGrid<DATA_TYPE>::ToString() const {
    return "ToSting is not implemented";
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::Recalculate() {
    //no need to recalculate right now
   // UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
    bool interp;
    Double factor;
    UInt curStep = GetStepNum(interp,factor);

    if(interp){
      ReadSolution(curStep,solVecFuture_);

      if(curStep!=1)
        ReadSolution(curStep-1,solVec_);

      for(UInt i = 0;i<numEqns_;++i){
        solVec_[i] += factor*solVecFuture_[i];
      }
      this->curStep_ = curStep;
    }else{
      this->curStep_ = curStep;
      ReadSolution(this->curStep_,solVec_);
    }
  }

  template<class DATA_TYPE>
  UInt CoefFunctionNodalGrid<DATA_TYPE>::GetStepNum(bool & interpolateT,Double & iFactor){
    Double aTimeFreq = this->mp_->Eval(this->mHandleStep_);
    curTStep_ = aTimeFreq;
    UInt step = 0;
    //ok find makes no sense here, we need to iterate over it and
    // apply some tolerance..
    std::map<UInt,Double>::iterator stepIter = this->stepValueMap_.begin();
    for(;stepIter != this->stepValueMap_.end(); ++stepIter){
      if(abs(stepIter->second - aTimeFreq) < 1e-8){
        break;
      }
    }
    if(stepIter==this->stepValueMap_.end()){
      //ok there we have to trigger temporal interpolation
      //for now we just throw an exception
      stepIter = this->stepValueMap_.begin();
      Double oldTime = stepIter->second;
      UInt pos = 1;
      if(curTStep_ > oldTime){
        for(;stepIter != this->stepValueMap_.end(); ++stepIter,++pos){
          if(oldTime < curTStep_ && stepIter->second > curTStep_){
            break;
          }else{
            oldTime = stepIter->second;
          }
        }
      }
      if(pos==1){
        interpolateT = true;
        Double dt = this->stepValueMap_[1];
        iFactor = (curTStep_  - this->stepValueMap_[1])/dt;
      }else{
        interpolateT = true;
        Double dt = this->stepValueMap_[this->curStep_+1] - this->stepValueMap_[this->curStep_];
        iFactor = (curTStep_  - this->stepValueMap_[this->curStep_])/dt;
      }
      if(stepIter==this->stepValueMap_.end()){
        step = this->curStep_;
        iFactor = 1;
      }else{
        step = stepIter->first;
      }
    }else{
      step = stepIter->first;
      interpolateT = false;
      iFactor = 0.0;
    }
    return step;
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::ReadXmlNode(PtrParamNode configNode){

    std::string solString;
    configNode->Get("input")->GetValue("inputId", inputId_);
    configNode->Get("input")->GetValue("gridId", gridId_);
    configNode->Get("input")->GetValue("quantity", solString);
    solType_ = SolutionTypeEnum.Parse(solString );
    srcGrid_ = domain->GetGrid(gridId_);
    std::string factorString = "1.0";
    configNode->GetValue("factor",factorString);
    factorFnc_ = CoefFunction::Generate(Global::REAL,factorString,"");

    this->dependType_ = CoefFunction::GENERAL;

    //configNode->GetValue("factor",factor_);

    PtrParamNode regionListNode = configNode->Get("regionList", ParamNode::PASS );

    ParamNodeList regionNodes;

    //configNode->GetValue("globalTolerance",globalTolerance_);
    //configNode->GetValue("localTolerance",localTolerance_);

    if( regionListNode)
      regionNodes = regionListNode->GetList("region");

    UInt numRegions = regionNodes.GetSize();
    // iterate over all regions
    for(UInt i = 0; i < numRegions; ++i ) {
      // get data from node
      std::string region;
      regionNodes[i]->GetValue( "srcRegionName", region );
      srcRegions_.insert(srcGrid_->GetRegion().Parse(region) );
      this->srcRegionNames_.insert(region);
    }


    UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
    //obtain availResults
    StdVector<shared_ptr<ResultInfo> > results;
    domain->GetResultHandler()->GetResultTypes(this->inputId_,aSeqStep,results,false);
    //now we search for the appropriate result
    for(UInt i = 0;i<results.GetSize();i++){
    if( results[i]->resultType == this->solType_ ) {
      this->resultInfo_ = results[i];
      break;
    }
    }
    if(!this->resultInfo_.get())
      EXCEPTION("Could not find result " << solString);



    domain->GetResultHandler()->GetStepValues(this->inputId_,aSeqStep,this->resultInfo_,this->stepValueMap_,false);
    MapEqns();
    ReadSolution(this->stepValueMap_.begin()->first,solVec_);
    //feFunct_ = domain->GetResultHandler()->GetFeFunction<DATA_TYPE>(inputId_,aSeqStep,stepValueMap_.begin()->first,solType_,srcRegionNames_,factorFnc_);
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::MapEqns(){
    UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
    shared_ptr<BaseResult> Bres = domain->GetResultHandler()->GetResult( this->inputId_, aSeqStep ,this->stepValueMap_.begin()->first, this->solType_, *(this->srcRegionNames_.begin()) );
    Grid * ptGrid = dynamic_cast<Result<DATA_TYPE>&>(*Bres).GetEntityList()->GetGrid();

    shared_ptr<ResultInfo> actInfo = dynamic_cast<Result<DATA_TYPE>&>(*Bres).GetResultInfo();

    std::set<std::string>::iterator regIter = srcRegionNames_.begin();
    UInt pos = 0;
    for( ; regIter != srcRegionNames_.end(); ++regIter) {
      StdVector<UInt> nList;
      RegionIdType curId = ptGrid->GetRegion().Parse(*regIter);
      ptGrid->GetNodesByRegion(nList,curId);
      for(UInt i=0; i<nList.GetSize(); i++){
        nodeIdxMap_[nList[i]] = pos++;
      }
    }
    //determine number of dofs
    if(actInfo->dofNames.GetSize() == 0){
      dofDim_ = 1;
    }else{
      dofDim_ = actInfo->dofNames.GetSize();
    }

    eqnNumbers_.Resize(pos,StdVector<UInt>(dofDim_));
    std::map<UInt,UInt>::iterator idxIter = nodeIdxMap_.begin();
    pos = 0;
    UInt eqnNr = 0;
    for(;idxIter!=nodeIdxMap_.end();++idxIter,++pos){
      for(UInt d = 0; d < dofDim_;d++){
        eqnNumbers_[pos][d] = eqnNr++;
      }
    }
    numEqns_ = eqnNr;
    solVec_.Resize(eqnNr);
    solVec_.Init();
  }

  template<class DATA_TYPE>
  void CoefFunctionNodalGrid<DATA_TYPE>::ReadSolution(UInt step,Vector<DATA_TYPE> & sol){
    UInt aSeqStep = domain->GetDriver()->GetActSequenceStep();
    std::set<std::string>::iterator regIter = this->srcRegionNames_.begin();
    sol.Resize(numEqns_);
    sol.Init();
    for( UInt i = 0; regIter != srcRegionNames_.end(); ++i,++regIter) {
      shared_ptr<BaseResult> Bres = domain->GetResultHandler()->GetResult( this->inputId_, aSeqStep , step , this->solType_, *(regIter) );
      Grid * ptGrid = dynamic_cast<Result<DATA_TYPE>&>(*Bres).GetEntityList()->GetGrid();
      shared_ptr<EntityList> regionList = Bres->GetEntityList();
      shared_ptr<EntityList> nodeList = ptGrid->GetEntityList(EntityList::NODE_LIST, regionList->GetName());
      EntityIterator it= nodeList->GetIterator();

      Vector<DATA_TYPE> & resVec = dynamic_cast<Result<DATA_TYPE>&>(*Bres).GetVector();
      StdVector<UInt> eqns;
      UInt locPos = 0;

      for( it.Begin(); !it.IsEnd(); it++ ) {
        UInt idx = nodeIdxMap_[it.GetNode()];
        eqns = eqnNumbers_[idx];
        for(UInt d = 0; d<eqns.GetSize();++d){
          solVec_[eqns[d]] = resVec[locPos++];
        }
      }
    }
  }

  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class CoefFunctionNodalGrid<Double>;
    template class CoefFunctionNodalGrid<Complex>;
  #endif


}
