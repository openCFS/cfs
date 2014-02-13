// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGridNodal.cc 
 *       \brief    Implementation of nodal grid inteprolation 
 *
 *       \date     11/21/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#include <def_expl_templ_inst.hh>

#include "CoefFunctionGridNodal.hh"

#include "FeBasis/FeSpace.hh"
#include <boost/lexical_cast.hpp>
#include <boost/tr1/type_traits.hpp>
#include <set>
#include <string>

namespace CoupledField{


  template<class DATA_TYPE>
  CoefFunctionGridNodal<DATA_TYPE>::CoefFunctionGridNodal(Domain* ptDomain,
                                                          PtrParamNode configNode)
                                   :CoefFunctionGrid(ptDomain, configNode){
    eqnMapComplete_ = false;

    //Set sequence Step according to XML definition
    this->aSeqStep_ = configNode->Get("sequenceStep")->As<UInt>();

    std::string solString = configNode->Get("quantity")->As<std::string>();
    solType_ = SolutionTypeEnum.Parse(solString );
    isComplex_ =  std::tr1::is_same<DATA_TYPE,Complex>::value;
    this->numEqns_ = 0;
    // obtain handle from internal variable coefficient function
    mp_ = domain_->GetMathParser();
    mHandleStep_ = mp_->GetNewHandle(true);

    std::string var;
    if(isComplex_)
      var = "f";
    else
      var = "t";
    mp_->SetExpr(mHandleStep_, var);

    std::string factorString;
    if(configNode->Has("globalFactor")){
      configNode->GetValue("globalFactor",factorString);

      if(isComplex_)
        factorFnc_ = CoefFunction::Generate(mp_,Global::COMPLEX,factorString);
      else
        factorFnc_ = CoefFunction::Generate(mp_, Global::REAL,factorString);
    }else{

      if(isComplex_)
        factorFnc_ = CoefFunction::Generate(mp_,Global::COMPLEX,"1.0");
      else
        factorFnc_ = CoefFunction::Generate(mp_,Global::REAL,"1.0");
    }


  }

  template<class DATA_TYPE>
  CoefFunctionGridNodal<DATA_TYPE>::~CoefFunctionGridNodal(){
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetTensor is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetVector is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetScalar is not implemented here"); 
  }


  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::CreateOperator(UInt spaceDim, UInt dofDim){
    if(dofDim == 1){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,1,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,1,DATA_TYPE>());
    }else if(dofDim == 2){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,2,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,2,DATA_TYPE>());
    }else if(dofDim == 3){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,3,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,3,DATA_TYPE>());
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::CreateDivOperator(UInt spaceDim, UInt dofDim){

    if(spaceDim != dofDim)
      EXCEPTION("CoefFunctionGridNodal<DATA_TYPE>: Divergence need vectorial data!");

    if(spaceDim == 2)
      this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,2,DATA_TYPE>());
    else if(spaceDim == 3)
      this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,3,DATA_TYPE>());
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::MapEqns(){
    //Be careful we determine the current sequence step according to the 
    //current simulation run. This could fail in a multisequence analysis!!!
    //the user should give an argument where to find the results!
    shared_ptr<BaseResult> Bres = domain_->GetResultHandler()->GetResult( inputId_, aSeqStep_ ,stepValueMap_.begin()->first, solType_, *srcRegions_.begin() );

    std::set<std::string>::iterator regIter = srcRegions_.begin();
    UInt pos = 0;
    for( ; regIter != srcRegions_.end(); ++regIter) {
      StdVector<UInt> nList;
      RegionIdType curId = srcGrid_->GetRegion().Parse(*regIter);
      srcGrid_->GetNodesByRegion(nList,curId);

      for(UInt i=0; i<nList.GetSize(); i++){
        nodeIdxMap_[nList[i]] = pos++;
      }
    }

    //catch the case in which the dimDof_ varable is zero
    assert(dimDof_ != 0);

    eqnNumbers_.Resize(pos,StdVector<UInt>(dimDof_));
    std::map<UInt,UInt>::iterator idxIter = nodeIdxMap_.begin();
    pos = 0;
    UInt eqnNr = 0;
    for(;idxIter!=nodeIdxMap_.end();++idxIter,++pos){
      for(UInt d = 0; d < dimDof_;d++){
        eqnNumbers_[pos][d] = eqnNr++;
      }
    }
    numEqns_ = eqnNr;
    solVec_.Resize(eqnNr);
    solVec_.Init();

    eqnMapComplete_ = true;
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum){
    StdVector<UInt> connect;
    srcGrid_->GetElemNodes(connect,eNum);
    sol.Resize(connect.GetSize()*dimDof_);
    for(UInt aNode=0;aNode < connect.GetSize();aNode++){
      UInt nodeIdx = nodeIdxMap_[connect[aNode]];
      for(UInt d = 0; d<dimDof_;++d){
        sol[aNode*dimDof_+d] = solVec_[eqnNumbers_[nodeIdx][d]];
      }
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::ReadSolution(UInt step,Vector<DATA_TYPE> & sol){
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    sol.Resize(numEqns_);
    sol.Init(0.0);
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      shared_ptr<BaseResult> Bres = domain_->GetResultHandler()->GetResult( this->inputId_, this->aSeqStep_ , step , this->solType_, *regIter );
      shared_ptr<EntityList> regionList = Bres->GetEntityList();
      shared_ptr<EntityList> nodeList = srcGrid_->GetEntityList(EntityList::NODE_LIST, regionList->GetName());
      EntityIterator it= nodeList->GetIterator();
      Vector<DATA_TYPE> resVec;
      try{
        Result<DATA_TYPE>* myResult = dynamic_cast<Result<DATA_TYPE>* >(Bres.get());
        resVec =   myResult->GetVector();
      }catch(...){
        EXCEPTION("Cannot cast to desired vector type. Are you trying to load real data into a harmonic computation?");
      }

      //Vector<DATA_TYPE> curVec = resVec; // = dynamic_cast< Vector<DATA_TYPE>* >(Bres-);
      StdVector<UInt> eqns;
      UInt locPos = 0;

      //TODO> In case of multiple source regions, node results on the interface between the regions appear multiple times
      //The current implementation will overwrite the results obtained from reg1 with the one from reg2
      //this should not be problematic as the results in region 1 and 2 should be identical.
      //but we should be careful!
      if(factorFnc_->GetDependency() == CoefFunction::CONSTANT ||
         factorFnc_->GetDependency() == CoefFunction::TIMEFREQ){
        //determine factor
        DATA_TYPE factor;
        LocPointMapped lpm;
        factorFnc_->GetScalar(factor,lpm);
        //std::cout << "Computed Factor for timestep: " << factor << std::endl;
        for( it.Begin(); !it.IsEnd(); it++ ) {
          UInt idx = nodeIdxMap_[it.GetNode()];
          eqns = eqnNumbers_[idx];
          for(UInt d = 0; d<eqns.GetSize();++d){
            sol[eqns[d]] = resVec[locPos++] * factor;
            //curVec->GetEntry(locPos++,sol[eqns[d]]);
          }
        }
      }else{
        StdVector<Vector<Double> > CoordVec(nodeList->GetSize());
        StdVector< DATA_TYPE > values(nodeList->GetSize());
        UInt node=0;
        for( it.Begin() ; !it.IsEnd() ; it++ , node++ ) {
          UInt curNodeNum = it.GetNode();
          srcGrid_->GetNodeCoordinate(CoordVec[node],curNodeNum,true);
        }


        // CURRENTLY NOT WORKING!
        //factorFnc_->GetScalarValuesAtCoords(CoordVec,values);
        factorFnc_->GetScalarValuesAtCoords(CoordVec,values,this->domain_->GetGrid());
        //values.Init(1.0);
        node = 0;
        for( it.Begin(); !it.IsEnd(); it++,node++ ) {
          UInt idx = nodeIdxMap_[it.GetNode()];
          eqns = eqnNumbers_[idx];
          for(UInt d = 0; d<eqns.GetSize();++d,++locPos){
            sol[eqns[d]] = resVec[locPos] * values[locPos];

            //curVec->GetEntry(locPos++,sol[eqns[d]]);
          }
        }
      }
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::UpdateSolution(){
    if(this->GetDependency() != CoefFunction::CONSTANT){
      UInt stepnumber=0;
      bool needTinterp=false;
      Double factor1,factor2;
      stepnumber = GetStepNum(needTinterp,factor1,factor2);
      if(needTinterp){
        EXCEPTION("StepValue interpolation not supported right now")
      }else{
        //we just read the solution vector
        if(lastStepRead_ != stepnumber){
          this->ReadSolution(stepnumber,this->solVec_);
          lastStepRead_ = stepnumber;
        }
      }
    }
  }

  template<class DATA_TYPE>
  UInt CoefFunctionGridNodal<DATA_TYPE>::GetStepNum(bool & interpolateT,Double & iFactor1, Double & iFactor2){
    //first we cover the Const Case
    if(dependType_ == CoefFunction::CONSTANT){
      interpolateT = false;
      iFactor1 = 0.0;
      iFactor2 = 0.0;
      //we return just the first step from the value map
      //as an future expansion the user could supply
      //the desired time/freqeuncy step...
      return stepValueMap_.begin()->first;
    }

    //Ok, this is the case for which it is possible that we need
    //temporal interpolation
    Double aTimeFreq = this->mp_->Eval(mHandleStep_);
    curTStep_ = aTimeFreq;
    UInt step = 0;
  
    //std::cout << "timestep: " << curTStep_ << std::endl;
    //ok find makes no sense here, we need to iterate over it and
    // apply some tolerance..
    std::map<UInt,Double>::iterator stepIter = stepValueMap_.begin();
    for(;stepIter != stepValueMap_.end(); ++stepIter){
      if(abs(stepIter->second - aTimeFreq) < 1e-4){
        break;
      }
    }
    if(stepIter==stepValueMap_.end()){
      //ok there we have to trigger temporal interpolation
      //for now we just throw an exception
      stepIter = stepValueMap_.begin();
      Double oldTime = stepIter->second;
      UInt pos = 1;
      if(curTStep_ > oldTime){
        for(;stepIter != stepValueMap_.end(); ++stepIter,++pos){
          if(oldTime < curTStep_ && stepIter->second > curTStep_){
            break;
          }else{
            oldTime = stepIter->second;
          }
        }
      }
      interpolateT = true;
//      std::cout << "oldTime = " << oldTime << std::endl;
//      std::cout << "stepIter->second = " << stepIter->second << std::endl;
      Double dt = stepIter->second - oldTime;
      iFactor1 = (stepIter->second  - curTStep_)/dt;
      iFactor2 = (curTStep_  - oldTime)/dt;
  
      if(stepIter==stepValueMap_.end()){
        step = curStep_;
        iFactor1 = 1;
        iFactor2 = 0;
      }else{
        step = stepIter->first;
      }
    }else{
      step = stepIter->first;
      interpolateT = false;
      iFactor1 = 0.0;
      iFactor2 = 0.0;
    }
    return step;
  }
}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionGridNodal<Double>;
  template class CoefFunctionGridNodal<Complex>;
#endif
