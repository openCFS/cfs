// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefault.cc 
 *       \brief    Implementation File
 *
 *       \date     Jan. 20, 2013
 *       \author   Andreas Hueppe
 */
//================================================================================================


#include "CoefFunctionGridNodalDefault.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridNodalDefault<DATA_TYPE>::CoefFunctionGridNodalDefault(PtrParamNode configNode,PtrParamNode curInfo)
                                        :CoefFunctionGridNodal<DATA_TYPE>(configNode){

  //====================================================
  // Determine information about source grid and result
  //====================================================                                          
  this->inputId_ = "default";
  this->gridId_ = "default";
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;
  this->conservativeReady_ = false;
  this->srcIsSurface_ = false;



  this->extDataInfo_ = curInfo->Get("defaultGrid");
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue("noInterpolation");

  std::string factorString = this->factorFnc_->ToString();
  this->extDataInfo_->Get("interpolation")->Get("factor",ParamNode::INSERT)->SetValue(factorString);

  this->srcGrid_ = domain->GetGrid();

  //lets determine the destination region and set it to our source regions
  std::string destreg = configNode->GetParent()->GetParent()->Get("name")->As<std::string>();
  //obtain entitylist from grid and add it
  shared_ptr<EntityList> curList = this->srcGrid_->GetEntityList(EntityList::ELEM_LIST,destreg);
  this->DetermineResult(this->inputId_,this->aSeqStep_);
  this->dimDof_ = this->resultInfo_->dofNames.GetSize();
  // Determine which steps are available
  domain->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);
  this->curStep_ = this->stepValueMap_.begin()->first;
  this->curTStep_ = this->stepValueMap_.begin()->second;
  //this->AddEntityList(curList);





  //====================================================
  // Create interpolation in time and space
  //====================================================                                          
  UInt dimGrid = this->srcGrid_->GetDim();
  if(this->resultInfo_->entryType == ResultInfo::SCALAR){
    this->dimDof_ = 1;
    this->dimType_ = CoefFunction::SCALAR;
  }else if(this->resultInfo_->entryType == ResultInfo::VECTOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::VECTOR;
  }else if(this->resultInfo_->entryType == ResultInfo::TENSOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::TENSOR;
  }
  this->CreateOperator(dimGrid,this->dimDof_);

  //check what kind of dependency is there...

  //GridcoefFunctions are always general!
  //std::string dependString;
  //configNode->GetValue("dependtype",dependString);
  //this->dependType_ = CoefFunction::CoefDependType_.Parse(dependString);
  this->dependType_ = CoefFunction::GENERAL;

}

// ========================
//  ACCESS METHODS
// ========================
template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  EXCEPTION("TENSORIAL RESULTS NOT SUPPORTED FOR INTERPOLATION FROM EXTERNAL GRIDS");
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  //ok we always take the volume element
  //because even if we are coping here with surfaces, e.g. for boundary conditions,
  //the volume element connectivity should give also the boundary nodes and
  //we do not cover special, nonconforming cases here...
  if(lpm.isSurface && !this->srcIsSurface_)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;

  Vector<DATA_TYPE> elemSol;
  if(this->dimType_ == CoefFunction::SCALAR)
    CoefMat.Resize(1);
  else if (this->dimType_ == CoefFunction::VECTOR)
    CoefMat.Resize(this->dimDof_);
  
  if(!sourceElem){
    EXCEPTION("Could not determine source element.")
  }

  this->GetElemSolution( elemSol, sourceElem->elemNum);

  shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  BaseFE * ptFe = esm->GetBaseFE();
  if(lpm.isSurface)
    this->myOperator_->ApplyOp(CoefMat,(*lpm.lpmVol),ptFe,elemSol);
  else
    this->myOperator_->ApplyOp(CoefMat,lpm,ptFe,elemSol);



}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);
  //this is really simple. we just take the nodal result and
  //interpolate it to lpm

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  if(lpm.isSurface && !this->srcIsSurface_)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;
  Vector<DATA_TYPE> elemSol;
  Vector<DATA_TYPE> ptSol(1);
  ptSol.Init();

  sourceElem = lpm.ptEl;
  this->GetElemSolution( elemSol, sourceElem->elemNum);
  shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  BaseFE * ptFe = esm->GetBaseFE();
  if(lpm.isSurface)
    this->myOperator_->ApplyOp(ptSol,(*lpm.lpmVol),ptFe,elemSol);
  else
    this->myOperator_->ApplyOp(ptSol,lpm,ptFe,elemSol);

  CoefMat = ptSol[0];
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::AddEntityList(shared_ptr<EntityList> ent){
    //this function may not be called if eqnMapping has already
    //been performed because there will be exactly one srcRegion to one target 
    //region and both are identical
    if(this->eqnMapComplete_ && !this->entities_.Contains(ent))
      EXCEPTION("Call to AddEntities after eqnMapping has been performed. This should not happen!");
    
    //first we check if its a region list.
    //otherwise this will be shit
    if(ent->GetDefineType() != EntityList::REGION)
      EXCEPTION("CoefFunctionGridNodalDefault can only be used with regions");

    if(!this->entities_.Contains(ent)){
      this->entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }

    //in this special coeffunction,
    //the sourceRegion is equal to the destination
    this->srcRegions_.insert(ent->GetName());

    if(ent->GetType() == EntityList::SURF_ELEM_LIST){
      this->srcIsSurface_ = true;
    }


    //====================================================
    // Create Data structures for easy solution access
    //====================================================                                          
    //in this special class we start with the equation mapping right
    //after the first call to get entities as this method should not be called twice!
    this->MapEqns();
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    this->extDataInfo_->Get("RegionList")->Get("SourceRegion")->Get("name")->SetValue(ent->GetName());
  }

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetScalarValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<DATA_TYPE >  & vals){
  assert(this->dimType_ == CoefFunction::SCALAR);

  vals.Resize(points.GetSize());
  vals.Init();
  StdVector< const Elem* > elements;
  StdVector<LocPoint> locals;
  this->GetElemsForPoints(points,elements,locals);
  shared_ptr<ElemShapeMap> esm;
  LocPointMapped lpm;
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    const Elem* sourceElem = elements[curPoint];
    esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
    LocPoint lp = locals[curPoint];
    lpm.Set(lp,esm,1.0);
    this->GetScalar(vals[curPoint],lpm);
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetVectorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                      StdVector<Vector<DATA_TYPE> >  & vals){
  vals.Resize(points.GetSize(),Vector<DATA_TYPE>(this->dimDof_));

  vals.Init();
  StdVector< const Elem* >  elements;
  StdVector<LocPoint>  locals;
  this->GetElemsForPoints(points,elements,locals);
  shared_ptr<ElemShapeMap> esm;
  LocPointMapped lpm;
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    const Elem* sourceElem = elements[curPoint];
    esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
    LocPoint lp = locals[curPoint];
    lpm.Set(lp,esm,1.0);
    this->GetVector(vals[curPoint],lpm);
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetTensorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<Matrix<DATA_TYPE> >  & vals){

  assert(this->dimType_ == CoefFunction::TENSOR);
  EXCEPTION("GetTensorValuesAtPoints is not implemented");
  //vals.Resize(points.GetSize(),Matrix<DATA_TYPE>(this->numRows_,this->numCols_));
  //vals.Init();
  //StdVector< const Elem* > elements;
  //StdVector<LocPoint> locals;
  //this->GetElemsForPoints(points,elements,locals);
  //shared_ptr<ElemShapeMap> esm;
  //LocPointMapped lpm;
  //for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
  //  const Elem* sourceElem = elements[curPoint];
  //  esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  //  LocPoint lp = locals[curPoint];
  //  lpm.Set(lp,esm,1.0);
  //  this->GetTensor(vals[curPoint],lpm);
  //}
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::MapConservative( shared_ptr<FeSpace> targetSpace,
                                                                 Vector<DATA_TYPE>& feFncVec){

  //if the targetSpace is also using grid ordering, we just take the nodal values and map them to the
  //target equations
  //if the target space has arbitrary order, this does not really make sense, so we issue an exception
  //the sources then need more some sort of  distribution than a clean mapping
  //this is the old problem with conservative interpolation. what happens if the target grid is finer than the source grid?
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    FeSpace::MappingType curtype = targetSpace->GetMapType(curId);
    if(curtype!=FeSpace::GRID){
      EXCEPTION("You are trying to map a nodal grid to a polynomial space. This is not supported!")
    }
  }

  //First, update the solution vector
  this->UpdateSolution();

  //IDEA> create vector of pairs and associate each eqnNumber of feFncVec a index of own solvec
  //a little more memory but very efficient updates
  if(!this->conservativeReady_){
    BuildNodeIdxAssoc(targetSpace);
  }
  for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i){
    const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
    feFncVec[curP.first] += this->solVec_[curP.second];
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::BuildNodeIdxAssoc(shared_ptr<FeSpace> targetSpace){
  //loop over the entitylist and obtain for each element the equation numbers
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    shared_ptr<NodeList> actSDList( new NodeList(this->srcGrid_ ) );
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    actSDList->SetNodesOfRegion( curId );
    EntityIterator ents = actSDList->GetIterator();
    this->fctSolAssoc_.Reserve(this->fctSolAssoc_.GetSize()+actSDList->GetSize());
    while(!ents.IsEnd()){
      //obtain eqn and node number
      std::pair<UInt,UInt> curP;
      StdVector<Integer> eqns;
      targetSpace->GetEqns(eqns,ents);
      UInt curNode = ents.GetNode();
      UInt solIdx = this->nodeIdxMap_[curNode];
      //safety check
      if(this->eqnNumbers_[solIdx].GetSize() != eqns.GetSize()){
        WARN("Detected incompatibility during external data mapping...")
        ents++;
        continue;
      }

      for(UInt i=0; i<eqns.GetSize(); i++){
        if(eqns[i] > 0){
          curP.first = eqns[i]-1;
          curP.second = this->eqnNumbers_[solIdx][i];
          this->fctSolAssoc_.Push_back(curP);
        }
      }
      ents++;
    }
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetElemsForPoints(const StdVector<Vector<Double> >  & points,
                                                                  StdVector< const Elem* > & elements,
                                                                  StdVector<LocPoint> & locals){
  //build up set of source regions
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  std::set<RegionIdType> scrRegIds;
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    scrRegIds.insert(curId);
  }

  //this is unfortunate but we need to search even for the defaultGridCase...
  this->srcGrid_->GetElemsAtGlobalCoords( points,
                                          locals,
                                          elements,
                                          scrRegIds);

}

}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionGridNodalDefault<Double>;
  template class CoefFunctionGridNodalDefault<Complex>;
#endif
