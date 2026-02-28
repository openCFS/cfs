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
CoefFunctionGridNodalDefault<DATA_TYPE>::CoefFunctionGridNodalDefault(Domain* ptDomain,
                                                                      PtrParamNode configNode,PtrParamNode curInfo,
                                                                      shared_ptr<RegionList> regions,
                                                                      ResultInfo::EntryType type)
                                        :CoefFunctionGridNodal<DATA_TYPE>(ptDomain, configNode, regions){
  //====================================================
  // Determine information about source grid and result
  //====================================================                                          
  this->inputId_ = "default";
  this->gridId_ = "default";
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;
  this->conservativeReady_ = false;
  this->srcIsSurface_ = false;
  
  this->extDataInfo_ = curInfo->Get("defaultGrid",ParamNode::APPEND);
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue("noInterpolation");

  this->srcGrid_ = this->domain_->GetGrid();

  //lets determine the destination region and set it to our source regions
  this->DetermineResult(this->inputId_,this->aSeqStep_);
  this->dimDof_ = this->resultInfo_->dofNames.GetSize();
  // Determine which steps are available
  this->domain_->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);

  //====================================================
  // Create interpolation in time and space
  //====================================================                                          
  UInt dimGrid = this->srcGrid_->GetDim();
  //if(this->resultInfo_->entryType == ResultInfo::SCALAR){
  if(type == ResultInfo::SCALAR){
    this->dimDof_ = 1;
    this->dimType_ = CoefFunction::SCALAR;
  }else if(type == ResultInfo::VECTOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::VECTOR;
  }else if(type == ResultInfo::TENSOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::TENSOR;
    if ( ptDomain->GetGrid()->IsAxi() ) {
      EXCEPTION("Axisymmtric case not supported yet!");
    }
    //4 entries: (2x2) tensor;  9 entries (3x3) tensor
    if ( this->dimDof_ == 4) {
      this->numRowTensor_ = 2;
      this->numColTensor_ = 2;
    }
    else if( this->dimDof_ == 9) {
      this->numRowTensor_ = 3;
      this->numColTensor_ = 3;
    }
    else
      EXCEPTION("CoefFunctionGridNodalDefault: Just tensor with 9 or 6 entries known");
  }
  this->CreateOperator(dimGrid,this->dimDof_);

  this->SetRegions(regions);
  this->InitSolVec();
  this->WriteGlobalFactorsToXML(configNode);
}

// ========================
//  ACCESS METHODS
// ========================
template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  this->UpdateSolution();
  //ok we always take the volume element
  //because even if we are coping here with surfaces, e.g. for boundary conditions,
  //the volume element connectivity should give also the boundary nodes and
  //we do not cover special, nonconforming cases here...
  if(lpm.isSurface && !this->srcIsSurface_)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;

  Vector<DATA_TYPE> elemSol;
  Vector<DATA_TYPE> CoefVec;
  CoefVec.Resize(this->dimDof_);
  
  if(!sourceElem){
    EXCEPTION("Could not determine source element.")
  }

  this->GetElemSolution( elemSol, sourceElem->elemNum);

  shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  BaseFE *ptFe = esm->GetBaseFE();

  if(lpm.isSurface)
    this->myOperator_->ApplyOp(CoefVec,(*lpm.lpmVol),ptFe,elemSol);
  else
    this->myOperator_->ApplyOp(CoefVec,lpm,ptFe,elemSol);
 
  //convert interpolated vector to tensor
  UInt k = 0;
  CoefMat.Resize(this->numRowTensor_,this->numColTensor_);
  for (UInt i=0; i<this->numRowTensor_; i++){
    for (UInt j=0; j<this->numColTensor_; j++){
      CoefMat[i][j] = CoefVec[k];
      k++;
    }
  }

  //EXCEPTION("TENSORIAL RESULTS NOT SUPPORTED FOR INTERPOLATION FROM EXTERNAL GRIDS");
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  this->UpdateSolution();
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

  this->UpdateSolution();

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
void CoefFunctionGridNodalDefault<DATA_TYPE>::SetRegions(shared_ptr<RegionList> regions){
  Grid* grid = regions->GetGrid();
  StdVector<RegionIdType> regIDs = regions->GetRegionIds();
  std::stringstream ss;
  for (UInt i = 0; i < regIDs.GetSize(); i++) {
    std::string name = grid->GetRegionName(regIDs[i]);
    this->srcRegions_.insert(name);
    ss << name;
    if (i < regIDs.GetSize() - 1) {
      ss << ",";
    }
  }
  this->allSrcRegionNames_ = ss.str();
  
  this->srcIsSurface_ = false;
  for (UInt i = 0; i < this->entities_.GetSize(); i++) {
    if(this->entities_[i]->GetType() == EntityList::SURF_ELEM_LIST){
      this->srcIsSurface_ = true;
    }
    std::string name = this->entities_[i]->GetName();
    this->extDataInfo_->Get("RegionList")->Get("SourceRegion",ParamNode::APPEND)->Get("name")->SetValue(name);
  }
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
  
  // map values
  this->BuildNodeIdxAssoc(targetSpace);
#pragma omp parallel for if (this->numEqns_ > OMP_THRESHOLD/2)
  for(Integer i=0;i< (Integer) this->numEqns_;++i){
    if (this->copyValueIndex_[i]) {
      feFncVec[this->valueTargetIndex_[i]] = this->solVec_[i];
    }
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::BuildNodeIdxAssoc(shared_ptr<FeSpace> targetSpace){
  if(this->conservativeReady_){
    return;
  }
  this->copyValueIndex_.resize(this->numEqns_, false);
  this->valueTargetIndex_.Resize(this->numEqns_, 0);
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    shared_ptr<NodeList> actSDList( new NodeList(this->srcGrid_ ) );
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    actSDList->SetNodesOfRegion( curId );
    EntityIterator ents = actSDList->GetIterator();
    while(!ents.IsEnd()){
      StdVector<Integer> eqns;
      targetSpace->GetEqns(eqns,ents);
      UInt curNode = ents.GetNode();
      //safety check
      if(this->dimDof_ != eqns.GetSize()){
        WARN("Detected incompatibility during external data mapping...")
        ents++;
        continue;
      }

      for(UInt i=0; i<eqns.GetSize(); i++){
        if(eqns[i] > 0){
          this->copyValueIndex_[curNode*this->dimDof_+i] = true;
          this->valueTargetIndex_[curNode*this->dimDof_+i] = eqns[i]-1;
        }
      }
      ents++;
    }
  }
  this->conservativeReady_ = true;
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefault<DATA_TYPE>::GetElemsForPoints(const StdVector<Vector<Double> >  & points,
                                                                  StdVector< const Elem* > & elements,
                                                                  StdVector<LocPoint> & locals){
  //build up set of source regions
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  StdVector<shared_ptr<EntityList> > lists; 
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    shared_ptr<ElemList> newList(new ElemList(this->srcGrid_));
    newList->SetRegion(curId);
    lists.Push_back(newList);
  }

  //this is unfortunate but we need to search even for the defaultGridCase...
  this->srcGrid_->GetElemsAtGlobalCoords( points,
                                          locals,
                                          elements,
                                          lists);
}

template class CoefFunctionGridNodalDefault<Double>;
template class CoefFunctionGridNodalDefault<Complex>;

} // namespace CoupledField
