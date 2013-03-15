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

namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridNodalDefault<DATA_TYPE>::CoefFunctionGridNodalDefault(PtrParamNode configNode)
                                        :CoefFunctionGridNodal<DATA_TYPE>(configNode){

  //====================================================
  // Determine information about source grid and result
  //====================================================                                          
  this->inputId_ = "default";
  this->gridId_ = "default";
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;


  //obtain the sequence step for result
  //plz note we take here the same value as for the current simulation
  //TODO: Parametrize by XML
  this->aSeqStep_ = domain->GetDriver()->GetActSequenceStep();

  DetermineResult(this->inputId_,this->aSeqStep_);


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
  CreateOperator(dimGrid,this->dimDof_);

  //check what kind of dependency is there...
  std::string dependString;
  configNode->GetValue("dependtype",dependString);
  this->dependType_ = CoefFunction::CoefDependType_.Parse(dependString);
 

  this->curStep_ = 1;
  this->curTStep_ = 0;
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
  if(lpm.isSurface)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;

  Vector<DATA_TYPE> elemSol;
  if(this->dimType_ == CoefFunction::SCALAR)
    CoefMat.Resize(1);
  else if (this->dimType_ == CoefFunction::VECTOR)
    CoefMat.Resize(this->dimDof_);
  
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

  if(lpm.isSurface)
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
    if(this->eqnMapComplete_)
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

    //====================================================
    // Create Data structures for easy solution access
    //====================================================                                          
    //in this special class we start with the equation mapping right
    //after the first call to get entities as this method should not be called twice!
    this->MapEqns();
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
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
