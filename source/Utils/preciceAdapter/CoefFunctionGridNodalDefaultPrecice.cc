// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefaultPrecice.cc
 *       \brief    Implementation File
 *
 *       \date     Feb. 18, 2025
 *       \author   Klaus Roppert
 */
//================================================================================================

#include "CoefFunctionGridNodalDefaultPrecice.hh"
#include "Utils/preciceAdapter/IPreciceAdapter.hh"
#include "FeBasis/FeSpace.hh"

#ifdef USE_PRECICE



namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridNodalDefaultPrecice<DATA_TYPE>::CoefFunctionGridNodalDefaultPrecice(Domain* ptDomain,
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

  this->extDataInfo_ = curInfo->Get("defaultGrid",ParamNode::APPEND);
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue("noInterpolation");

  this->srcGrid_ = this->domain_->GetGrid();

  //lets determine the destination region and set it to our source regions
  this->DetermineResult(this->inputId_,this->aSeqStep_);

  if(type == ResultInfo::SCALAR){
    this->dimDof_ = 1;
    this->dimType_ = CoefFunction::SCALAR;
  }else if(type == ResultInfo::VECTOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::VECTOR;
  }

  //this->dimDof_ = this->resultInfo_->dofNames.GetSize();

  // Determine which steps are available
  //this->domain_->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);

  this->SetRegions(regions);
  this->InitSolVec();

  preciceAdapter_ = ptDomain->GetPreciceAdapter();

}

// ========================
//  ACCESS METHODS
// ========================

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefaultPrecice<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  assert(this->dimType_ != CoefFunction::TENSOR);

  const Elem* sourceElem = lpm.ptEl;
  if (!sourceElem)
    EXCEPTION("CoefFunctionGridNodalDefaultPrecice::GetVector: could not determine source element.");

  CoefMat.Resize(this->dimDof_);
  CoefMat.Init();

  const StdVector<UInt>& nodes = sourceElem->connect;
  const UInt nNodes = nodes.GetSize();
  for (UInt i = 0; i < nNodes; ++i) {
    Vector<Double> nodeSol = preciceAdapter_->GetNodeResult(this->solType_, nodes[i]);
    for (UInt k = 0; k < this->dimDof_; ++k)
      CoefMat[k] += DATA_TYPE(nodeSol[k]);
  }

  DATA_TYPE factor = this->hasConstantFactor_ ? this->constantFactor_ : DATA_TYPE(1);
  const DATA_TYPE inv = factor / DATA_TYPE(nNodes);
  for (UInt k = 0; k < this->dimDof_; ++k)
    CoefMat[k] *= inv;
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefaultPrecice<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);
  //this is really simple. we just take the nodal result and
  //interpolate it to lpm

  //this->UpdateSolution();

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  sourceElem = lpm.ptEl;
  StdVector<UInt> nodes = sourceElem->connect;
  Vector<DATA_TYPE> nodeSol;
  Vector<DATA_TYPE> solAv(1);
  solAv.Init();
  for( UInt i = 0; i < nodes.GetSize(); i++  ) {
    nodeSol = preciceAdapter_->GetNodeResult(this->solType_, nodes[i]);
    solAv[0] += nodeSol[0];
  }
  solAv[0] /= DATA_TYPE(nodes.GetSize());
  DATA_TYPE factor = this->hasConstantFactor_ ? this->constantFactor_ : DATA_TYPE(1);
  CoefMat = factor * solAv[0];
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalDefaultPrecice<DATA_TYPE>::SetRegions(shared_ptr<RegionList> regions){
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

}



template<typename DATA_TYPE>
void CoefFunctionGridNodalDefaultPrecice<DATA_TYPE>::DetermineResult(std::string inputID, UInt seqStep){
  //obtain the registered result contexts and search for the requested one
  std::map<shared_ptr<BaseResult>, shared_ptr<ResultHandler::ResultContext> > resultContexts = *this->domain_->GetResultHandler()->GetResultContexts();
  for (const auto &entry : resultContexts) {
    if( entry.first->GetResultInfo()->resultType == this->solType_ ) {
      this->resultInfo_ = entry.first->GetResultInfo();
      break;
    }
  }

  std::string solString = SolutionTypeEnum.ToString(this->solType_);
  if((!this->resultInfo_.get()) && !this->domain_->GetPreciceAdapter())
    EXCEPTION("Can not find result " << solString << " in given file!");
}


template class CoefFunctionGridNodalDefaultPrecice<Double>;
template class CoefFunctionGridNodalDefaultPrecice<Complex>;

} // namespace CoupledField

#else

namespace CoupledField {

template<typename T>
CoefFunctionGridNodalDefaultPrecice<T>::CoefFunctionGridNodalDefaultPrecice(
    Domain* domain,
    shared_ptr<ParamNode> configNode,
    shared_ptr<ParamNode> curInfo,
    shared_ptr<RegionList> region,
    ResultInfo::EntryType /*entryType*/)
  : CoefFunctionGridNodal<T>(domain, configNode, region)
{
    // Dummy constructor: no additional functionality.
}

template<typename T>
void CoefFunctionGridNodalDefaultPrecice<T>::GetVector(Vector<T>& vec, const LocPointMapped& /*lpm*/)
{
    // Dummy implementation: simply resize and fill with zeros
    vec.Resize(1); // or an appropriate dimension
    vec[0] = T(); // default-constructed value
}

template<typename T>
void CoefFunctionGridNodalDefaultPrecice<T>::GetScalar(T& scalar, const LocPointMapped& /*lpm*/)
{
    scalar = T(); // default value
}

template<typename T>
void CoefFunctionGridNodalDefaultPrecice<T>::DetermineResult(std::string /*name*/, UInt /*id*/)
{
    return ;
}

template<typename T>
void CoefFunctionGridNodalDefaultPrecice<T>::SetRegions(shared_ptr<RegionList> /*region*/)
{
    // Do nothing.
}

// Explicit instantiation for the types used in the project.
template class CoefFunctionGridNodalDefaultPrecice<Double>;
template class CoefFunctionGridNodalDefaultPrecice<Complex>;

} // namespace CoupledField

#endif