// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridElemDefaultPrecice.cc 
 *       \brief    Implementation File
 *
 *       \date     Feb. 18, 2025
 *       \author   Klaus Roppert
 */
//================================================================================================

#include "CoefFunctionGridElemDefaultPrecice.hh"
#include "Utils/preciceAdapter/IPreciceAdapter.hh"
#include "FeBasis/FeSpace.hh"

#ifdef USE_PRECICE



namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridElemDefaultPrecice<DATA_TYPE>::CoefFunctionGridElemDefaultPrecice(Domain* ptDomain,
                                                                      PtrParamNode configNode,PtrParamNode curInfo,
                                                                      shared_ptr<RegionList> regions,
                                                                      ResultInfo::EntryType type)
                                                    :CoefFunctionGridElem<DATA_TYPE>(ptDomain, configNode, regions){
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
  this->dimDof_ = this->resultInfo_->dofNames.GetSize();
  // Determine which steps are available
  //this->domain_->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);

  this->SetRegions(regions);
  this->InitSolVec();

  preciceAdapter_ = ptDomain->GetPreciceAdapter();

  if(type == ResultInfo::SCALAR){
    this->dimDof_ = 1;
    this->dimType_ = CoefFunction::SCALAR;
  }else if(type == ResultInfo::VECTOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::VECTOR;
  }
}

// ========================
//  ACCESS METHODS
// ========================

template<typename DATA_TYPE>
void CoefFunctionGridElemDefaultPrecice<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);

  const Elem* sourceElem = NULL;

  //this->UpdateSolution();
  
  sourceElem = lpm.ptEl;

  Vector<DATA_TYPE> elemSol;
  
  if(this->dimType_ == CoefFunction::SCALAR)
    CoefMat.Resize(1);
  else if (this->dimType_ == CoefFunction::VECTOR)
    CoefMat.Resize(this->dimDof_);
  
  if(!sourceElem){
    EXCEPTION("Could not determine source element.")
  }

  // handle this with the preciceAdapter_
  elemSol = preciceAdapter_->GetElemResult(this->solType_, sourceElem->elemNum);
  UInt nComp = std::min(static_cast<UInt>(elemSol.GetSize()), this->dimDof_);
  for (UInt i = 0; i < nComp; ++i) {
    CoefMat[i] = elemSol[i];
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridElemDefaultPrecice<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);
  //this is really simple. we just take the nodal result and
  //interpolate it to lpm

  //this->UpdateSolution();

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;
  
  sourceElem = lpm.ptEl;
  Vector<DATA_TYPE> elemSol;
  Vector<DATA_TYPE> ptSol(1);
  ptSol.Init();

  // handle this with the preciceAdapter_
  elemSol = preciceAdapter_->GetElemResult(this->solType_, sourceElem->elemNum);

  CoefMat = elemSol[0];
}

template<typename DATA_TYPE>
void CoefFunctionGridElemDefaultPrecice<DATA_TYPE>::SetRegions(shared_ptr<RegionList> regions){
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
void CoefFunctionGridElemDefaultPrecice<DATA_TYPE>::DetermineResult(std::string inputID,UInt seqStep){
  //obtain availResults and search for the requested one
  StdVector<shared_ptr<ResultInfo> > results;
  std::map<shared_ptr<BaseResult>, shared_ptr<ResultHandler::ResultContext> > test = *this->domain_->GetResultHandler()->GetResultContexts();
  for (const auto &entry : test) {
    std::cout << "Result name: " << entry.first->GetResultInfo()->resultName << std::endl;
      
    if( entry.first->GetResultInfo()->resultType == this->solType_ ) {
      this->resultInfo_ = entry.first->GetResultInfo();
      break;
    }
  }

  std::string solString = SolutionTypeEnum.ToString(this->solType_);
  if((!this->resultInfo_.get()) && !this->domain_->GetPreciceAdapter())
	  EXCEPTION("Can not find result " << solString << " in given file!");
}


template class CoefFunctionGridElemDefaultPrecice<Double>;
template class CoefFunctionGridElemDefaultPrecice<Complex>;

} // namespace CoupledField



#else

namespace CoupledField {

template<typename T>
CoefFunctionGridElemDefaultPrecice<T>::CoefFunctionGridElemDefaultPrecice(
    Domain* domain,
    boost::shared_ptr<ParamNode> configNode,
    boost::shared_ptr<ParamNode> curInfo,
    boost::shared_ptr<RegionList> region,
    ResultInfo::EntryType /*entryType*/)
  : CoefFunctionGridElem<T>(domain, configNode, region)
{
    // Dummy constructor: no additional functionality.
}

template<typename T>
void CoefFunctionGridElemDefaultPrecice<T>::GetVector(Vector<T>& vec, const LocPointMapped& /*lpm*/)
{
    // Dummy implementation: simply resize and fill with zeros
    vec.Resize(1); // or an appropriate dimension
    vec[0] = T(); // default-constructed value
}

template<typename T>
void CoefFunctionGridElemDefaultPrecice<T>::GetScalar(T& scalar, const LocPointMapped& /*lpm*/)
{
    scalar = T(); // default value
}

template<typename T>
void CoefFunctionGridElemDefaultPrecice<T>::DetermineResult(std::string /*name*/, UInt /*id*/)
{
    return ;
}

template<typename T>
void CoefFunctionGridElemDefaultPrecice<T>::SetRegions(boost::shared_ptr<RegionList> /*region*/)
{
    // Do nothing.
}

// Explicit instantiation for the types used in the project.
template class CoefFunctionGridElemDefaultPrecice<Double>;
template class CoefFunctionGridElemDefaultPrecice<Complex>;

} // namespace CoupledField

#endif