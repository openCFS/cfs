// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PostLighthillSource.cc
 *       \brief    <Description>
 *
 *       \date     Feb 16, 2017
 *       \author   kroppert
 */
//================================================================================================

#include "PostLighthillSource.hh"
#include <vector>
#include <algorithm>

namespace CFSDat{

PostLighthillSource::PostLighthillSource(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                     :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;

  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  std::string outRes = params_->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  component_ = params_->Get("extractComponent")->As<int>();
  filtResNames.insert(outRes);
  upResNames.insert(inRes);
}

PostLighthillSource::~PostLighthillSource(){

}

ResultIdList PostLighthillSource::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}


void PostLighthillSource::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Interpolator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things

  if(resultManager_->GetDefOn(upResIds[0]) != ExtendedResultInfo::ELEMENT){
    EXCEPTION("PostLighthillSource.cc: input is not defined on elements, this filter must be placed right after LighthillSourceTerm!");
  }

  resultManager_->SetDefOn(filterResIds[0], resultManager_->GetDefOn(upResIds[0]));

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::SCALAR);
  CF::StdVector<std::string> dofnames;
  dofnames.Push_back("x");
  resultManager_->SetDofNames(filterResIds[0],dofnames);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}

bool PostLighthillSource::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

  // check if input can be 2D-data
  if ( inVec.GetSize() % 2 != 0){
    EXCEPTION("PostLigthillSource.cc: Is the input of this filter really 2D? Because it's size isn't even");
  }

  Vector<Double> r;
  r.Resize(inVec.GetSize() / 2);
  for(UInt i = 0; i < r.GetSize(); ++i){
    r[i] = inVec[2 * i + (component_)];
  }

  returnVec = r;

  return true;
}


}
