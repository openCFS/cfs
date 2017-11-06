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

PostLighthillSource::PostLighthillSource(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;

  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  std::string outRes = params_->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outRes);
  upResNames.insert(inRes);
}

PostLighthillSource::~PostLighthillSource(){

}

ResultIdList PostLighthillSource::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

  //add input result to manager
  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  uuids::uuid newId = resultManager_->AddResult(inRes,this->filterTag_);

  //set the timeline of upstream data if already set
  resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(newId);

  return generated;
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

bool PostLighthillSource::Run(){
  // we deactivate every result, except for our own
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();

  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1){
      WARN(" There are still active results when reaching the interpolation filter. This indicates an unexpected use of the pipeline.")
    }
    resultManager_->DeactivateResult(*aIter);
  }
  Double aTF = resultManager_->GetStepValue(filterResIds[0]);
  resultManager_->SetTimeValue(upResIds[0],aTF);
  // now we deactivate our own result and activate the others
  resultManager_->ActivateResult(upResIds[0]);

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }

  CF::StdVector<UInt> eqnNumsRet;
  CF::StdVector<UInt> eqnNumsIn;
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNumsRet);
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNumsIn);

  // check if input can be 2D-data
  if ( inVec.GetSize() % 2 != 0){
    EXCEPTION("PostLigthillSource.cc: Is the input of this filter really 2D? Because it's size isn't even");
  }

  Vector<Double> r;
  r.Resize(inVec.GetSize() / 2);
  for(UInt i = 0; i < r.GetSize(); ++i){
    r[i] = inVec[2 * i];
  }

  returnVec = r;

  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}


}
