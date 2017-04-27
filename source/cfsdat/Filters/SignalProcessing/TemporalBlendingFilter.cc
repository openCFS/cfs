// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TemporalBlendingFilter.cc
 *       \brief    <Description>
 *
 *       \date     Apr 18, 2017
 *       \author   kroppert
 */
//================================================================================================

#include "TemporalBlendingFilter.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CFSDat{

TemporalBlendFilter::TemporalBlendFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
: BaseFilter(numWorkers,config, resMan){

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  this->resName = config->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  this->outName = config->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
  upResNames.insert(resName);

  expr_ = config->Get("blendingFunction")->Get("temporalFunction")->As<std::string>();
  mp_ = new MathParser();
}

TemporalBlendFilter::~TemporalBlendFilter(){
  delete mp_;
  mp_ = NULL;
}

bool TemporalBlendFilter::Run(){
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

  CF::StdVector<UInt> eqnNums;
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  Vector<Double> inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);

  Double currentStep = resultManager_->GetStepValue(upResIds[0],0);

  MathParser::HandleType hand;
  hand = mp_->GetNewHandle();
  mp_->SetExpr(hand, expr_);
  // get variables in xml-string, if there is anything different than "t" then exception
  StdVector<std::string> vars;
  mp_->GetExprVars(hand, vars);
  if(vars.GetSize() != 1 || vars[0] != "t"){
    EXCEPTION("Wrong variable in temporalBlending...only time dependent functions f(t) are allowed!!");
  }
  mp_->SetValue(hand, "t", currentStep);
  Double ev = mp_->Eval(hand);
  mp_->ReleaseHandle(hand);

  inVec.ScalarMult(ev);
  returnVec = inVec;

  resultManager_->ActivateResult(filterResIds[0]);
  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }
  return true;
}


ResultIdList TemporalBlendFilter::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  outId = *aIt;
  //we only have one filter result, add input result to manager
  resId = resultManager_->AddResult(resName,this->filterTag_);
  //set the timeline of upstream data if already set
  resultManager_->SetTimeLine(resId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(resId);
  return generated;
}

void TemporalBlendFilter::AdaptFilterResults(){
  //we would need to copy related data to our filter results
  resultManager_->CopyResultData(resId,outId);
  resultManager_->SetValid(outId);
}

}

