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

TemporalBlendFilter::TemporalBlendFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
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

bool TemporalBlendFilter::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Double aTF = resultManager_->GetStepValue(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);


  unsigned int mp_hand = mp_->GetNewHandle();
  mp_->SetExpr(mp_hand, expr_);
  // get variables in xml-string, if there is anything different than "t" then exception
  StdVector<std::string> vars;
  mp_->GetExprVars(mp_hand, vars);
  if(vars.GetSize() != 1 || vars[0] != "t"){
    EXCEPTION("Wrong variable in temporalBlending...only time dependent functions f(t) are allowed!!");
  }
  mp_->SetValue(mp_hand, "t", aTF);
  const Double ev = mp_->Eval(mp_hand);
  mp_->ReleaseHandle(mp_hand);
  const UInt size = inVec.GetSize();
  
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt i = 0; i < size; i++) {
    returnVec[i] = inVec[i] * ev;
  }
  
  return true;
}


ResultIdList TemporalBlendFilter::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  resId = upResNameIds[*upResNames.begin()];
  outId = filterResIds[0];
  return generated;
}

void TemporalBlendFilter::AdaptFilterResults(){
  //we would need to copy related data to our filter results
  resultManager_->CopyResultData(resId,outId);
  resultManager_->SetValid(outId);
}

}

