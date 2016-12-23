// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VecOpFilter.cc
 *       \brief    <Description>
 *
 *       \date     Oct 19, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "VecOpFilter.hh"

namespace CFSDat{

FilterPtr VecOpFilter::GenerateVectorOperator(PtrParamNode vecOpNode, PtrResultManager resMana){
  FilterPtr newFilter;

 if(vecOpNode->Get("opType")->As<std::string>() == "plus"){
   newFilter = FilterPtr(new CFSDat::GenericVecOpFilter<PlusOp>(0,vecOpNode, resMana));
 } else if(vecOpNode->Get("opType")->As<std::string>() == "minus"){
   newFilter = FilterPtr(new CFSDat::GenericVecOpFilter<MinusOp>(0,vecOpNode, resMana));;
 } else if(vecOpNode->Get("opType")->As<std::string>() == "div"){
   newFilter = FilterPtr(new CFSDat::GenericVecOpFilter<DivOp>(0,vecOpNode, resMana));;
 } else if(vecOpNode->Get("opType")->As<std::string>() == "mult"){
   newFilter = FilterPtr(new CFSDat::GenericVecOpFilter<MultOp>(0,vecOpNode,resMana));;
 }
 return newFilter;
}

template<class Operator>
GenericVecOpFilter<Operator>::GenericVecOpFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                             :VecOpFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  this->res1Name = config->Get("quantity1")->Get("resultName")->As<std::string>();
  this->res2Name = config->Get("quantity2")->Get("resultName")->As<std::string>();
  this->outName = config->Get("output")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
  upResNames.insert(res1Name);
  upResNames.insert(res2Name);
}

template<class Operator>
bool GenericVecOpFilter<Operator>::Run(){
  Double aTF = resultManager_->GetStepValue(outId);
  //set time value and activate the input results
  resultManager_->SetTimeValue(res1Id,aTF);
  resultManager_->SetTimeValue(res2Id,aTF);
  resultManager_->ActivateResult(res1Id);
  resultManager_->ActivateResult(res2Id);
  resultManager_->DeactivateResult(outId);

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }

  //equation numbers for this result
  //here, those equation numbers are equal for all in/out results
  CF::StdVector<UInt> eqnNums;
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(outId,eqnNums);
  returnVec.Init();

  Vector<Double>& res1V = resultManager_->GetResultVector<Double>(res1Id,eqnNums);
  Vector<Double>& res2V = resultManager_->GetResultVector<Double>(res2Id,eqnNums);

  UInt last = (eqnNums.GetSize() == 0)? returnVec.GetSize() : eqnNums.GetSize();

  for(UInt i=0;i<last;++i){
    returnVec[i] = Operator::Apply(res1V[i],res2V[i]);
  }
  resultManager_->ActivateResult(outId);
  return true;
}

template<class Operator>
ResultIdList GenericVecOpFilter<Operator>::SetUpstreamResults(){
  ResultIdList generated;
  //sanity check
  if(filterResIds.GetSize()>1){
    EXCEPTION("Filter can only have one result. This indicates a Bug!")
  }
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  outId = *aIt;
  //we only have one filter result
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
  res1Id = resultManager_->AddResult(res1Name,this->filterTag_);
  res2Id = resultManager_->AddResult(res2Name,this->filterTag_);

  //copy the timeline from input result
  resultManager_->SetTimeLine(res1Id,(*resultManager_->GetExtInfo(outId)->timeLine.get()));
  resultManager_->SetTimeLine(res2Id,(*resultManager_->GetExtInfo(outId)->timeLine.get()));
  generated.Push_back(res1Id);
  generated.Push_back(res2Id);

  return generated;
}

template<class Operator>
void GenericVecOpFilter<Operator>::AdaptFilterResults(){
  //We should check some validity...
  //we can almost copy everything from time input (also scalar, etc.)
  resultManager_->CopyResultData(res1Id,outId);
  resultManager_->SetValid(outId);
}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class GenericVecOpFilter<VecOpFilter::PlusOp>;
  template class GenericVecOpFilter<VecOpFilter::MultOp>;
  template class GenericVecOpFilter<VecOpFilter::DivOp>;
  template class GenericVecOpFilter<VecOpFilter::MinusOp>;
#endif
}
