// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BinOpFilter.cc
 *       \brief    <Description>
 *
 *       \date     Apr 13, 2016
 *       \author   ahueppe
 */
//================================================================================================


#include "BinOpFilter.hh"

#include <algorithm> 

namespace CFSDat{

FilterPtr BinOpFilter::GenerateOperator(PtrParamNode binOpNode, PtrResultManager resMana){
  FilterPtr newFilter;

 if(binOpNode->Get("opType")->As<std::string>() == "plus"){
   newFilter = FilterPtr(new CFSDat::GenericBinOpFilter<PlusOp>(0,binOpNode,resMana));
 } else if(binOpNode->Get("opType")->As<std::string>() == "minus"){
   newFilter = FilterPtr(new CFSDat::GenericBinOpFilter<MinusOp>(0,binOpNode,resMana));;
 } else if(binOpNode->Get("opType")->As<std::string>() == "div"){
   newFilter = FilterPtr(new CFSDat::GenericBinOpFilter<DivOp>(0,binOpNode,resMana));;
 } else if(binOpNode->Get("opType")->As<std::string>() == "mult"){
   newFilter = FilterPtr(new CFSDat::GenericBinOpFilter<MultOp>(0,binOpNode,resMana));;
 }
 return newFilter;
}

template<class Operator>
GenericBinOpFilter<Operator>::GenericBinOpFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                             :BinOpFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  this->res1Name = config->Get("quantity1")->Get("resultName")->As<std::string>();
  this->res2Name = config->Get("quantity2")->Get("resultName")->As<std::string>();
  this->outName = config->Get("output")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
  upResNames.insert(res1Name);
  upResNames.insert(res2Name);
}

template<class Operator>
bool GenericBinOpFilter<Operator>::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  UInt last = returnVec.GetSize();
  Double aTF = resultManager_->GetStepValue(outId);

  Vector<Double>& res1V = GetUpstreamResultVector<Double>(res1Id, aTF);
  Vector<Double>& res2V = GetUpstreamResultVector<Double>(res2Id, aTF);
  if (res1V.GetSize() != res2V.GetSize()) {
    WARN("Input vectors for BinOp Filter are of different size");
    last = std::min(res1V.GetSize(),res2V.GetSize());
  }
  
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(UInt i=0;i<last;++i){
    returnVec[i] = Operator::Apply(res1V[i],res2V[i]);
  }
  
  return true;
}

template<class Operator>
ResultIdList GenericBinOpFilter<Operator>::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  res1Id = upResNameIds[res1Name];
  res2Id = upResNameIds[res2Name];
  outId = filterResIds[0];
  return generated;
}

template<class Operator>
void GenericBinOpFilter<Operator>::AdaptFilterResults(){
  //We should check some validity...
  //we can almost copy everything from time input (also scalar, etc.)
  resultManager_->CopyResultData(res1Id,outId);
  resultManager_->SetValid(outId);
}


// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class GenericBinOpFilter<BinOpFilter::PlusOp>;
  template class GenericBinOpFilter<BinOpFilter::MultOp>;
  template class GenericBinOpFilter<BinOpFilter::DivOp>;
  template class GenericBinOpFilter<BinOpFilter::MinusOp>;
#endif
}
