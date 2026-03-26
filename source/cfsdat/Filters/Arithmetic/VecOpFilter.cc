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
GenericVecOpFilter<Operator>::GenericVecOpFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
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
bool GenericVecOpFilter<Operator>::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  Integer stepIndex = resultManager_->GetStepIndex(outId);

  Vector<Double>& res1V = GetUpstreamResultVector<Double>(res1Id, stepIndex);
  Vector<Double>& res2V = GetUpstreamResultVector<Double>(res2Id, stepIndex);
  
  const UInt last = returnVec.GetSize();
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(UInt i=0;i<last;++i){
    returnVec[i] = Operator::Apply(res1V[i],res2V[i]);
  }
  
  return true;
}

template<class Operator>
ResultIdList GenericVecOpFilter<Operator>::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  res1Id = upResNameIds[res1Name];
  res2Id = upResNameIds[res2Name];
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
template class GenericVecOpFilter<VecOpFilter::PlusOp>;
template class GenericVecOpFilter<VecOpFilter::MultOp>;
template class GenericVecOpFilter<VecOpFilter::DivOp>;
template class GenericVecOpFilter<VecOpFilter::MinusOp>;
}
