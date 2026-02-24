// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BinOpFilter.cc
 *       \brief    <Description>
 *
 *       \date     Mai 13, 2018
 *       \author   mtautz
 */
//================================================================================================


#include "BinOpFilter.hh"
#include "BinOpFunctions.hh"

#include <algorithm> 

namespace CFSDat{

BinOpFilter::BinOpFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                             :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  
  opType = config->Get("opType")->As<std::string>();
  multType = config->Get("multType")->As<std::string>();
  
  if (config->Get("quantity1")->Has("resultName")) {
    this->resAName = config->Get("quantity1")->Get("resultName")->As<std::string>();
    upResNames.insert(resAName);
    isConstantA = false;  
  } else if (config->Get("quantity1")->Has("value")) {
    constantA.Resize(1);
    constantA[0] = config->Get("quantity1")->Get("value")->As<double>();
    isConstantA = true;  
  } else {
    EXCEPTION("Either resultName or value should be given for quantity1");
  }
  
  if (config->Get("quantity2")->Has("resultName")) {
    this->resBName = config->Get("quantity2")->Get("resultName")->As<std::string>();
    upResNames.insert(resBName);
    isConstantB = false;  
  } else if (config->Get("quantity2")->Has("value")) {
    constantB.Resize(1);
    constantB[0] = config->Get("quantity2")->Get("value")->As<double>();
    isConstantB = true;  
  } else {
    EXCEPTION("Either resultName or value should be given for quantity2");
  }
  
  if (isConstantA && isConstantB) {
    EXCEPTION("At least one argeument for the BinOpFilter should not be a constant");
  }
  
  this->outName = config->Get("output")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
}

BinOpFilter::~BinOpFilter() {
  
}

bool BinOpFilter::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  Integer stepIndex = resultManager_->GetStepIndex(outId);

  Vector<Double>& resAV = isConstantA ? constantA : GetUpstreamResultVector<Double>(resAId, stepIndex);
  Vector<Double>& resBV = isConstantB ? constantB : GetUpstreamResultVector<Double>(resBId, stepIndex);
  
  applyFctPtr(returnVec, resAV, resBV, size);
  
  
  return true;
}

ResultIdList BinOpFilter::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  if (!isConstantA) {
    resAId = upResNameIds[resAName];
  }
  if (!isConstantB) {
    resBId = upResNameIds[resBName];
  }
  outId = filterResIds[0];
  return generated;
}

void BinOpFilter::AdaptFilterResults(){
  // copy result data from given (non consant and possibly non static) input to output
  bool isStaticA = isConstantA || resultManager_->IsStatic(resAId);
  bool isStaticB = isConstantB || resultManager_->IsStatic(resBId);
  if (!isStaticA) {
    resultManager_->CopyResultData(resAId,outId);  
  } else if (!isStaticB) {
    resultManager_->CopyResultData(resBId,outId);  
  } else if (!isConstantA) {
    resultManager_->CopyResultData(resAId,outId);  
  } else {
    resultManager_->CopyResultData(resBId,outId);  
  }
  resultManager_->SetStatic(outId, isStaticA && isStaticB); 
  
  // Getting Operator Function
  BinOpFunctions::Operand<Double> opA(resultManager_.get(), resAId, constantA);
  BinOpFunctions::Operand<Double> opB(resultManager_.get(), resBId, constantB);
  Vector<Double> dummyConstant;
  BinOpFunctions::Operand<Double> opOut(resultManager_.get(), outId, dummyConstant);
  applyFctPtr = BinOpFunctions::GetBinOpFunction<Double>(opType, multType, opOut, opA, opB);
  
  resultManager_->SetValid(outId);
}

void BinOpFilter::PrepareCalculation() {
  size = resultManager_->GetEqnMap(outId)->GetNumEntities();
}

}
