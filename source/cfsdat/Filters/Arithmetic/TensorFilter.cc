// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TensorFilter.cc
 *       \brief    <Description>
 *
 *       \date     Apr 13, 2016
 *       \author   ahueppe
 */
//================================================================================================


#include "TensorFilter.hh"

#include <algorithm> 

namespace CFSDat{

TensorFilter::TensorFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                             :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  
  resNames.Resize(9);
  for(UInt i = 0; i < 9;i++) {
    resNames[i] = config->Get("quantity" + std::to_string(i+1))->Get("resultName")->As<std::string>();
    upResNames.insert(resNames[i]);
  }
  
  this->outName = config->Get("output")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
}

TensorFilter::~TensorFilter() {
  
}

bool TensorFilter::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  Integer stepIndex = resultManager_->GetStepIndex(outId);
  
  for(UInt i = 0; i < 9;i++) {
    Vector<Double>& upResVec = GetUpstreamResultVector<Double>(resIds[i], stepIndex);
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt j = 0; j < size; j++) {
      returnVec[j * 9 + i] = upResVec[j];
    }
  }
  return true;
}

ResultIdList TensorFilter::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  resIds.Resize(9);
  for(UInt i = 0; i < 9;i++) {
    resIds[i] = upResNameIds[resNames[i]];
  }
  outId = filterResIds[0];
  return generated;
}

void TensorFilter::AdaptFilterResults(){
  bool isStatic = true;
  UInt copyResultIdx = 0;
  for(UInt i = 0; i < 9;i++) {
    isStatic = isStatic && resultManager_->IsStatic(resIds[i]);
    if (!isStatic) {
      copyResultIdx = i;
      break;
    }
  }
  
  resultManager_->CopyResultData(resIds[copyResultIdx],outId);
  StdVector<std::string> newDofNames;
  newDofNames.Resize(9);
  newDofNames[0] = "xx";
  newDofNames[1] = "xy";
  newDofNames[2] = "xz";
  newDofNames[3] = "yx";
  newDofNames[4] = "yy";
  newDofNames[5] = "yz";
  newDofNames[6] = "zx";
  newDofNames[7] = "zy";
  newDofNames[8] = "zz";
  resultManager_->SetDofNames(outId, newDofNames);
  resultManager_->SetEntryType(outId, ExtendedResultInfo::TENSOR);

  resultManager_->SetValid(outId);
}

void TensorFilter::PrepareCalculation() {
  size = resultManager_->GetEqnMap(outId)->GetNumEntities();
}

}
