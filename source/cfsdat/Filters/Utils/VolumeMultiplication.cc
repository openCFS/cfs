// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VolumeMultiplication.cc
 *       \brief    Filter to multiply values given on elements by the elements volume. This does not work for values given on nodes.
 *
 *       \date     Apr 23, 2018
 *       \author   mtautz
 */
//================================================================================================

#include "VolumeMultiplication.hh"

namespace CFSDat{
  
VolumeMultiplication::VolumeMultiplication(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
: BaseFilter(numWorkers,config, resMan){

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  this->resName = config->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  this->outName = config->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
  upResNames.insert(resName);

}

VolumeMultiplication::~VolumeMultiplication(){
}

void VolumeMultiplication::PrepareCalculation() {
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(resId);
  Grid* grid = inInfo->ptGrid;
  
  str1::shared_ptr<EqnMapSimple> scrMap = resultManager_->GetEqnMap(resId);
  StdVector<UInt> revMap = resultManager_->GetEntityNumbers(resId);
  scrMap->GetReverseEntityMap(revMap);
  usedElems_ = scrMap->GetNumEntities();
  volume_.Resize(usedElems_);
  volume_.Init(0.0);
  
  #pragma omp parallel for num_threads(CFS_NUM_THREADS) 
  for (UInt i = 0; i < usedElems_; i++) {
    UInt elemNum = revMap[i];
    const Elem* elem = grid->GetElem(elemNum);
    // shared_ptr<ElemShapeMap> esm = grid->GetElemShapeMap(elem);
    // grid->UpdateIndividualElemShapeMap(elem, true); // const
    // shared_ptr<ElemShapeMap> esm(elem->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = ((elem))->GetElemShapeMap(grid);
    Elem::ShapeType st = Elem::GetShapeType(elem->type);
    ElemShape& eShape = Elem::GetShape(st);
    if (eShape.dim == 3) {
      volume_[i] = esm->CalcVolume();
    }
  }
  
  CF::StdVector<std::string> dofNames = resultManager_->GetDofNames(resId);
  numDofs_ = dofNames.GetSize();
}

bool VolumeMultiplication::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  Integer stepIndex = resultManager_->GetStepIndex(outId);
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(resId, stepIndex);
  
  if (numDofs_ > 0) {
    #pragma omp parallel for num_threads(CFS_NUM_THREADS) 
    for (UInt i = 0; i < usedElems_; i++) {
      Double factor = volume_[i];
      Double* iIn = &inVec[i * numDofs_];
      Double* iOut = &returnVec[i * numDofs_];
      for (UInt d = 0; d < numDofs_; d++) {
        iOut[d] = factor * iIn[d];
      }
    }
  } else {
    #pragma omp parallel for num_threads(CFS_NUM_THREADS) 
    for (UInt i = 0; i < usedElems_; i++) {
      returnVec[i] =  volume_[i] * inVec[i];
    }
  }
  
  //VerboseSum(outId);
  
  return true;
}

ResultIdList VolumeMultiplication::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  resId = upResNameIds[*upResNames.begin()];
  outId = filterResIds[0];
  return generated;
}

void VolumeMultiplication::AdaptFilterResults() {
  //we would need to copy related data to our filter results
  // getting some basic information
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(resId);
  if (inInfo->definedOn != ExtendedResultInfo::ELEMENT) {
    EXCEPTION("Volume multiplication only works with data provided on elements");
  }
  Grid* grid = inInfo->ptGrid;
  if (grid->GetDim() != 3) {
    EXCEPTION("Volume multiplication only works with three dimensional meshes");
  }
  
  resultManager_->CopyResultData(resId,outId);
  resultManager_->SetValid(outId);
}

}

