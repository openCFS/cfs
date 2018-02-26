// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ResultCache.cc
 *       \brief    <Description>
 *
 *       \date     Jan 18, 2015
 *       \author   mtautz
 */
//================================================================================================

#include "ResultCache.hh"
#include "MatVec/Vector.hh"

namespace CFSDat {

bool EqualStepValues(Double valueA, Double valueB) {
  return valueA == valueB;
}

GenericResultCache::GenericResultCache(UInt cacheSize) {
  vectorSize_ = 0;
  cacheSize_ = cacheSize;
  noCache_ = cacheSize_ == 1;
  isInit_ = false;
}

GenericResultCache::~GenericResultCache() {
  adapters_.Clear();
}

void GenericResultCache::SetVectorSize(UInt vectorSize) {
  vectorSize_ = vectorSize;
}

UInt GenericResultCache::GetVectorSize() {
  return vectorSize_;
}

void GenericResultCache::SetStepValue(Double stepValue) {
  Init();
  // value is alright already
  if (EqualStepValues(stepValue, actualAdapter_->stepValue)) {
    return;
  }
  
  // only one cached value exists, so it is overwritten
  if (noCache_) {
    actualAdapter_->stepValue = stepValue;
    actualAdapter_->isUpToDate = false;
    return;
  }

  // looking if wanted value is still in cache
  UInt i = 0;
  for(; i < cacheSize_; i++) {
    if (EqualStepValues(stepValue, adapters_[i]->stepValue)) {
      actualAdapter_ = adapters_[i];
      return;
    }
  }

  // replacing old cached value by new value
  Double maxDistance = abs(adapters_[0]->stepValue - stepValue);
  uint replaceIndex = 0;
  i = 1;
  for(; i < cacheSize_; i++) {
    Double distance = abs(adapters_[i]->stepValue - stepValue);
    if (distance > maxDistance) {
      maxDistance = distance;
      replaceIndex = i;
    }
  }
  
  actualAdapter_ = adapters_[replaceIndex];
  actualAdapter_->stepValue = stepValue;
  actualAdapter_->isUpToDate = false;
}

Double GenericResultCache::GetStepValue() {
  Init();
  return actualAdapter_->stepValue;
}

void GenericResultCache::SetUpToDate(bool upToDate) {
  Init();
  actualAdapter_->isUpToDate = upToDate;
}

bool GenericResultCache::IsUpToDate() {
  Init();
  return actualAdapter_->isUpToDate;
}

SingleVector* GenericResultCache::GetSingleVector() {
  Init();
  return actualAdapter_->GetSingleVector();
}

void GenericResultCache::Init() {
  if (isInit_) {
    return;
  }
  InitAdapterVector();
  actualAdapter_= adapters_[0];
  isInit_ = true;
}


}

