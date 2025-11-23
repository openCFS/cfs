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

GenericResultCache::GenericResultCache(UInt cacheSize, bool isStatic) {
  vectorSize_ = 0;
  isStatic_ = isStatic;
  cacheSize_ = isStatic ? 1 : cacheSize;
  noCache_ = cacheSize_ == 1;
  isInit_ = false;
}

GenericResultCache::~GenericResultCache() {
  adapters_.Clear();
}

void GenericResultCache::SetVectorSize(UInt vectorSize) {
  vectorSize_ = vectorSize;
  isInit_ = false;
}

UInt GenericResultCache::GetVectorSize() {
  return vectorSize_;
}

void GenericResultCache::SetStepIndex(UInt stepIndex) {
  Init();
  // value is alright already
  if (stepIndex == actualAdapter_->stepIndex) {
    return;
  }
  
  // only one cached value exists, so it is overwritten
  if (noCache_) {
    actualAdapter_->stepIndex = stepIndex;
    if (!isStatic_) {
      actualAdapter_->isUpToDate = false;
    }
    return;
  }

  // looking if wanted value is still in cache
  for(UInt i = 0; i < cacheSize_; i++) {
    if (stepIndex == adapters_[i]->stepIndex) {
      actualAdapter_ = adapters_[i];
      return;
    }
  }

  // replacing old cached value by new value
  UInt maxDistance = 0;
  UInt replaceIndex = 0;
  for(UInt i = 0; i < cacheSize_; i++) {
    UInt iDistance = 0;
    if (adapters_[i]->stepIndex > stepIndex) {
      iDistance = adapters_[i]->stepIndex - stepIndex;
    } else {
      iDistance = stepIndex - adapters_[i]->stepIndex;
    }
    if (iDistance > maxDistance) {
      maxDistance = iDistance;
      replaceIndex = i;
    }
  }
  
  actualAdapter_ = adapters_[replaceIndex];
  actualAdapter_->stepIndex = stepIndex;
  actualAdapter_->isUpToDate = false;
}

UInt GenericResultCache::GetStepIndex() {
  Init();
  return actualAdapter_->stepIndex;
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

