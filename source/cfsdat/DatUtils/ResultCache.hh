// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ResultCache.hh
 *       \brief    <Description>
 *
 *       \date     Jan 18, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef RESULTCACHE_HH_
#define RESULTCACHE_HH_

#include "General/Environment.hh"
#include "Domain/Results/BaseResults.hh"
#include "cfsdat/DatUtils/Defines.hh"
#include "EqnNumberingSimple.hh"

#include <limits>

namespace CFSDat{

struct GenericResultAdapter{

  GenericResultAdapter(){
    stepIndex = std::numeric_limits<UInt>::max();
    isUpToDate = false;
  }

  virtual ~GenericResultAdapter(){

  }

  virtual SingleVector* GetSingleVector() = 0;

  //shared_ptr< EqnMapSimple > mapping;

  //! result struct for export, one result for each entity list/region
  //CF::StdVector< shared_ptr<CF::BaseResult> > baseResultVector;

  UInt stepIndex;

  bool isUpToDate;
};

template<typename T>
struct ResultAdaptor : public GenericResultAdapter{

  ResultAdaptor(){
  }

  virtual ~ResultAdaptor(){
    resultVector.Clear();
  }

  SingleVector* GetSingleVector(){
    return &resultVector;
  }

  Vector<T> resultVector;
};

typedef shared_ptr<GenericResultAdapter> ResAdaptPtr;

typedef shared_ptr<const GenericResultAdapter> ConstResAdaptPtr;



class GenericResultCache {

public:

  GenericResultCache(UInt cacheSize=1, bool isStatic=false);

  virtual ~GenericResultCache();

  void SetStepIndex(UInt stepIndex);

  UInt GetStepIndex();

  void SetVectorSize(UInt vectorSize);

  UInt GetVectorSize();
  
  void SetUpToDate(bool upToDate);

  bool IsUpToDate();

  virtual SingleVector* GetSingleVector();

  shared_ptr< EqnMapSimple > mapping;

  //! result struct for export, one result for each entity list/region
  CF::StdVector< shared_ptr<CF::BaseResult> > baseResultVector;

protected:

  virtual void InitAdapterVector() = 0; 

  ResAdaptPtr actualAdapter_;

  CF::StdVector<ResAdaptPtr> adapters_;
  
  UInt cacheSize_;
  
  UInt vectorSize_;
  
  bool isStatic_;
  
private:

  bool noCache_;
  
  bool isInit_;
  
  void Init();

};

template<typename T>
class ResultCache : public GenericResultCache {

public:  

  ResultCache(UInt cacheSize=1, bool isStatic=false):
    GenericResultCache(cacheSize,isStatic) {
  }

  virtual ~ResultCache() {
    adapters_.Clear();
  }
  
  CF::Vector<T>& GetResultVector() {
    CF::Vector<T>& resVec = dynamic_cast< ResultAdaptor<T>* >(actualAdapter_.get())->resultVector;
    resVec.Resize(vectorSize_);
    return resVec;
  }
  
protected:

  virtual void InitAdapterVector() {
    adapters_.Resize(cacheSize_);
    for (UInt i = 0; i < cacheSize_; i++) {
      adapters_[i].reset(new ResultAdaptor<T>());
    }
  }
};

}


#endif /* RESULTCACHE_HH_ */
