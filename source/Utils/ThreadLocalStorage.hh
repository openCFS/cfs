// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ThreadLocalStorage.hh
 *       \brief    <Description>
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef THREADLOCALSTORAGE_HH_
#define THREADLOCALSTORAGE_HH_

#include "StdVector.hh"
#include "Domain/ElemMapping/Elem.hh"
#include <map>
#include <type_traits>

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField{

//! Interface class for the definition of cloneable Classes
class CfsCopyable{
public:
  CfsCopyable() { } 

  virtual ~CfsCopyable() { }

  //! Method for copying a pointer
  virtual CfsCopyable* Clone()=0;
};




//!  Base class for thread-safe containers in CFS.
/*!
  Basic interface definition. The amount of available slots is
  always the maximum number of OMP threads as specified by the user.
*/
class BaseTLS{
public:
  inline UInt GetNumSlots() const{
    return numSlots_;
  }
protected:
  BaseTLS(){
     numSlots_ = CFS_NUM_THREADS;
  }
  UInt numSlots_;
};

/*! This class stores and administrates thread local copies of
 *! the given type. Addressing the access is done via a Mine function
 *! in future releases, we might consider a Acquire-Release functionality
 *! ASSUMPTION: T has a default and/or a copy constructor
 */
template<class T>
class CfsTLS : public BaseTLS{
public:
  CfsTLS(){
    tlsContainer_.Resize(numSlots_);
  }

  CfsTLS(const T& serialObjToCopy){
    Set(serialObjToCopy);
  }

  ~CfsTLS(){
    tlsContainer_.Clear();
  }

   inline T& Mine(Integer tNum = -1){
#ifdef USE_OPENMP
    // both asserts should test the same, but sometimes tlsContainer_ is not of size numSlots_
    assert((int) omp_get_thread_num() < (int) numSlots_);
    assert((int) omp_get_thread_num() < (int) tlsContainer_.GetSize());
    return (tNum>=0)? tlsContainer_[tNum] : tlsContainer_[omp_get_thread_num()];
#else
    return tlsContainer_[0];
#endif
  }

   inline const T& Mine(Integer tNum = -1) const {
#ifdef USE_OPENMP
    // both asserts should test the same, but sometimes tlsContainer_ is not of size numSlots_
    assert((int) omp_get_thread_num() < (int) numSlots_);
    assert((int) omp_get_thread_num() < (int) tlsContainer_.GetSize());
    return (tNum>=0)? tlsContainer_[tNum] : tlsContainer_[omp_get_thread_num()];
#else
    return tlsContainer_[0];
#endif
  }


   inline const T& ConstMine(Integer tNum = -1) const{
#ifdef USE_OPENMP
    // both asserts should test the same, but sometimes tlsContainer_ is not of size numSlots_
    assert((int) omp_get_thread_num() < (int) numSlots_);
    assert((int) omp_get_thread_num() < (int) tlsContainer_.GetSize());
    return (tNum>=0)? tlsContainer_[tNum] : tlsContainer_[omp_get_thread_num()];
#else
    return tlsContainer_[0];
#endif
  }

  void Set(const T& serialObjToCopy){
    tlsContainer_.Resize(numSlots_,serialObjToCopy);
  }

protected:
  StdVector<T> tlsContainer_;
};


/*! Thread Local Storage Container
 *! This class stores and administrates thread local copies of
 *! the given pointer type. Addressing the access is done via a Mine function
 *! in future releases, we might consider an Acquire-Release functionality
 *! ASSUMPTION: T has a default constructor and a clone method
 */
template<class T>
class CfsTLS<T*> : public BaseTLS{
public:
  CfsTLS(){
    tlsContainer_.Resize(numSlots_);
    for(UInt i=0;i<numSlots_;i++){
      tlsContainer_[i] = NULL;
    }
  }

  CfsTLS(const T* serialObjToCopy){
    tlsContainer_.Resize(numSlots_);
    for(UInt i=0;i<numSlots_;i++){
      tlsContainer_[i] = serialObjToCopy->Clone();
    }
  }

  ~CfsTLS(){
    for(UInt i=0;i<numSlots_;i++){
      if(tlsContainer_[i]){
        delete tlsContainer_[i];
      }
    }
    tlsContainer_.Clear();
  }

  /** distributes the serializable object to all entries */
  void Set(const T* serialObjToCopy){
    for(UInt i=0;i<numSlots_;i++){
      if(tlsContainer_[i]){
        delete tlsContainer_[i];
      }
    }
    tlsContainer_.Resize(numSlots_);
    for(UInt i=0;i<numSlots_;i++){
      tlsContainer_[i] = serialObjToCopy->Clone();
    }
  }

  /** distribute the object via copy constructor to all entries */
  void Set(const T& ref)
  {
    tlsContainer_.Resize(numSlots_);
    for(unsigned int i = 0; i < numSlots_; i++)
      tlsContainer_[i] = ref;
  }

  inline T*& Mine(Integer tNum = -1)
  {
#ifdef USE_OPENMP
     return (tNum>=0)? tlsContainer_[tNum] : tlsContainer_[omp_get_thread_num()];
#else
    return tlsContainer_[0];
#endif
  }

  

private:
  StdVector<T*> tlsContainer_;
};

/*! Thread Local Storage for maps
 *! This class serves as a base class for the thread-local storage of maps.
 *! This map can be initialized with a standard map from a serial region
 *! or from a parallel region. Either the threadlocal storage slot
 *! or all slots will be set.
 */
template<typename K,typename V>
class TLMapBase : public BaseTLS{
public:

  typedef typename std::map<K,V>::iterator tl_iterator;
  typedef typename std::map<K,V>::const_iterator tl_const_iterator;

  TLMapBase(): isCleared_(true){
#ifdef USE_OPENMP
    tlsContainer_.Resize(numSlots_);
#endif
  }

  inline V& operator[](const K& a){
   isCleared_ = false;
#ifdef USE_OPENMP
    return tlsContainer_[omp_get_thread_num()][a];
#else
    return dummyContainer_[a];
#endif
  }

  inline std::map<K,V>& Mine(Integer tNum = -1){
   isCleared_ = false;
#ifdef USE_OPENMP
   return (tNum>=0)? tlsContainer_[tNum] : tlsContainer_[omp_get_thread_num()];
#else
    return dummyContainer_;
#endif
  }

  bool empty() const {
#ifdef USE_OPENMP
    return tlsContainer_[omp_get_thread_num()].empty();
#else
    return dummyContainer_.empty();
#endif
  }

  inline tl_iterator begin(){
   isCleared_ = false;
#ifdef USE_OPENMP
    return tlsContainer_[omp_get_thread_num()].begin();
#else
    return dummyContainer_.begin();
#endif
  }

  inline tl_iterator end(){
   isCleared_ = false;
#ifdef USE_OPENMP
    return tlsContainer_[omp_get_thread_num()].end();
#else
    return dummyContainer_.end();
#endif
  }

  inline tl_iterator find(K key){
   isCleared_ = false;
#ifdef USE_OPENMP
    return tlsContainer_[omp_get_thread_num()].find(key);
#else
    return dummyContainer_.find(key);
#endif
  }


protected:
#ifdef USE_OPENMP
  StdVector< std::map<K,V> > tlsContainer_;
#else
  std::map<K,V> dummyContainer_;
#endif
  //! indicates if the map is supposed to be initialized
  bool isCleared_;
  //NOTE: Initialize the variable with false on construction
  //as we store (if any) only NULL pointers.
  //A call to the [] operator,
  //obtaining any non-const iterator via find,begin,end or a call
  //to Mine() will set the variable to true as it would be possible
  //that memory is allocated.
  //A more consistent class behaviour would require a more restrictive
  //interface of the class which would require more changes in the
  //methods using this class. Currently, this is not what you want to do.
};

/*! General TLS map
 *!
 *! Implementation for maps of objects which have an assignment operator
 */
template<typename K,typename V>
class TLMap : public TLMapBase<K,V>{
public:
  TLMap(){

  }

  ~TLMap(){
    if(!this->isCleared_){
      Clear();
    }
  }

  //! assignment operator to implicitly convert a stl map into the TLS
  TLMap<K,V>& operator=(const std::map<K,V>&  serialObjToCopy){
#ifdef USE_OPENMP
#pragma omp master
    {
    this->tlsContainer_.Resize(this->numSlots_);
    for(UInt i=0;i<this->numSlots_;i++){
      typename TLMapBase<K,V>::tl_iterator mIter = serialObjToCopy.begin();
      for(;mIter!=serialObjToCopy.end();++mIter){
        this->tlsContainer_[i][mIter->first] = mIter->second;
      }
    }
    }
#else
    typename TLMapBase<K,V>::tl_iterator mIter = serialObjToCopy.begin();
    for(;mIter!=serialObjToCopy.end();++mIter){
      this->dummyContainer_[mIter->first] = mIter->second;
    }
#endif
    this->isCleared_ = false;
    return *this;
  }

  void Clear(){
    //if any thread than the master will clear that,
    // we will have a problem. Still, in our current setup this
    //should not happen and if, the master will also be here...
#ifdef USE_OPENMP
#pragma omp master
    {
    this->tlsContainer_.Clear();
    this->isCleared_ = true;
    }
#else
    this->dummyContainer_.clear();
    this->isCleared_ = true;
#endif
  }
};

/*! TLS map of pointers
 *!
 *! Implementation for maps which store pointers to objects as a value.
 *! Thereby, on Clear the pointer values will be deleted for all threads
 *! unless called from a parallel region
 */
template<typename K,typename V>
class TLMap<K,V*> : public TLMapBase<K,V*>{

public:

  TLMap(){
  }

  ~TLMap(){
    if(!this->isCleared_){
      Clear();
    }
  }

  //! assignment operator to implicitly convert a stl map into the TLS
  TLMap<K,V*>& operator=(const std::map<K,V*>&  serialObjToCopy){
#ifdef USE_OPENMP
    //complete assignments should only occur from serial regions...
#pragma omp master
    {
    this->tlsContainer_.Resize(this->numSlots_);
    for(UInt i=0;i<this->numSlots_;i++){
      typename TLMapBase<K,V*>::tl_const_iterator mIter = serialObjToCopy.begin();
      for(;mIter!=serialObjToCopy.end();++mIter){
        this->tlsContainer_[i][mIter->first] = mIter->second->Clone();
      }
    }
    }
#else
    typename TLMapBase<K,V*>::tl_iterator mIter = serialObjToCopy.begin();
    for(;mIter!=serialObjToCopy.end();++mIter){
      this->dummyContainer_[mIter->first] = mIter->second->Clone();
    }
#endif
    this->isCleared_ = false;
    return *this;
  }

  void Clear(){
    //if any thread than the master will clear that,
    // we will have a problem. Still, in our current setup this
    //should not happen and if, the master will also be here...
#ifdef USE_OPENMP
#pragma omp master
    {
    for(UInt i=0;i<this->numSlots_;i++){
      typename TLMapBase<K,V*>::tl_iterator mIter = this->tlsContainer_[i].begin();
      for(;mIter!=this->tlsContainer_[i].end();++mIter){
        if(mIter->second)
          delete mIter->second;
      }
    }
    this->tlsContainer_.Clear();
    this->isCleared_ = true;
    }
#else
    typename TLMapBase<K,V*>::tl_iterator mIter = this->dummyContainer_.begin();
    for(;mIter!=this->dummyContainer_.end();++mIter){
      if(mIter->second)
        delete mIter->second;
    }
    this->dummyContainer_.clear();
    this->isCleared_ = true;
#endif
  }
};

}



#endif /* THREADLOCALSTORAGE_HH_ */
