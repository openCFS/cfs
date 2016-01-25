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

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField{

//! Interface class for the definition of cloneable Classes
class CfsCopyable{
public:
  CfsCopyable(){

  }

  virtual ~CfsCopyable(){

  }

  //! Method for copying a pointer
  virtual CfsCopyable* Clone()=0;
};




//!  Base class for thread-safe containers in CFS.
/*!
  Basic interface definition. The amount of available slots is
  always the maximum number of OMP threads.
  Needs to be adapted in the future to have a variable?
*/
class BaseTLS{
public:
  UInt GetNumSlots(){
    return numSlots_;
  }
protected:
  BaseTLS(){
     numSlots_ = NUM_CFS_THREADS;
  }
  UInt numSlots_;
};

/*! This class stores and administrates thread local copies of
 *! the given type. Addessing the access is done via a Mine function
 *! in future releases, we might consider a Aqquire-Release functionality
 *! ASSUMPTION: T has a default and/or a copy constructor
 */
template<class T>
class CfsTLS : public BaseTLS{
public:
  CfsTLS(){
    tlsContainer_.Resize(numSlots_);
  }

  CfsTLS(const T& serialObjToCopy){
    tlsContainer_.Resize(numSlots_,serialObjToCopy);
  }

  ~CfsTLS(){
    tlsContainer_.Clear();
  }

  T& Mine(Integer tNum = -1){
#ifdef USE_OPENMP
    if(tNum>=0){
      return tlsContainer_[tNum];
    }else{
      return tlsContainer_[omp_get_thread_num()];
    }
#else
    return tlsContainer_[0];
#endif
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
      tlsContainer_[i] = new T();
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

  T* Mine(Integer tNum = -1){
#ifdef USE_OPENMP
    if(tNum>=0){
      return tlsContainer_[tNum];
    }else{
      return tlsContainer_[omp_get_thread_num()];
    }
#else
    return tlsContainer_[0];
#endif
  }
private:
  StdVector<T*> tlsContainer_;
};


}



#endif /* THREADLOCALSTORAGE_HH_ */
