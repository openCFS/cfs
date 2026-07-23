// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     SyncAccess.hh
 *       \brief    Implementation structures for syncronized data access functions
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SYNCACCESS_HH_
#define SYNCACCESS_HH_

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#define SYNC_DATA true
#else
#define SYNC_DATA false
#endif

namespace CoupledField {

  //! Synchronize access functions to data entries
  //! we need the specialization in case we have openMP not activated
  //! in future we might want to change the critical stuff to use locks
  //! which are more favorable
  template<bool S>
  struct SyncAccess{
  public:
    static inline void AddTo(Double & target   , const Double  & src);
    static inline void AddTo(Complex & target  , const Complex & src);
    static inline void Set(Double & target     , const Double  & src);
    static inline void Set(Complex & target    , const Complex & src);
  };

  template<>
  struct SyncAccess<true>{
  public:
    static inline void AddTo(Double & target   , const Double  & src){
      #pragma omp atomic
        target += src;
    }
    static inline void AddTo(Complex & target  , const Complex & src){
      // std::complex<Double> is guaranteed to be layout-compatible with
      // Double[2], so the real and imaginary parts can be updated with two
      // independent atomic adds. This avoids the former global critical
      // section, which serialized all threads on every complex-valued update.
      Double* t = reinterpret_cast<Double*>(&target);
      const Double re = src.real();
      const Double im = src.imag();
      #pragma omp atomic
        t[0] += re;
      #pragma omp atomic
        t[1] += im;
    }
    static inline void Set(Double & target   , const Double  & src){
      #pragma omp atomic write
        target = src;
    }
    static inline void Set(Complex & target  , const Complex & src){
      // see AddTo: component-wise atomic writes instead of a critical section
      Double* t = reinterpret_cast<Double*>(&target);
      const Double re = src.real();
      const Double im = src.imag();
      #pragma omp atomic write
        t[0] = re;
      #pragma omp atomic write
        t[1] = im;
    }
  };

  template<>
  struct SyncAccess<false>{
  public:
    static inline void AddTo(Double & target   , const Double  & src){
      target += src;
    }
    static inline void AddTo(Complex & target  , const Complex & src){
      target += src;
    }
    static inline void Set(Double & target   , const Double  & src){
      target = src;
    }
    static inline void Set(Complex & target  , const Complex & src){
      target = src;
    }
  };
}

#endif

