// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMSHHELPER_2009
#define FILE_GMSHHELPER_2009

#include <map>

#include <boost/detail/endian.hpp>
#include <boost/bimap.hpp>

#include "General/Environment.hh"
#include "Domain/ElemMapping/Elem.hh"

#include "def_config.hh"

#if defined(__INTEL_COMPILER)
#include "immintrin.h"
#endif

#if defined(__GNUC__) && not defined(__INTEL_COMPILER)
#ifdef HAVE_BYTESWAP_H
#include <byteswap.h>
#endif
#endif

namespace CoupledField {

  // Solution for endian conversion in C++ from:
  // http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c

  // We define some constant for little, big and host endianess. Here I use
  // BOOST_LITTLE_ENDIAN/BOOST_BIG_ENDIAN to check the host indianess. If you
  // don't want to use boost you will have to modify this part a bit.
  enum EEndian
  {
    LITTLE_ENDIAN_ORDER,
    BIG_ENDIAN_ORDER,
 #if defined(BOOST_LITTLE_ENDIAN)
    HOST_ENDIAN_ORDER = LITTLE_ENDIAN_ORDER
 #elif defined(BOOST_BIG_ENDIAN)
    HOST_ENDIAN_ORDER = BIG_ENDIAN_ORDER
 #else
 #error "Impossible to determine endianness of system."
 #endif
  };

  // this function swap the bytes of values given it's size as a template
  // parameter (could sizeof be used?).
  template <class T, unsigned int size>
  inline T SwapBytes(T value)
  {
    union
    {
      T value;
      char bytes[size];
    } in, out;
    
    in.value = value;

#if defined(__GNUC__) && not defined(__INTEL_COMPILER)
#define GCC_VERSION (__GNUC__ * 10000 \
                     + __GNUC_MINOR__ * 100 \
                     + __GNUC_PATCHLEVEL__)

#ifdef HAVE_BYTESWAP_H
    out = in;

    switch(size)
    {
      case 2:
        bswap_16(out.value);
        break;
      case 4:
        bswap_32(out.value);
        break;
      case 8:
        {
          unsigned long long int* t = (unsigned long long int*)&out.value;
          bswap_64(*t);
        }
        break;
      default:
        break;
    }

#else // HAVE_BYTESWAP_H

    switch(size)
    {
/* Test for GCC > 4.1.0 */
#if GCC_VERSION >= 40200
      case 4:
        out.value = __builtin_bswap32(in.value);
        break;
#endif
      default:
        for (unsigned int i = 0; i < size / 2; ++i)
        {
          out.bytes[i] = in.bytes[size - 1 - i];
          out.bytes[size - 1 - i] = in.bytes[i];
        }
        break;
    }

#endif // HAVE_BYTESWAP_H    

#elif defined(__INTEL_COMPILER)
#if 0
    std::cout << "_bswap 1: " << _bswap(1) << std::endl;
    in.value = 1;
    for (unsigned int i = 0; i < size / 2; ++i)
    {
      out.bytes[i] = in.bytes[size - 1 - i];
      out.bytes[size - 1 - i] = in.bytes[i];
    }
    std::cout << "1: " << out.value << std::endl;
#endif

    // Intel compiler
    switch(size)
    {
      case 4:
        *((int*)&out.value) = _bswap(*((int*)&in.value));
        break;
      case 8:
        *((__int64*)&out.value) = _bswap64(*((__int64*)&in.value));
        break;
      default:
        for (unsigned int i = 0; i < size / 2; ++i)
        {
          out.bytes[i] = in.bytes[size - 1 - i];
          out.bytes[size - 1 - i] = in.bytes[i];
        }
        break;
    }
#else
    for (unsigned int i = 0; i < size / 2; ++i)
    {
      out.bytes[i] = in.bytes[size - 1 - i];
      out.bytes[size - 1 - i] = in.bytes[i];
    }
#endif

    return out.value;
  }
  
  // Here is the function you will use. Again there is two compile-time assertion
  // that use the boost librarie. You could probably comment them out, but if you
  // do be cautious not to use this function for anything else than integers
  // types. This function need to be calles like this :
  //
  //     int x = someValue;
  //     int i = EndianSwapBytes<HOST_ENDIAN_ORDER, BIG_ENDIAN_ORDER>(x);
  //
  
  template<class T>
  struct EndianSwapper 
  {
    bool swap;
    EndianSwapper() : swap(false) {};
    
    inline T EndianSwapBytes(T value)
    {
      // A : La donnee a swapper a une taille de 2, 4 ou 8 octets
      BOOST_STATIC_ASSERT(sizeof(T) == 2 || sizeof(T) == 4 || sizeof(T) == 8);
      // A : La donnee a swapper est d'un type arithmetic
      BOOST_STATIC_ASSERT(boost::is_arithmetic<T>::value);
      // Si from et to sont du meme type on ne swap pas.
      if (!swap)
      {
        return value;
      }
      
      return SwapBytes<T, sizeof(T)>(value);
    }
  };

  class GmshHelper
  {
  public:
    GmshHelper();
    virtual ~GmshHelper();
    
    static void ElemType2FEType(UInt eType, Elem::FEType& feType, UInt& numNodes);
    static void FEType2ElemType(Elem::FEType feType, UInt& eType, UInt& numNodes);

    static void InitElemNodeMap();
    
    struct NodeStruct 
    {
      UInt nodeId;
      Double x,y,z;
    };
    
    typedef boost::bimap< UInt, UInt > ElemNodeMapType;
    static std::map<Elem::FEType, ElemNodeMapType> elemNodeMap_;
  };
} 

#endif
