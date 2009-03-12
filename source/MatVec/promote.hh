// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef PROMOTE_HH
#define PROMOTE_HH

#include "General/defs.hh"
using namespace CoupledField;
template<typename T1, typename T2> 
struct Trait2 {
};

#define DECLARE_TRAIT2(A,B,C)                    \
  template<> struct Trait2<A,B> {                 \
    typedef C T_PROMOTE;\
};
DECLARE_TRAIT2(Complex,Double,Complex)
DECLARE_TRAIT2(Double,Complex,Complex)
DECLARE_TRAIT2(Double,Double,Double)
DECLARE_TRAIT2(Complex,Complex,Complex)

// DEFINE for type promotion
#define PROMOTE(P,P2) typename Trait2<P,P2>::T_PROMOTE

#endif
