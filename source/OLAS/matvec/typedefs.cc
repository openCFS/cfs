// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "matvec/typedefs.hh"

namespace OLAS {

  // Scalar matrices and vectors
template<> std::string assocType<Double >::tagM = "Double";
template<> std::string assocType<Double >::tagV = "Double";
template<> std::string assocType<Double >::tagS = "Double";
template<> std::string assocType<Complex>::tagM = "Complex";
template<> std::string assocType<Complex>::tagV = "Complex";
template<> std::string assocType<Complex>::tagS = "Complex";


}

