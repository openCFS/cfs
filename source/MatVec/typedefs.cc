// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "MatVec/typedefs.hh"

namespace CoupledField {

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
// Scalar matrices and vectors
template<> std::string AssocType<Double >::tagM = "Double";
template<> std::string AssocType<Double >::tagV = "Double";
template<> std::string AssocType<Double >::tagS = "Double";
template<> std::string AssocType<Complex>::tagM = "Complex";
template<> std::string AssocType<Complex>::tagV = "Complex";
template<> std::string AssocType<Complex>::tagS = "Complex";
template<> std::string AssocType<Integer>::tagM = "Integer";
template<> std::string AssocType<Integer>::tagV = "Integer";
template<> std::string AssocType<Integer>::tagS = "Integer";
template<> std::string AssocType<UInt>::tagM = "UInt";
template<> std::string AssocType<UInt>::tagV = "UInt";
template<> std::string AssocType<UInt>::tagS = "UInt";
template<> std::string AssocType<bool>::tagM = "bool";
template<> std::string AssocType<bool>::tagV = "bool";
template<> std::string AssocType<bool>::tagS = "bool";
#endif

}

