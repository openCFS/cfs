// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "boost/serialization/export.hpp"

#include "vector.cc"
#include "vectorSerialization.hh"

namespace CoupledField {
class SingleVector;
template <class TYPE> class Vector;
}  // namespace CoupledField

BOOST_CLASS_EXPORT_GUID(CoupledField::SingleVector, "CoupledField_SingleVector")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Double>, "CoupledField_Vector_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Complex>, "CoupledField_Vector_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Integer>, "CoupledField_Vector_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::UInt>, "CoupledField_Vector_UInt")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<bool>, "CoupledField_Vector_Bool")

