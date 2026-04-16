#include <algorithm>
#include "vectorSerialization.hh"
#include "Vector.cc"
#include "Utils/boost-serialization.hh"
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT_GUID(CoupledField::SingleVector, "CoupledField_SingleVector")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Double>, "CoupledField_Vector_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Complex>, "CoupledField_Vector_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::Integer>, "CoupledField_Vector_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Vector<CoupledField::UInt>, "CoupledField_Vector_UInt")

