#include "matrixSerialization.hh"
#include <boost/serialization/export.hpp>

#ifndef USE_ADOLC
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Double>, "CoupledField_Matrix_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Complex>, "CoupledField_Matrix_Complex")
#endif
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Integer>, "CoupledField_Matrix_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::UInt>,    "CoupledField_Matrix_UInt" )
