// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include "boost/serialization/export.hpp"

#include "matrixSerialization.hh"
#include "matrix.hh" // IWYU pragma: keep

BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Double>, "CoupledField_Matrix_Double")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Complex>, "CoupledField_Matrix_Complex")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::Integer>, "CoupledField_Matrix_Integer")
BOOST_CLASS_EXPORT_GUID(CoupledField::Matrix<CoupledField::UInt>,    "CoupledField_Matrix_UInt" )
