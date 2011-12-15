// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <utility>

#include "MatVec/SingleVector.hh"
#include "MatVec/vectorSerialization.hh" // IWYU pragma: keep
#include "pdememento.hh"


namespace CoupledField{

  PDEMemento::PDEMemento()
  {

    isSet_ = false;
    stepNum_ = 0;
    isIterCoupled_ = false;
  }

  PDEMemento::~PDEMemento()
  {
    Clear();
  }

  void PDEMemento::Clear()
  {
    std::map<std::string,SingleVector*>::iterator it;
    
    for( it = solution_.begin(); it != solution_.end(); it++ )
      delete it->second;
    
    solution_.clear();
  }


} // namespace

#include "boost/serialization/export.hpp"

BOOST_CLASS_EXPORT_GUID(CoupledField::PDEMemento, "CoupledField_PDEMemento")

