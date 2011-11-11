// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "CouplingMemento.hh"

namespace CoupledField{

  CouplingMemento::CouplingMemento()
  {
  
    isSet_ = false;
  
  }

  CouplingMemento::~CouplingMemento()
  {
  
  }


  void CouplingMemento::Clear()
  {
  
    inputTypes_.Clear();
    inputQuantities_.Clear();
    inputInterfaces_.Clear();

  }

}

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::CouplingMemento, 
                        "CoupledField_CouplingMemento")

