// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "pdememento.hh"

namespace CoupledField{

  PDEMemento::PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::PDEMemento");

    isSet_ = false;
    stepNum_ = 0;
    isIterCoupled_ = false;
  }

  PDEMemento::~PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::~PDEMemento");
    solution_.clear();
  }

  void PDEMemento::Clear()
  {
    ENTER_FCN( "PDEMemento::Clear");
    std::map<std::string,CFSVector*>::iterator it;
    
    for( it = solution_.begin(); it != solution_.end(); it++ ) {
      delete it->second;
    }
    solution_.clear();
  }

  std::string PDEMemento::Enum2String( ValueUsageType type ) {

    std::string ret;
    switch( type ) {
    case NO_TYPE:
      ret = "noType";
      break;
    case START_VALUE :
      ret = "startValue";
      break;
    case DIRICHLET_VALUE:
      ret = "dirichletValue";
      break;
    default:
      *error << "ValueUsageType " << type << " not known!";
      Error( __FILE__, __LINE__ );
    }

    return ret;
  }

  PDEMemento::ValueUsageType PDEMemento::String2Enum( std::string typeString ) {
    ValueUsageType ret = NO_TYPE;
    
    if( typeString == "noType" ) {
      ret = NO_TYPE;
    } else if ( typeString == "startValue" ) {
      ret = START_VALUE;
    } else if ( typeString == "dirichletValue" ) {
      ret = DIRICHLET_VALUE;
    } else {
      *error << "ValueUsageType '" << typeString << "' not known!";
      Error( __FILE__, __LINE__ );
    }

    return ret;
  }
  
  
  


} // namespace

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::PDEMemento, "CoupledField_PDEMemento")

