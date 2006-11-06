#include "pdememento.hh"

#include "Utils/boost-serialization.hh"

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
  
  
  
  template <class Archive>
  void PDEMemento::serialize(Archive & ar, const unsigned int version) {

    ar & isSet_;
    ar & analysisType_;
    ar & gridFileName_;
    ar & stepNum_;
    ar & freq_;
    ar & solution_;
    ar & solDeriv1_;
    ar & solDeriv2_;
    ar & isIterCoupled_;
    if( isIterCoupled_ ) {
      ar & couplingMemento_;
    }

  }

} // namespace

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::PDEMemento, "CoupledField_PDEMemento")

