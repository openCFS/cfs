#include "couplingmemento.hh"

#include "Utils/boost-serialization.hh"

namespace CoupledField{

  CouplingMemento::CouplingMemento()
  {
    ENTER_IFCN( "CouplingMemento::CouplingMemento");
  
    isSet_ = false;
  
  }

  CouplingMemento::~CouplingMemento()
  {
    ENTER_IFCN( "CouplingMemento::~CouplingMemento");
  
  }


  void CouplingMemento::Clear()
  {
    ENTER_FCN( "CouplingMemento::Clear");
  
    inputTypes_.Clear();
    inputQuantities_.Clear();
    inputInterfaces_.Clear();

  }

  template<class Archive>
  void CouplingMemento::serialize(Archive & ar, 
                                  const unsigned int version) {
    Error( "Not imeplemented at the momement", __FILE__, __LINE__ );
    
    // For further details, why this is currently not imeplemented,
    // see PDECoupling::CouplingInterface::serialize
  }
}

#include <boost/serialization/export.hpp>
BOOST_CLASS_EXPORT_GUID(CoupledField::CouplingMemento, 
                        "CoupledField_CouplingMemento")

