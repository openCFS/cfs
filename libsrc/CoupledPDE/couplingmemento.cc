#include "couplingmemento.hh"

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
}
