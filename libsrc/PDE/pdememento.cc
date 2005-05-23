#include "pdememento.hh"

namespace CoupledField{

  PDEMemento::PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::PDEMemento");
  
    sol_       = NULL;
    isSet_ = FALSE;
  }

  PDEMemento::~PDEMemento()
  {
    ENTER_IFCN( "PDEMemento::~PDEMemento");
  
    if (sol_)
      delete sol_;

    sol_       = NULL;
  }

  void PDEMemento::Clear()
  {
    ENTER_FCN( "PDEMemento::Clear");
    if (sol_)
      delete sol_;

    sol_       = NULL;
    solDeriv1_.Clear();
    solDeriv2_.Clear();
  }

} // namespace
