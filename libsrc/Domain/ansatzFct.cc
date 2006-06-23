#include "ansatzFct.hh"
#include "General/environment.hh"

namespace CoupledField {


  AnsatzFct::AnsatzFct() {
    ENTER_FCN( "AnsatzFct::AnsatzFct" );

  }

  ConstFct::ConstFct() {
    ENTER_FCN( "ConstFct::ConstFct" );
    type_ = CONST;
  }

  LagrangeFct::LagrangeFct() {
    ENTER_FCN( "LagrangeFct::LagrangeFct" );
    type_ = LAGRANGE;
  }

  LegendreFct::LegendreFct() {
    ENTER_FCN( "LegendreFct::LegendreFct" );
  }


}
