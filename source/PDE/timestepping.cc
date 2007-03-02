#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"
#include "olas.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(BaseSystem * algebraicsystem )
  {
    ENTER_FCN( "TimeStepping::TimeStepping" );

    algsys_  = algebraicsystem;
    rhsSize_ = 0 ;
    isDeriv1Set_ = false;
    isDeriv2Set_ = false;
  }

  const std::map<FEMatrixType,Double>  &
  TimeStepping::GetEffSysMatFactors( ) const {
    return matrix_factors_;
  }
                                              

  TimeStepping::~TimeStepping()
  {
    ENTER_FCN( "TimeStepping::~TimeStepping" );
  }

  bool TimeStepping::FeMatrixPresent( FEMatrixType type) {
    std::set<FEMatrixType> matTypes;    
    algsys_->GetFEMatrixTypes(matTypes);
    
    if ( matTypes.find(type) != matTypes.end() ) {
      return true;
    } else {
      return false;
    }
  }
     


} // end of namespace
