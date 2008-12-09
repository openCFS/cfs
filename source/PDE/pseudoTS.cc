#include <fstream>
#include <iostream>
#include <string>

//#include "DataInOut/WriteInfo.hh"
#include "pseudoTS.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{
  // declare logging stream
//  DECLARE_LOG(pseudoTS)
//  DEFINE_LOG(pseudoTS, "pseudoTS")

  PseudoTS::PseudoTS(BaseSystem * algebraicsystem )
    :TimeStepping(algebraicsystem )
  {
  }

  PseudoTS::~PseudoTS()
  {
  }

  void PseudoTS::Init( Double dt, UInt rhsSize ) {

    rhsSize_ = rhsSize;

    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
        
    matrix_factors_[STIFFNESS]  = 1.0;
    matrix_factors_[MASS]       = 0.0;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING]    = 0.0;

    if( !isSolTN1Set_){
      sol_tn_1_.Resize(rhsSize_);
      sol_tn_1_.Init();
    }

    sol_tn_2_.Resize(rhsSize_);
    sol_tn_2_.Init();

    //get the memory
    if( !isDeriv1Set_ ) {
      solderiv1_.Resize(rhsSize_);
      solderiv1_.Init();
    }
    //get the memory
    if( !isDeriv2Set_ ) {
      solderiv2_.Resize(rhsSize_);
      solderiv2_.Init();
    }
  }

  void PseudoTS::Predictor(Vector<Double>& solold)
  {
    sol_tn_2_ = sol_tn_1_;
    sol_tn_1_ = solold;
  }

  void PseudoTS::Corrector(Vector<Double>& solnew)
  {
    solderiv1_ = (solnew*3.0 - sol_tn_1_*4.0 + sol_tn_2_ )/(dt_*2.0);
  }

} // end of namespace
