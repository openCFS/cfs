#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "bdf2.hh"

namespace CoupledField
{
  DECLARE_LOG(bdf2)
  DEFINE_LOG(bdf2, "bdf2")

  Bdf2 :: Bdf2( BaseSystem * algebraicsystem)
    :TimeStepping( algebraicsystem )
  {

    // Commented out the warning, since defaults are not bad at all and the 
    // average student user gets not disturbed by any warnings
    //check if integration parameters are defined in conf-file
    //Info->Warning( "Bdf2: Using defaults for gamma!" );
    firstTime_=true;

  }

  Bdf2 :: ~Bdf2()
  {
  }

  void Bdf2::Init( Double dt, UInt rhsSize ) {
    dt_ = dt;
    rhsSize_ = rhsSize;
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0; //dt_*2.0/3.0
    matrix_factors_[MASS] = 3.0/(2.0*dt_); // 1.0

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 0.0;       

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
    if( !isDeriv2Set_ ) {
      solderiv2_.Resize(rhsSize_);
      solderiv2_.Init();
    }

//     solpred_.Resize(rhsSize_); 
//     solpred_.Init();
  }


  void Bdf2::Predictor(Vector<Double>& solold)
  {
    //std::cout << "solderiv1_\n" << solderiv1_ << std::endl;
    //if(firstTime_){
      //use old values from restart file; currently hard coded
      //for bdf2 solderiv1_ and solderiv2_ doesnot contain derivatives!!!!
      //sol_tn_2_ = solderiv2_; //see corrector
      //sol_tn_1_ = solold;
      //sol_tn_2_ = solold;
      //sol_tn_1_ = solold;

      //firstTime_=false;
      //std::cerr << "*\n";
    //}
    //else{
    sol_tn_2_ = sol_tn_1_;
    sol_tn_1_ = solold;
      //std::cerr << ".\n";
    //}
    //std::cout << "\n Bdf2 Pred sol_tn_1_: \n" << sol_tn_1_ << std::endl;
    //std::cout << "\n Bdf2 Pred sol_tn_2_: \n" << sol_tn_2_ << std::endl;
    //LOG_DBG(bdf2) << "\n Pred sol_tn_1_: \n" << sol_tn_1_ << std::endl;
    //LOG_DBG(bdf2) << "\n Pred sol_tn_2_: \n" << sol_tn_2_ << std::endl;

  }


  void Bdf2::UpdateRHS()
  {
    Vector<Double> coeffMass;

    // mass part
    coeffMass = (sol_tn_1_*a0_) - (sol_tn_2_*a1_);
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());
  }


  void Bdf2::Corrector(Vector<Double>& solnew)
  {
    //in solderiv1_ and solderiv2_ for bdf2 not the derivatives are stored.
    // in order to realize restart functionality the last solutions [t_(n-1) t_(n-2)} are stored to retrace a simulation
    //solderiv1_=solnew;
    //solderiv2_=sol_tn_1_;

    //std::cout << "\n Bdf2 Cor solnew: \n" << solnew << std::endl;
    //std::cout << "\n Bdf2 Cor solderiv2_: \n" << solderiv2_ << std::endl;
    //LOG_DBG(bdf2) << "\n Cor solnew: \n" << solnew << std::endl;
    //LOG_DBG(bdf2) << "\n Cor solderiv2_: \n" << solderiv2_ << std::endl;
  }

  void Bdf2 :: CalcParameters(Double dt)
  {
    //for u_{n-1}
    a0_ = 2.0/dt_;
    //a0_ = 4.0/3.0;

    //for u_{n-2}
    a1_ = 1.0/(2.0*dt_);
    //a1_ = 1.0/3.0;
  }


} // end of namespace
