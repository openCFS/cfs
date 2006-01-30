#include <fstream>
#include <iostream>
#include <string>

#include "solveStepStokesFluid.hh"
#include "Forms/forms_header.hh"

#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepStokesFluid::SolveStepStokesFluid(StdPDE& apde) : StdSolveStep(apde) 
  {
    ENTER_FCN( "SolveStepStokesFluid::SolveStepStokesFluid" );
  }

  SolveStepStokesFluid::~SolveStepStokesFluid() {
    ENTER_FCN( "SolveStepStokesFluid::~SolveStepStokesFluid" );
  }

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepStokesFluid::PreStepStatic( const Boolean reset ) {

    ENTER_FCN( "SolveStepStokesFluid::PreStepStatic" );

    // Set right-hand side to zero
    //
    // Note: Though this is PreStepStatic, this is most important
    //       for transient analysis ;-)
    algsys_->InitRHS();

    if (isIterCoupled_)     
      algsys_->InitSol();

    if (geoUpdate_) {
      algsys_->InitSol();
      algsys_->InitMatrix();
      assemble_->SetReassemble();   
    }
  }

  void SolveStepStokesFluid::StepStaticNonLin( const Boolean reset )
  {
    ENTER_FCN( "SolveStepStokesFluid::SolveStepStaticNonLin" );

    Boolean performOneMoreStep;
    UInt iterationCounter=0;
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
    Vector<Double>  actSol = solHelp.GetAlgSysVector();

     Vector<Double> solIncrement;
    solIncrement.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());

    PDE_.SetBCs(0);

    // store linear part of RHS
    Double loadFactor = 1.0;
    Double extForcesL2Norm = SetLinRHS(loadFactor); 

    PDE_.AssembleNLRHS();

   std::cout << "In SolveStepStaticNonLin" << std::endl;
    do
      {
    iterationCounter++;

    std::cout << std::endl << "Nonlinear Convective Stokes Fluid: Perform internal loop nr. " 
          << iterationCounter << std::endl;      

#ifdef DEBUG
    *debug << std::endl
           << "====================================================== "
           << std::endl
           << "Nonlinear Convective Stokes Fluid: Perform internal loop nr. "
           << iterationCounter << std::endl;
#endif


    // setup and solve new system (rhs is already set) =====================
    PDE_.InitNonLinMatrices();
    PDE_.AssembleMatrices();
      
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
    algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver();

    algsys_->Solve();   

    // new solution is only an increment of the full solution =============
    Double *solPtr;
    algsys_->GetSolutionVal( solPtr );
    StoreAlgsysToVec(solIncrement, solPtr);

    Double residualL2Norm = 0;
    Double etaLineSearch=0;

    if ( lineSearch_ == "none" )
      actSol = solIncrement;
    else
      // TRUE is for transient simulation
      residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch);
      
    sol_->SetAlgSysVector(actSol);

    // recalculate RHS with new values to get new residual (f^(k+1))========
    algsys_->InitRHS(RhsLinVal_.GetPointer());

    PDE_.AssembleNLRHS();  // inner forces due to nonlin formulation


    // =====================================================================
    // calculation of error norms
    // =====================================================================
    if ( lineSearch_ == "none" )
      {
        Vector<Double> actRHS;
        algsys_->GetRHSVal(solPtr);
        StoreAlgsysToVec(actRHS, solPtr );       
          
        // calculation of residual error =======================================
        residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
      }


    // calculation of residual error =======================================
    //Double residualErr = residualL2Norm / extForcesL2Norm;
    Double residualErr = residualL2Norm;


    // calculate incremental error ========================================
    Double solIncrL2Norm = solIncrement.NormL2();
    Double actSolL2Norm  = actSol.NormL2();
      
    Double incrementalErr;
      
    if (actSolL2Norm)
      incrementalErr = solIncrL2Norm / actSolL2Norm;
    else
      {
        incrementalErr = solIncrL2Norm;
        Warning("Zero solution vector!! ", __FILE__,__LINE__);      
      }
      


    // =====================================================================
    // output of norms and data
    // =====================================================================


    WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
            solIncrL2Norm, actSolL2Norm, incrementalErr);

      
    Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, 
                  etaLineSearch);


    // boolean variable, holds condition if another iteration step is necessary
    performOneMoreStep = 
      (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
      
    if (!(performOneMoreStep && iterationCounter < nonLinMaxIter_))
      mycout << "incrementalErr " << incrementalErr << myendl 
         << "incStopCrit_ " << incStopCrit_ << myendl
         << "residualErr " << residualErr  << myendl
         << "residualStopCrit_ " << residualStopCrit_ << myendl;

      }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  }

  void SolveStepStokesFluid::PostStepStatic()
  {
    ENTER_FCN( "SolveStepStokesFluid::PostStepStatic" );

    if (isIterCoupled_)
      (*iterCoupledCounter_)++;
  }


 
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepStokesFluid::StepTransNonLin( const Boolean reset ) {
    ENTER_FCN( "SolveStepStokesFluid::StepTransNonLin" );
  }


  void SolveStepStokesFluid::AddNonLinRHS() {
    ENTER_FCN( "SolveStepStokesFluid::AddNonLinRHS" );

  }


} // end of namespace

