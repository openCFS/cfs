// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepStokesFluid.hh"

#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepStokesFluid::SolveStepStokesFluid(StdPDE& apde) : StdSolveStep(apde) 
  {
  }

  SolveStepStokesFluid::~SolveStepStokesFluid() {
  }

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepStokesFluid::StepStaticNonLin()
  {

    bool performOneMoreStep;
    UInt iterationCounter=0;
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
    Vector<Double>  actSol = solHelp.GetAlgSysVector();

    Vector<Double> newSol;
    Vector<Double> solIncrement;
    newSol.Resize( numEqns_ );
    solIncrement.Resize( numEqns_ );

    PDE_.SetBCs();

    // store linear part of RHS
    Double loadFactor = 1.0;
    SetLinRHS(loadFactor); 

    //PDE_.AssembleNLRHS();

    std::cout << "In SolveStepStaticNonLin" << std::endl;
    do
      {
        iterationCounter++;

        std::cout << std::endl << "Nonlinear Convective Stokes Fluid: Perform internal loop nr. " 
                  << iterationCounter << std::endl;      

      /*  
        *debug << std::endl
               << "====================================================== "
               << std::endl
               << "Nonlinear Convective Stokes Fluid: Perform internal loop nr. "
               << iterationCounter << std::endl;
       */

        // setup and solve new system (rhs is already set) =====================
        //InitNonLinMatrices();
        assemble_->AssembleMatrices();
      
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver();

        algsys_->Solve();   

        // new solution is only an increment of the full solution =============
        Double *solPtr;
        algsys_->GetSolutionVal( solPtr );
        StoreAlgsysToVec(newSol, solPtr);

        solIncrement = actSol-newSol;

        Double residualL2Norm = 0;

        actSol = newSol;
      
        sol_->SetAlgSysVector(actSol);

        // recalculate RHS with new values to get new residual (f^(k+1))========
        //algsys_->InitRHS(RhsLinVal_.GetPointer());

        //PDE_.AssembleNLRHS();  // inner forces due to nonlin formulation


        // =====================================================================
        // calculation of error norms
        // =====================================================================

        // calculation of residual error =======================================
        //Double residualErr = residualL2Norm / extForcesL2Norm;
        Double residualErr = residualL2Norm;


        // calculate incremental error ========================================
        Double solIncrL2Norm = solIncrement.NormL2();
        Double actSolL2Norm  = actSol.NormL2();
      
        Double incrementalErr;
        Double solAuxU, incAuxU;

        solAuxU=0.0;
        incAuxU=0.0;
        
        for (UInt i=0; i<actSol.GetSize(); i+=4)
          {
            solAuxU+=(actSol[i]*actSol[i]);
            incAuxU+=(solIncrement[i]*solIncrement[i]);
          }
        solAuxU = sqrt(solAuxU);  
        incAuxU = sqrt(incAuxU);  

        Double relAuxU = incAuxU/solAuxU;  

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


//        WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
//                        solIncrL2Norm, actSolL2Norm, incrementalErr);

      
//        Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, 
//                              etaLineSearch);


        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep = 
          (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
      
//        if (!(performOneMoreStep && iterationCounter < nonLinMaxIter_))
          mycout << "solAuxU " << solAuxU << myendl 
                 << "incAuxU " << incAuxU << myendl
                 << "relAuxU " << relAuxU << myendl;

      }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
  }


  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  //  void SolveStepStokesFluid::StepTransNonLin() {
  //  }


  void SolveStepStokesFluid::AddNonLinRHS() {
  }


} // end of namespace

