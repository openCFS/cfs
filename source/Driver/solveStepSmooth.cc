#include <fstream>
#include <iostream>
#include <string>

#include "solveStepSmooth.hh"
#include "assemble.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepSmooth::SolveStepSmooth(StdPDE& apde) : StdSolveStep(apde) {
  }
  
  
  SolveStepSmooth::~SolveStepSmooth() {
  }

//   void SolveStepSmooth::PreStepTrans() {
//     PreStepStatic();
//   }


//   void SolveStepSmooth::SolveStepTrans() {

//     ENTER_FCN( "SolveStepSmooth::SolveStepTrans" );

//     if ( nonLin_ ) {
//       StepStaticNonLin();
//     }
//     else {
//       StepStaticLin();
//     }
//   }
  
//   void SolveStepSmooth::PostStepTrans(){
//     PostStepStatic();
//   }

  void SolveStepSmooth::StepTransNonLin() {

	UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
    bool isIterCoupled    = PDE_.IsIterCoupled();

    std::string pdeNameLong(pdename_);
    pdeNameLong += "-PDE: ";

    Vector<Double> newSol;
    newSol.Resize( numEqns_ );

    UInt iterationCounter=0;
    // TODO: Check if this is still needed
    // Double iterDeltaT=0.1;
    // Double iterDeltaT2=0.01;

    bool performOneMoreStep;

    if ( isIterCoupled == false || iterCoupledCounter == 0 ) {        
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      TS_alg_->Predictor(solHelp);
    }

    PDE_.InitStabParams();

    do{
      iterationCounter++;
      

      assemble_->CalcMinMaxStrain();

      assemble_->AssembleMatrices();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      // set boundary conditions
      PDE_.SetBCs( );
      algsys_->BuildInDirichlet();

      algsys_->SetupPrecond( );
      algsys_->SetupSolver( );
      algsys_->Solve();
      
      algsys_->GetSolutionVal( newSol );
      sol_->SetAlgSysVector(newSol);

      Vector<Double> solHelp;
      algsys_->GetSolutionVal(solHelp);
      PDE_.SaveSolution(solHelp.GetPointer(),solHelp.GetSize());

      //if(iterationCounter < nonLinMaxIter_ ){
      //  PDE_.PostProcess();
      //  PDE_.WriteResultsInFile(iterationCounter, actTime_-((nonLinMaxIter_-iterationCounter)*iterDeltaT), 0, 0);
      //}

      Vector<Double> actRHS;
      Double RhsL2Norm;
      algsys_->GetRHSVal( actRHS );
      RhsL2Norm = actRHS.NormL2();
      if(RhsL2Norm<1e-9)
        performOneMoreStep=false;
      else
        performOneMoreStep=true;

      if ( nonLinLogging_ == true ) {
        *(Info->GetInfoStreamPointer()) << std::endl << pdeNameLong << "NONLINEAR ITERATION "
                                        << iterationCounter << " of "
                                        << pdeNameLong << " =============================== "
                                        << " RhsL2Norm=" << RhsL2Norm
                                        << std::endl;
      }
              
    }while( performOneMoreStep && iterationCounter < nonLinMaxIter_ );  


    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    TS_alg_->Corrector(solHelp);

    if ( isIterCoupled ) {
      iterCoupledCounter++;
    }

  }

} // end of namespace

