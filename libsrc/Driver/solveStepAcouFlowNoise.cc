#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcouFlowNoise.hh"
#include "Forms/forms_header.hh"
#include "PDE/acouFlowNoise.hh"

namespace CoupledField {

  SolveStepAcouFlowNoise::SolveStepAcouFlowNoise(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepAcouFlowNoise::SolveStepAcouFlowNoise" );
  }

  
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepAcouFlowNoise::SolveStepTrans( const Boolean reset ) {
    ENTER_FCN( "SolveStepAcouFlowNoise::SolveStepTrans" );
    
    //    Boolean Recalc=FALSE;

    //if (laststepcalc_==kstep && kstep!=0) Recalc=TRUE;
    //else laststepcalc_= kstep;

    Double * ptsol;
    UInt job = 3;

    //perform predictor step
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    if ( actStep_ == 1 )
      {
        job = 1;
        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);

        algsys_->InitRHS();
        assemble_->AssembleSrcRHS( actTime_ );

        PDE_.ComputeRHS( actTime_ );
        TS_alg_->UpdateRHS();
      }
    else if (reset)
      {
        job    = 1;

        algsys_->InitMatrix(SYSTEM);
        algsys_->ConstructEffectiveMatrix(matrix_factor_);

        algsys_->InitRHS();
        assemble_->AssembleSrcRHS( actTime_ );
        PDE_.ComputeRHS( actTime_ );
        TS_alg_->UpdateRHS();
      }
    else
      {
        job    = 3;
        algsys_->InitRHS();
        assemble_->AssembleSrcRHS( actTime_ );
        PDE_.ComputeRHS( actTime_ );
        TS_alg_->UpdateRHS();
      };

    SetBCs( actTime_ );

    if ( job == 1 ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }

    algsys_->Solve();

    // Save solution
    algsys_->GetSolutionVal( ptsol );
    sol_->CopyFromAlgSysDataPointer(ptsol);
  
    //perform corrector step 
    TS_alg_->Corrector(solhelp->GetAlgSysVector());
  }


  //   Default Destructor
  // **********************
  SolveStepAcouFlowNoise::~SolveStepAcouFlowNoise() {

    ENTER_FCN( "SolveStepAcouFlowNoise::~SolveStepAcouFlowNoise" );
 
  }

} // end of namespace

