#include <fstream>
#include <iostream>
#include <string>

#include "solveStepSmooth.hh"
#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"


namespace CoupledField {

  SolveStepSmooth::SolveStepSmooth(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepSmooth::SolveStepSmooth" );
  }

  
  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepSmooth::StepStaticNonLin( const Boolean reset )
  {
    ENTER_FCN( "SolveStepSmooth::StepStaticNonLin" );

    UInt job = 1;
    Double * ptsol;

    assemble_->AssembleMatrices();
    assemble_->AssembleSrcRHS();
  
    PDE_.SetBCs( 0.0 );
    
    algsys_->ConstructEffectiveMatrix(matrix_factor_);
    algsys_->BuildInDirichlet();

    if (job == 1) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }

    algsys_->Solve();

    algsys_->GetSolutionVal( ptsol );

    // save solution
    sol_->CopyFromAlgSysDataPointer(ptsol);
  }


  //   Default Destructor
  // **********************
  SolveStepSmooth::~SolveStepSmooth() {

    ENTER_FCN( "SolveStepSmooth::~SolveStepSmooth" );
 
  }

} // end of namespace

