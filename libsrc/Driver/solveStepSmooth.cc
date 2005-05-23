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
  void SolveStepSmooth::PreStepStatic(const Integer kstep, const Double asteptime,
                                      const Boolean reset)
  {
    ENTER_FCN( "SolveStepSmooth::PreStepStatic" );

    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();
    assemble_->SetReassemble();

  }

  void SolveStepSmooth::StepStaticNonLin(const Integer kstep, const Double aTime,
                                         const Boolean reset)
  {
    ENTER_FCN( "SolveStepSmooth::StepStaticNonLin" );

    Integer job = 1;
    Double * ptsol;

    assemble_->AssembleMatrices();
    assemble_->AssembleSrcRHS();
  
    SetBCs(0);

    algsys_->BuildInDirichlet();
    if (job == 1) {
      algsys_->SetupPrecond(job);
      algsys_->SetupSolver(job);
    }

    algsys_->Solve();

    algsys_->GetSolutionVal( ptsol );

    // save solution
    sol_->CopyFromAlgSysDataPointer(ptsol);
  }


  void SolveStepSmooth:: PostStepStatic(const Integer kstep, const Double asteptime)
  {
    ENTER_FCN( "SolveStepSmooth::PostStepStatic" );
  
    iterCoupledCounter_++;
  }

  //   Default Destructor
  // **********************
  SolveStepSmooth::~SolveStepSmooth() {

    ENTER_FCN( "SolveStepSmooth::~SolveStepSmooth" );
 
  }

} // end of namespace

