#include <fstream>
#include <iostream>
#include <string>

#include "solveStepSmooth.hh"


namespace CoupledField {

  SolveStepSmooth::SolveStepSmooth(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepSmooth::SolveStepSmooth" );
  }

  
  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================
  void SolveStepSmooth::PreStepStatic(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepSmooth::PreStepStatic" );

    algsys_->InitRHS();
    algsys_->InitSol();
    assemble_->InitMatrices();

    assemble_->SetReassemble();

  }

  void SolveStepSmooth::StepStaticNonLin(const Integer kstep, const Double aTime,
				   const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepSmooth::StepStaticNonLin" );

    Integer job = 1;
    Double * ptsol;

    assemble_->AssembleMatrices(level);
    assemble_->AssembleSrcRHS(level);
  
    SetBCs(level,0);

    algsys_->BuildInDirichlet();
    if (job == 1) {
      algsys_->SetupPrecond(job);
      algsys_->SetupSolver(job);
    }

    algsys_->Solve();

    ptsol = algsys_->GetSolutionVal();

    // save solution
    sol_->CopyFromAlgSysDataPointer(ptsol);
    //  sol_->SetAlgSysDataPointer(ptsol);
  }


  void SolveStepSmooth:: PostStepStatic(const Integer kstep, const Double asteptime,
				  const Integer level)
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

