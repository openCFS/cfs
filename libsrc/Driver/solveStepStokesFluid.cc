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

