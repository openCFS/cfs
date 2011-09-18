// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMagEdge.hh"
#include "PDE/magEdgePDE.hh"
#include "assemble.hh"
#include "Forms/linearForm.hh"
#include "OLAS/algsys/algebraicSys.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  SolveStepMagEdge::SolveStepMagEdge(StdPDE& apde) : 
      StdSolveStep(apde)
  {
   strategy_ = STRAT_STANDARD;
    magEdgePDE_ = dynamic_cast<MagEdgePDE*>(&apde);
  }
  
  
  SolveStepMagEdge::~SolveStepMagEdge() {
  }
  
  void SolveStepMagEdge::SetSolStrategy( SolStrategyType type ) {
    strategy_ = type;
  }

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMagEdge::PreStepStatic() {
    StdSolveStep::PreStepStatic();
  }


  void SolveStepMagEdge::SolveStepStatic(PtrParamNode analysis_id) {



    // Check, if we have a two-level solving strategy
//    if( strategy_ == STRAT_STANDARD ) {
      std::ofstream clockFile_;
      clockFile_.open("time.txt");
      Timer clock;
      clock.Start();
           
      if (nonLin_) 
        StepStaticNonLin(analysis_id);
      else 
        StepStaticLin(analysis_id);
      

      Double wall = clock.GetWallTime();
      Double user = clock.GetCPUTime();
      clock.Stop();
      clockFile_ << "1 " << wall << "\t" << user << std::endl;
      clockFile_.close();
      
//    } else if( strategy_ == STRAT_TWO_LEVEL ){
//
//      // create hard cocded timer file
//      std::ofstream clockFile_;
//      clockFile_.open("time.txt");
//      Timer clock;
//      clock.Start();
//      
//      std::cerr << " *********************************\n";
//      std::cerr << "        TWO LEVEL APPROACH \n";
//      std::cerr << " *********************************\n";
//      // Rely on fact, that PDE recognizes itself, that the 
//      // first computation is performed on the "small" system
//      // Here we have to take care of the following issues:
//      // - make sure, that a direct solver is chosen (not the definition for
//      //   the second step, where we use an iterative solver)
//      // - The order of the FeSpace has to be hard-coded to 1
//
//
//      
//      
//      std::cerr << " *********************************\n";
//      std::cerr << "       1) STARTING 1ST STEP\n";
//      std::cerr << " *********************************\n";
//      // Thus we solve the first "coarse" level like for the standard case
//
//      if (nonLin_) 
//        StepStaticNonLin(analysis_id);
//      else 
//        StepStaticLin(analysis_id);
//      
//      
//      Double wall = clock.GetWallTime();
//      Double user = clock.GetCPUTime();
//      clock.ResetStart();
//      clockFile_ << "1 " << wall << "\t" << user << std::endl;
//      
//      std::cerr << " *********************************\n";
//      std::cerr << "       2) STARTING 2ND STEP\n";
//      std::cerr << " *********************************\n";
//
//      // Get nodestoresol and store solution vector
//      Vector<Double> oldSolVec;
//
//      dynamic_cast<NodeStoreSol<Double>* >(sol_)->GetAlgSysVector(oldSolVec);
//
//      // store also old RHS vector
//      Vector<Double> oldRhsLinVal = RhsLinVal_;
//
//
//      // Re-Init the PDE -> Have a look at Fabians implementation, 
//      // to find out what exactly has to be deleted
//      // => (Bi)LinearForms, Assemble-class, Algebraic System, FeSpace
//      //    FeFunction, NodeStoreSol etc.)
//      PDE_.SetSolutionStep(2);
//
//      //        pdename_      = PDE_.GetName();
//      numPDENodes_  = PDE_.GetNumPdeEquations();
//      numPDEElems_  = PDE_.getPDE_numElems();
//      //        isaxi_        = PDE_.GetIsaxi();
//      //        subdoms_      = PDE_.getPDE_subdoms();
//      //        materialData_ = PDE_.getPDEMaterialData();
//      //        ptgrid_       = PDE_.getPDE_grid();
//      algsys_       = PDE_.getPDE_algsys();
//      //        sol_          = PDE_.getPDESolution();
//      //        assemble_     = PDE_.getPDE_assemble();
//      //        TS_alg_       = PDE_.getTimeStepping();
//
//
//      // ==== Now we have the new setup, i.e. new integrators, the new 
//      // Fespace etc. ===
//
//      // Obtain again the nodestoresol object and copy solution of previous
//      // simulation run back ( => "projection" of coarse to fine space)
//
//      NodeStoreSol<Double> & newSol = dynamic_cast<NodeStoreSol<Double>& >(*sol_); 
//      Vector<Double> * newSolVec;
//      newSol.GetAlgSysVectorPointer(newSolVec);
//      
//      RhsLinVal_.Resize(newSolVec->GetSize());
//      // copy initial values and also rhs vector
//      for( UInt i = 0; i < oldSolVec.GetSize(); i++ ) {
//        (*newSolVec)[i] = oldSolVec[i];
//        RhsLinVal_[i] = oldRhsLinVal[i];
//      }
//      //std::cerr << "newSolVec is\n" << *newSolVec << std::endl;
//
//      // Initialize OLAS with the correct solution
//      algsys_->InitSol(*newSolVec);
//      // set also the correct RHS vector
//      //algsys_->InitRHS(RhsLinVal_);
//      // Careful: We also have to set up the RHS vector correctly, i.e. 
//      // the nonlinear terms have also te be present
//      //assemble_->AssembleNonLinRHS();
//      algsys_->InitRHS();
//      
//
//      // Perform again the StepStaticLin
//      // ATTENTION: Of course, we want now the better initial solution so 
//      // prevent re-setting of the solution vector in the StepStatic(Lin)
//      if (nonLin_) 
//        StepStaticNonLin(analysis_id);
//      else 
//        StepStaticLin(analysis_id);
//
//      wall = clock.GetWallTime();
//      user = clock.GetCPUTime();
//      clock.Stop();
//      clockFile_ << "2 " << wall << "\t" << user << std::endl;
//      clockFile_.close();
//    } else {
//      EXCEPTION("Solver strategy '" << SolStrategyEnum.ToString(strategy_)
//                << "' not yet implemented!");
//    }


  }


  void SolveStepMagEdge::StepStaticLin( PtrParamNode analysis_id ) {
    StdSolveStep::StepStaticLin(analysis_id);  
  }

  void SolveStepMagEdge::StepStaticNonLin( PtrParamNode analysis_id ) {
    
    // Note: currently hard-coded to section from StdSolveStep
    StdSolveStep::StepStaticNonLin(analysis_id);
    
  }
  
} // end of namespace

