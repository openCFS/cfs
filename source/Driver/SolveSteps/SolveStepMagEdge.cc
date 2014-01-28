#include <fstream>
#include <iostream>
#include <string>

#include "SolveStepMagEdge.hh"
#include "PDE/MagEdgePDE.hh"
#include "Driver/Assemble.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  SolveStepMagEdge::SolveStepMagEdge(StdPDE& apde) : 
      StdSolveStep(apde)
  {
    magEdgePDE_ = dynamic_cast<MagEdgePDE*>(&apde);
  }
  
  
  SolveStepMagEdge::~SolveStepMagEdge() {
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
    //StdSolveStep::StepStaticNonLin(analysis_id);

    bool performOneMoreStep;

    SBM_Vector solInc(BaseMatrix::DOUBLE);

    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;

    shared_ptr<SolStrategy> solStrat_ = algsys_->GetSolStrategy();

    // =================================
    //  Outer loop: Multilevel strategy
    // =================================
    UInt numLevels = solStrat_->GetNumSolSteps();
    for( UInt iLevel = 0; iLevel < numLevels; ++iLevel ) {

      // create new timer object and put it to related info element
      shared_ptr<Timer> timer(new Timer());
      PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
      iter->GetByVal("solStep","value",iLevel+1,ParamNode::INSERT)
          ->Get("timer")->SetValue(timer);
      timer->Start();

      // update the current solution step in a multilevel approach and
      // inform PDEs (containing the FeSpaces), as well as the AlgebraicSystem
      solStrat_->SetActSolStep(iLevel + 1);
      ReadNonLinData();
      PDE_.UpdateToSolStrategy();
      algsys_->UpdateToSolStrategy();

      // Inititalize RHS load vector
      algsys_->InitRHS();

      // set the boundary conditions
      PDE_.SetBCs();

      //perform the load-steps
      Double loadFactor = 1.0;
      PDE_.GetInfoNode()->Get("PDE")->Get(pdename_)->
          Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // assemble nonlinear parts to RHS
      assemble_->AssembleNonLinRHS();

      // set iteration counter
      UInt iterationCounter=0;

      // =================================
      //  Inner nonlinear loop
      // =================================
      do {
        iterationCounter++;
        // RHS is already set up!!

//        PtrParamNode child_id =
//            BaseDriver::CreateAnalysisIdChild(analysis_id, "nonLin", iterationCounter);

        // setup and solve new system (rhs is already set) =====================
        assemble_->AssembleMatrices();
        bool isNewtonPart = true;
        assemble_->AssembleMatrices(isNewtonPart);

        algsys_->ConstructEffectiveMatrix( NO_FCT_ID,
                                           matrix_factor_[NO_FCT_ID] );

        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond(analysis_id);
        algsys_->SetupSolver(analysis_id);

        bool setIDBC = false;
        if ( iterationCounter == 1 )
          setIDBC = true;

        algsys_->Solve(analysis_id, setIDBC);

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solInc, setIDBC );


        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;
        if ( lineSearch_ == "none" ) {
          actSol.Add(1.0, solInc);
        }
        else {
          // true is for transient simulation
          residualL2Norm = LineSearchMag(solInc, actSol, etaLineSearch);
        }

        // store the new solution
        solVec_ = actSol;

        if ( lineSearch_ == "none" ) {
          // recalculate RHS with new values to get new residual (f^(k+1))========
          algsys_->InitRHS(RhsLinVal_);
          assemble_->AssembleNonLinRHS();
          //Set special RHS Values
          PDE_.SetRhsValues();

          //get RHS vector
          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error =======================================
          residualL2Norm = actRHS.NormL2(); // L2Norm of  ( f_i^(k+1) - f_a )
        } else {
          algsys_->InitRHS(RhsLinVal_ );
          assemble_->AssembleNonLinRHS();
          //Set special RHS Values
          PDE_.SetRhsValues();
        }

        // calculation of residual error =======================================
        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr = residualL2Norm / RhsLinL2Norm;
        else
          residualErr = residualL2Norm;

        // calculate incremental error ========================================
        Double incrementalErr;
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm  = actSol.NormL2();

        if ( actSolL2Norm )
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else {
          incrementalErr = solIncrL2Norm;
          WARN("Zero solution vector!! ");
        }

        //std::cout << "Norm:\n" << "Residual: " << residualErr << "  Incr.Error: " <<  incrementalErr << std::endl << std::endl;
        // output of norms and data
        if ( nonLinLogging_ == true ) {
          WriteNonLinIterToInfoXML(pdename_, iLevel+1,iterationCounter, residualErr,
                                   incrementalErr, etaLineSearch);
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
              << residualErr << "\t"
              << incrementalErr << "\t"
              << etaLineSearch << std::endl;
        }

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
            (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

      // stop timer
      timer->Stop();
    } // loop over levels
    
  }
  





} // end of namespace

