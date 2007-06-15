// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

#include "stdSolveStep.hh"
#include "assemble.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/singleDriver.hh"


namespace CoupledField {

  StdSolveStep::StdSolveStep(StdPDE & apde) 
    :BaseSolveStep(),
     PDE_(apde)
  {

    ENTER_FCN( "StdSolveStep::StdSolveStep" );
 

    pdename_      = PDE_.GetName();
    numPDENodes_  = PDE_.getPDE_numPDENodes();
    numPDEElems_  = PDE_.getPDE_numElems();
    isaxi_        = PDE_.GetIsaxi();
    subdoms_      = PDE_.getPDE_subdoms();
    materialData_ = PDE_.getPDEMaterialData();
    ptgrid_       = PDE_.getPDE_grid();
    algsys_       = PDE_.getPDE_algsys();
    sol_          = PDE_.getPDESolution();
    assemble_     = PDE_.getPDE_assemble();
    TS_alg_       = PDE_.getTimeStepping();

    if( TS_alg_ != NULL ) {
      matrix_factor_ = TS_alg_->GetEffSysMatFactors();
    }
    
    eqnMap_       = PDE_.GetEqnMap();
    results_      = PDE_.GetResultInfos();
    numEqns_      = PDE_.GetSolutionVector()->GetSize();

    // nonlinear parameters
    incStopCrit_ = 1e-2; 
    residualStopCrit_ = 1e-3;
    nonLin_           = PDE_.IsNonLin();
    nonLinMaterial_   = PDE_.IsNonLinMaterial();
    isHyst_           = PDE_.IsHysteresis();
    regionNonLinType_ = PDE_.GetNonLinRegionTypes();

    // for direct coupled PDEs
    pdeId1_   = NO_PDE_ID;
    pdeId2_   = NO_PDE_ID;

    startStep_ = 1;

    // In the end, read nonlinear data from xml-file
    ReadNonLinData();
  }

  
  //! Destructor
  StdSolveStep::~StdSolveStep() {
    ENTER_FCN( "StdSolveStep::~StdSolveStep" );

  }


  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================

  void StdSolveStep::PreStepStatic( ) {
    ENTER_FCN( "StdSolveStep::PreStepStatic" );
    
    // init RHS at this place, because e.g. forces of other PDEs are added 
    // to RHS afterwards
    algsys_->InitRHS();        
    
  }


  void StdSolveStep::PostStepStatic() {
    ENTER_FCN( "StdSolveStep::PostStepStatic" );

    // increment coupling counter
    if ( PDE_.IsIterCoupled() ) {
      UInt& counter = PDE_.GetIterCoupledCounter();
      counter++;
    }
  }

  
  void StdSolveStep::SolveStepStatic() {
    ENTER_FCN( "StdSolveStep::SolveStepStatic" );
    
    if (nonLin_) {
      StepStaticNonLin();
    }
    else {
      StepStaticLin();
    }
  }


  void StdSolveStep::StepStaticLin() {
    ENTER_FCN( "StdSolveStep::StepStaticLin" );
    
    assemble_->AssembleMatrices();
    
    // The RHS-sources and boundary conditions 
    // have to be reassembled each time
    assemble_->AssembleLinRHS( actTime_ );
    PDE_.SetBCs( actTime_ );

    UInt length = 0;
    
    // store rhs vector back to PDE
    Double * rhsPt;
    length = algsys_->GetRHSVal(rhsPt);
    PDE_.SaveRHS( rhsPt, length );
    
    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }
    
    // Incorporate Boundary conitions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    SETPROFILE("Before SetupPrecond");
    
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      SETPROFILE("After SetupPrecond / Before SetupSolver");
      
      algsys_->SetupSolver( );
      SETPROFILE("After SetupSolver / Before Solve");
    }
    
    // Solve problem
    algsys_->Solve();
    SETPROFILE("After Solve");
    
    // Get the solution and store it
    Double * ptSol;
    UInt size = algsys_->GetSolutionVal(ptSol);

    PDE_.SaveSolution(ptSol,size);

  }


  void StdSolveStep::StepStaticNonLin()
  {
    ENTER_FCN( "StdSolveStep::SolveStepStaticNonLin" );

    bool performOneMoreStep;
 
    Vector<Double> solInc( numEqns_ );
    Vector<Double> actSol( numEqns_ );

    // get solution from algsys
    sol_->GetAlgSysVector(actSol);

    // set the boundary conditions
    PDE_.SetBCs(0);

    //perform the load-steps
    Double loadFactor = 0.0;

    // currently just for testing!!
    // loop over load factor
    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      Info->PrintF(pdename_, "\n");
      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side 
      Double RhsLinL2Norm = SetLinRHS(loadFactor); 

      // assemble nonlinear pars to RHS
      assemble_->AssembleNonLinRHS( actTime_ );

      // set iteration counter
      UInt iterationCounter=0;

      do
        {
          iterationCounter++;
          // RHS is already set up!!

          // setup and solve new system (rhs is already set) =====================
          //assemble_.InitNonLinMatrices();
          assemble_->AssembleMatrices();
      
          algsys_->ConstructEffectiveMatrix(matrix_factor_);
          algsys_->BuildInDirichlet();
          algsys_->SetupPrecond();
          algsys_->SetupSolver();
          algsys_->Solve();   

          // new solution is only an increment of the full solution =============
          Double *solPtr;
          algsys_->GetSolutionVal( solPtr );
          StoreAlgsysToVec(solInc, solPtr);

          Double residualL2Norm = 0.0;
          Double etaLineSearch  = 1.0;
          if ( lineSearch_ == "none" ) {
            actSol += solInc;
          }
          else {
            // true is for transient simulation
            residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
          }

          // store the new solution
          sol_->SetAlgSysVector(actSol);

          if ( lineSearch_ == "none" ) {
            // recalculate RHS with new values to get new residual (f^(k+1))========
            algsys_->InitRHS(RhsLinVal_.GetPointer());
            assemble_->AssembleNonLinRHS( actTime_ );  

            // compute the norm of the residual
            Vector<Double> actRHS;
            algsys_->GetRHSVal(solPtr);
            StoreAlgsysToVec(actRHS, solPtr );       
          
            // calculation of residual error =======================================
            residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
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
            Warning("Zero solution vector!! ", __FILE__,__LINE__);      
          }
      
          // output of norms and data
          if ( nonLinLogging_ == true ) {
            Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                                  incrementalErr, etaLineSearch);
          }	  

          // boolean variable, holds condition if another iteration step is necessary
          performOneMoreStep = 
            (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
        } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    } // load step loop

  }


  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void StdSolveStep::PreStepTrans() {

    ENTER_FCN( "StdSolveStep::PreStepTrans" );

    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();
    
  }


  void StdSolveStep::SolveStepTrans() {

    ENTER_FCN( "StdSolveStep::SolveStepTrans" );  

    // do a nonlinear material time step
    if ( nonLin_ && nonLinMaterial_ ) {
      StepTransNonLinMaterial();
    }
    else if ( nonLin_ && isHyst_ ) {
       StepTransNonLinHysteresis();
    }
    // do a nonlinear time step
    else if (nonLin_){
      StepTransNonLin();
    }
    // do a linear time step
    else {
      StepTransLin();
    }
  }


  void StdSolveStep::StepTransLin() {

    ENTER_FCN( "StdSolveStep::StepTransLin" );

    //account for RHS
    assemble_->AssembleLinRHS( actTime_ );
    PDE_.ComputeRHS( actTime_ );

    UInt length = 0;
    
    // store rhs vector back to PDE
    Double * rhsPt;
    length = algsys_->GetRHSVal(rhsPt);
    PDE_.SaveRHS( rhsPt, length );

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
    bool isIterCoupled    = PDE_.IsIterCoupled();
  

    // perform predictor step: if we have an iterative coupled 
    // PDE-system, we should perform the predictor state just 
    // in the first iteration
    if ( isIterCoupled == false || iterCoupledCounter == 0 ) {        
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      TS_alg_->Predictor(solHelp);
    }

    assemble_->AssembleMatrices();

    if (assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }


    // update Right Hand Side
    TS_alg_->UpdateRHS();

    // set boundary conditions
    PDE_.SetBCs( actTime_ );
    algsys_->BuildInDirichlet();

    if (assemble_->IsMatrixUpdated() ) {
      SETPROFILE("Before SetupPrecond");
      algsys_->SetupPrecond( );
      SETPROFILE("After SetupPrecond / Before SetupSolver");
      algsys_->SetupSolver( );
      SETPROFILE("After SetupSolver / Before Solve");  
    }

    algsys_->Solve(actStep_);
    SETPROFILE("After Solve");

    Double* ptsol;
    length = algsys_->GetSolutionVal(ptsol);
    PDE_.SaveSolution(ptsol,length);
   
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());


    TS_alg_->Corrector(solHelp);

    if ( isIterCoupled ) {
      iterCoupledCounter++;
    }
   
  }



  void StdSolveStep::StepTransNonLin() {

    ENTER_FCN( "StdSolveStep::StepTransNonLin" );

    Double *solPtr;

    bool performOneMoreStep;

    Vector<Double> solInc( numEqns_ );
    Vector<Double> actSol( numEqns_ );
    actSol.Init();

    UInt numNodes;
    numNodes = sol_->GetNumNodes();
  
    // Cast BaseStoreSol into StoreSol<Double>,
    // since this function is only called
    // in the transient case


    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //get actual solution  
    actSol = solhelp->GetAlgSysVector();

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Error( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin",
             __FILE__, __LINE__ );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( actTime_ );

    // currently just for testing!!
    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      Info->PrintF(pdename_, "\n");
      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side 
      Double RhsLinL2Norm = SetLinRHS(loadFactor); 

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS( actTime_ );  

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(solhelp->GetAlgSysVector());

      // set iteration counter
      UInt iterationCounter=0;


      do {
        iterationCounter++;

        // RHS is already set up!!

        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();
  
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        algsys_->Solve();

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solPtr );
        StoreAlgsysToVec(solInc, solPtr );

        Double residualL2Norm;
        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none" ) {
          actSol += solInc;
        }
        else {
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
        }
      
        //store A_(n+1) in the solution-object sol_
        sol_->SetAlgSysVector(actSol);

        if ( lineSearch_ == "none" ) {
          // calculation of error norms
          // recalculate RHS with new values to get new residual (f^(k+1))=======
          algsys_->InitRHS(RhsLinVal_.GetPointer());

          //Update RHS (mass matrix on right hand side)
          TS_alg_->UpdateRHS(actSol);

          // inner forces due to nonlin formulation
          assemble_->AssembleNonLinRHS( actTime_ );

          Vector<Double> actRHS;
          Double *rhsPtr;
          algsys_->GetRHSVal( rhsPtr );
          StoreAlgsysToVec( actRHS, rhsPtr );
          
          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
        }

        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr    = residualL2Norm /  RhsLinL2Norm;
        else
          residualErr    = residualL2Norm;

        // calculate incremental error
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm = actSol.NormL2();
        Double incrementalErr;
      
        if ( actSolL2Norm > 1.0)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        // --------------------------------------------------------------------
        // output of norms and data
        // --------------------------------------------------------------------
        if ( nonLinLogging_ == true ) {
          Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                                incrementalErr, etaLineSearch);
        }

        // boolean variable, holds condition if another iteration step
        // is necessary
        performOneMoreStep = 
          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
      
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    } // load step loop

    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }


  void StdSolveStep::StepTransNonLinMaterial() {

    ENTER_FCN( "StdSolveStep::StepTransNonLinMaterial" );

    Double *solPtr;
    Double *incPtr;

    bool performOneMoreStep;

    Vector<Double> solInc( numEqns_ );
    Vector<Double> actSol( numEqns_ );
    actSol.Init();

    if (false){

      UInt length = algsys_->GetSolutionVal( solPtr );
      
      StoreAlgsysToVec(actSol, solPtr );
      
    
      StdVector<SinglePDE*> ptSinglePDEs;
      PDECoupling *ptCoupling;     //!< pointer to coupling object
      ptCoupling = PDE_.GetCoupling();

      //account for RHS
      assemble_->AssembleLinRHS( actTime_ );
      PDE_.ComputeRHS( actTime_ );

      // perform predictor step
      if ( TS_alg_== NULL ) {
        Error( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin",
               __FILE__, __LINE__ );
      }
      else {
        //compute predictors
        TS_alg_->Predictor(actSol);
      }

      Double loadFactor = 1.0;

      // compute u_{n+1}^0

      // to incorporate loads 
      Double RhsLinL2Norm = SetLinRHS(loadFactor); 

      assemble_->AssembleMatrices();
      if (assemble_->IsMatrixUpdated() ) {
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
      }

      // move mass and damping to RHS
      TS_alg_->UpdateRHS(actSol);

      //! account for Dirichlet BCs
      PDE_.SetBCs( actTime_ );
      algsys_->BuildInDirichlet();

      if (assemble_->IsMatrixUpdated() ) {
        SETPROFILE("Before SetupPrecond");
        algsys_->SetupPrecond( );
        SETPROFILE("After SetupPrecond / Before SetupSolver");
        algsys_->SetupSolver( );
        SETPROFILE("After SetupSolver / Before Solve");  
      }
    
      //     algsys_->SetupPrecond();
      //     algsys_->SetupSolver();
      algsys_->Solve();

      // retrieve and store solution
      length = algsys_->GetSolutionVal( solPtr );
      PDE_.SaveSolution(solPtr,length);

      StoreAlgsysToVec(actSol, solPtr );

      // stores solution
      //     for (UInt i=0;i<actSol.GetSize();i++)
      //       *(solPtr+i)=actSol[i];
      //     PDE_.SaveSolution(solPtr,length);

    }
    
    // set iteration counter
    UInt iterationCounter=0;
    Double RhsLinL2Norm;
    Vector<Double> uOld;
    
    Double *rhsPtr;
    Vector<Double> actRHS;
     
    StepTransLin();

    UInt length = algsys_->GetSolutionVal( solPtr );
    PDE_.SaveSolution(solPtr,length);
    
    StoreAlgsysToVec(actSol, solPtr );
    
    // to incorporate loads 
    Double loadFactor = 1.0;
    RhsLinL2Norm = SetLinRHS(loadFactor); 
    
    
    do {
      uOld=actSol;
      // compute u_{n+1}^k+1
      iterationCounter++;

      // re initialize RHS and system matrix
      algsys_->InitRHS();

      assemble_->AssembleLinRHS( actFreq_ );
      //    algsys_->InitMatrix();

      assemble_->AssembleMatrices();

      // account for Dirichlet BCs
      PDE_.SetBCs( actTime_ );

      //   algsys_->RemoveIDBCInfoFromMatrix();      

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
           
      algsys_->BuildInDirichlet();

      // put mass and damping on RHS
      TS_alg_->UpdateRHS(actSol);

      //    algsys_->GetRHSVal( rhsPtr );

      algsys_->RemoveIDBCInfoFromMatrix();      

      // substract K^* u^k from RHS
      TS_alg_->SubstractStiffnessFromRHS(actSol);

      algsys_->SetupPrecond();
      algsys_->SetupSolver();
      algsys_->Solve();

      // new solution is only an increment of the full solution =============
      UInt length = algsys_->GetSolutionVal( incPtr );
      StoreAlgsysToVec(solInc, incPtr );

      Double residualL2Norm;
      Double etaLineSearch = 1.0;     

      residualL2Norm = PDE_.GetRhsL2Norm(solInc);

      if ( lineSearch_ == "none" ) {
        actSol += solInc;
      }
      else {
        residualL2Norm = LineSearchMaterial(solInc, actSol, etaLineSearch, RhsLinL2Norm);
        //        std::cout<<"-- eta from LineSearch = "<< etaLineSearch <<std::endl;
      }
         
      residualL2Norm = PDE_.GetRhsL2Norm(solInc);     

      PDE_.SaveSolution(actSol.GetPointer(),length);

      Vector<Double> actRHS;
      //  Double *rhsPtr;
      algsys_->GetRHSVal( rhsPtr );
      StoreAlgsysToVec( actRHS, rhsPtr );
      
      
      Vector<Double> u_uOld(uOld.GetSize());
      u_uOld.Init();
      for (UInt ii=0;ii<uOld.GetSize();ii++){
        if(uOld[ii]!=0)
          u_uOld[ii]=(actSol[ii]-uOld[ii])/uOld[ii];
        
      }
      Double incrementL2Norm = u_uOld.NormL2();
      std::cout<<"-- residual2Norm = " << residualL2Norm <<", incrementL2Norm = "<<incrementL2Norm<< std::endl;    
      
       
      Double residualErr;
      if ( RhsLinL2Norm > 1.0 )
        residualErr    = residualL2Norm /  RhsLinL2Norm;
      else
        residualErr    = residualL2Norm;
       
      // calculate incremental error
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm = actSol.NormL2();
      Double incrementalErr;

      if ( actSolL2Norm > 1.0)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;

            
      // --------------------------------------------------------------------
      // output of norms and data
      // --------------------------------------------------------------------
      if ( nonLinLogging_ == true ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                              incrementalErr, etaLineSearch);
      }
      
      // boolean variable, holds condition if another iteration step
      // is necessary
      performOneMoreStep = 
        (incrementL2Norm > incStopCrit_)||(residualErr > residualStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
  
    // perform corrector step  
    TS_alg_->Corrector(actSol);

  }


  void StdSolveStep::StepTransNonLinHysteresis() {

    ENTER_FCN( "StdSolveStep::StepTransNonLinHysteresis" );

    Double *solPtr;
    bool performOneMoreStep;

    Vector<Double> solInc( numEqns_ );
    Vector<Double> actSol( numEqns_ );
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> prevSol( numEqns_ );
    actSol.Init();
  
    // get actual solution
    UInt length = algsys_->GetSolutionVal( solPtr );      
    StoreAlgsysToVec(actSol, solPtr );

    //solution from previous time step
    prevSol = actSol;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Error( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin",
             __FILE__, __LINE__ );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( actTime_ );

    // currently just for testing!!
    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      Info->PrintF(pdename_, "\n");
      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side 
      Double RhsLinL2Norm = SetLinRHS(loadFactor); 

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS( actTime_ );  

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS();

      // set iteration counter
      UInt iterationCounter=0;

      do {
	//std::cout << "\n Iteration counter: " <<  iterationCounter << std::endl;

        iterationCounter++;
	oldSol = actSol;

	//	std::cout << "Old Solution:\n" << oldSol << std::endl;

	//        RHS is already set up!!
	if ( iterationCounter > 0 ) {
	  // setup linear part of right hand side 
	  algsys_->InitRHS(RhsLinVal_.GetPointer());

	  // inner forces due to nonlin formulation
	  assemble_->AssembleNonLinRHS( actTime_ );  

	  //Update RHS (mass matrix on right hand side)
	  TS_alg_->UpdateRHS();
	}

        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();
  
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        algsys_->Solve();

        // new solution is only an increment of the full solution =============
	length = algsys_->GetSolutionVal( solPtr );      
	StoreAlgsysToVec(solInc, solPtr );

        Double residualL2Norm;
        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none" ) {
          actSol = solInc;
        }
        else {
	  Error("Currently lineSreach not supported",__FILE__,__LINE__);
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
        }

	//	std::cout << "New Solution:\n" << actSol << std::endl;
      
        //store actual solution to the solution-object sol_
	PDE_.SaveSolution(actSol.GetPointer(),length); 

        if ( lineSearch_ == "none" ) {
          // calculation of error norms
          // recalculate RHS with new values to get new residual (f^(k+1))=======
          algsys_->InitRHS(RhsLinVal_.GetPointer());

          //substract from RHS: intFactro*MASS*acc - intFactor*DAMP*vel
	  TS_alg_->UpdateRHS(actSol);

	  // substracte stiff-matrix from RHS
	  TS_alg_->SubstractStiffnessFromRHS(actSol);

          // inner forces due to nonlin formulation
          //assemble_->AssembleNonLinRHS( actTime_ );

	  Vector<Double> actRHS;
          Double *rhsPtr;
	  algsys_->GetRHSVal( rhsPtr );
	  StoreAlgsysToVec( actRHS, rhsPtr );
	  //	  std::cout << "actRHS:\n" << actRHS << std::endl;

          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
        }

        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr    = residualL2Norm /  RhsLinL2Norm;
        else
          residualErr    = residualL2Norm;

	//	std::cout << "Residual error: " << residualErr << std::endl;

	residualErr = 0.0;

        // calculate incremental error

	// compute L2-Norm of error between last incremental solution and
	// actual incremental solution
	Double solIncrL2Norm=0;
	for (UInt i=0; i<actSol.GetSize(); i++)
	  solIncrL2Norm += (actSol[i]-oldSol[i])*(actSol[i]-oldSol[i]);
	
	solIncrL2Norm = sqrt(solIncrL2Norm);
	Double actSolL2Norm = actSol.NormL2();
	
	Double incrementalErr;
	if (actSolL2Norm > 1)
	  incrementalErr = solIncrL2Norm / actSolL2Norm;
	else
	  incrementalErr = solIncrL2Norm;

        // --------------------------------------------------------------------
        // output of norms and data
        // --------------------------------------------------------------------
        if ( nonLinLogging_ == true ) {
          Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                                incrementalErr, etaLineSearch);
        }

        // boolean variable, holds condition if another iteration step
        // is necessary
        performOneMoreStep = 
          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
      
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    } // load step loop

    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }


  void StdSolveStep::PostStepTrans( ) {

    ENTER_FCN( "StdSolveStep::PostStepTrans" );

    if ( PDE_.GetFracDamping() ) {
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

      // Following method is essential for fractional damping model
      TS_alg_->AdvanceTimestep(solHelp);
    }
  }



  // ======================================================
  // Solve Step Harmonic  SECTION  
  // ======================================================w

  void StdSolveStep::PreStepHarmonic() {
    
    ENTER_FCN( "StdSolveStep::PreStepHarmonic" );
  
    algsys_->InitRHS();
  }


  void StdSolveStep::SolveStepHarmonic() {

    ENTER_FCN( "StdSolveStep::SolveStepHarmonic" );

    if ( nonLin_ ) {
      StepHarmonicNonLin();
    }
    else {
      StepHarmonicLin();
    }
  }


  void StdSolveStep::StepHarmonicLin() {

    ENTER_FCN( "StdSolveStep::StepHarmonicLin" );

    //this has to be done each frequency!
    assemble_->AssembleLinRHS( actFreq_ );

    if ( PDE_.IsComputeRHS4HarmSet() ) {
      // Evaluating RHS with nodal srcs for harmonic flownoise problems
      std::cout<<"Calling ComputeRHS from StepHarmonicLin()"<<std::endl;
      PDE_.ComputeRHS(actFreq_);
    }

    assemble_->AssembleMatrices( );
    
    PDE_.SetBCs( actFreq_ );

    // store rhs vector back to PDE
    Complex * rhsPt;
    UInt length = algsys_->GetRHSVal(rhsPt);
    PDE_.SaveRHS( rhsPt, length );
    
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }
    
    algsys_->BuildInDirichlet();
 
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }

    algsys_->Solve();

    Complex* ptSol = NULL;
    length    =  algsys_->GetSolutionVal(ptSol);
    PDE_.SaveSolution(ptSol,length);
  }


  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================
  
  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {
    ENTER_FCN( "StdSolveStep::CalcEigenFrequencies" );
    
    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();
    
    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, false );

    // Calculate eigenfrequencies
    const Double * val = NULL;
    const Double * err = NULL;
    UInt numConverged = algsys_->CalcEigenFrequencies( val, err);

    // Copy eigenvalues into vector
    frequencies.Resize( numConverged );
    errBounds.Resize( numConverged );
    for ( UInt i = 0; i < numConverged; i++ ) {
      frequencies[i] = val[i];
      errBounds[i] = err[i];
    }

    return numConverged;
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {
    ENTER_FCN( "StdSolveStep::CalcEigenFrequencies<Complex>" );
    
    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();
    
    assemble_->AssembleMatrices();
    
    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, true );
    
    // Calculate eigenfrequencies
    const  Complex* val = NULL;
    const Double * err = NULL;
    UInt numConverged = algsys_->CalcEigenFrequencies( val, err);
    
    // Copy eigenvalues into vector
    frequencies.Resize( numConverged );
    errBounds.Resize( numConverged );
    for ( UInt i = 0; i < numConverged; i++ ) {
      frequencies[i] = val[i];
      errBounds[i] = err[i];
    }
    
    return numConverged;
  }

  void StdSolveStep::CalcEigenMode( UInt numMode ) {
    
    ENTER_FCN( "StdSolveStep::CalcEigenMode" );

    Double * ptSol = NULL;
    UInt size = 0;
    
    algsys_->CalcEigenMode( numMode );

    // Get the solution and store it
    size = algsys_->GetSolutionVal(ptSol);
    PDE_.SaveSolution(ptSol,size);
  }


  

  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor)
  {
    ENTER_FCN( "StdSolveStep::SetLinRHS" );

    Double RhsLinL2Norm;  


    // to incorporate loads
    assemble_->AssembleLinRHS(actTime_); 

    // Stores rhs vector into extForces and returns that L2-norm
 
    Double *solPtr;
    algsys_->GetRHSVal( solPtr );
    StoreAlgsysToVec(RhsLinVal_, solPtr );

    RhsLinVal_ *= loadFactor;

    RhsLinL2Norm = RhsLinVal_.NormL2();
 
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      Warning("Zero external force vector!! ", __FILE__,__LINE__);
    }

    return RhsLinL2Norm;
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void StdSolveStep::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "StdSolveStep::StoreAlgsysToVec" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    vec.Resize( numEqns_ );
    vec.Init();
      
    for (UInt i=0; i<numEqns_; i++) {
      vec[i] = pt[i];
    }
  }



  Double StdSolveStep::LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                                  Double& etaLineSearch, bool trans)
  {
    ENTER_FCN( "StdSolveStep::LineSearch" );

    Vector<Double> solOld(actSol);
    const UInt nrEtas = 5;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125, 0.1};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      actSol = solIncrement * eta[i];
      actSol += solOld;

      //store new solution
      sol_->SetAlgSysVector(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS(RhsLinVal_.GetPointer());

      if( trans ) {
        assemble_->AssembleNonLinRHS( actTime_ );
        TS_alg_->UpdateRHS(actSol);
      }
      else {
        assemble_->AssembleNonLinRHS( actTime_);
      }


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      Double *rhsPtr;
      algsys_->GetRHSVal( rhsPtr );
      StoreAlgsysToVec(actRHS, rhsPtr );

      // calculation of residual error =======================================
      Double residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;
    
    actSol  = solIncrement * etaOpt;
    actSol += solOld;

    return residualL2NormOpt;
  }


  Double StdSolveStep::LineSearchMaterial(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                                          Double& etaLineSearch, Double& RhsLinL2Norm, bool trans)
  {
    ENTER_FCN( "StdSolveStep::LineSearch" );

    Vector<Double> solOld(actSol);
    const UInt nrEtas = 3;
    const Double eta[nrEtas] = {0.9, 0.5, 0.3};
    //    const Double eta[nrEtas] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;
    Double incrementalErrOpt = 1e+15;

    Double *solPtr;

    UInt length = algsys_->GetSolutionVal( solPtr );
    
    for( UInt i=0; i<nrEtas; i++) {
//       if (i>0)
//         actSol=solOld;
      actSol = solOld + solIncrement * eta[i];
      //      actSol -= solOld;

  //     std::cout<<"actSol in LineSearch" <<std::endl;
//       std::cout<<actSol<<std::endl;
      //store new solution    
      for (UInt ii=0;ii<actSol.GetSize();ii++)
        *(solPtr+ii)=actSol[ii];
      PDE_.SaveSolution(solPtr,length);
      
      // Recalculate residual, f-Cu-Mu-K*u
      algsys_->InitRHS();

      assemble_->AssembleLinRHS( actFreq_ );

      // assemble!
      assemble_->AssembleMatrices();

      // account for Dirichlet BCs
      PDE_.SetBCs( actTime_ );

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      algsys_->BuildInDirichlet();

      TS_alg_->UpdateRHS(actSol);
      // substract K^* u^k from RHS

      algsys_->RemoveIDBCInfoFromMatrix();      

      TS_alg_->SubstractStiffnessFromRHS(actSol);
      
      
      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      Double *rhsPtr;
      algsys_->GetRHSVal( rhsPtr );
      StoreAlgsysToVec(actRHS, rhsPtr );

      // calculation of residual error =======================================
      Double residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  (f-Ku )

      Double residualErr;

       if ( RhsLinL2Norm > 1.0 )
         residualErr    = residualL2Norm /  RhsLinL2Norm;
       else
         residualErr    = residualL2Norm;
      
      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;
   
    actSol  = solOld + solIncrement * etaOpt;
    
    return residualL2NormOpt;
  }


  // returns that L2-norm of an algsys vector
  Double StdSolveStep::AlgsysL2Norm(Double * pt)
  {
    ENTER_FCN( "StdSolveStep::AlgsysL2Norm" );

    //const UInt numElems = numPDENodes_ * dofspernode_;
    Double quadSum = 0;
  
    for (UInt i=0; i<numEqns_; i++)   
      quadSum += pt[i]*pt[i];

    return sqrt(quadSum);
  }

  // read nonlinear parameters from xml file
  void StdSolveStep::ReadNonLinData() {
    ENTER_FCN( "StdSolveStep::ReadNonLinData" );
    
    // Get ParamNode of pde
    ParamNode  * nonLinNode = NULL;
    std::string pdeName = PDE_.GetName();

    // SPECIAL: If PDE-name contains the "-" sign,
    // we assume that we have a coupledPDE
    if(boost::find_first(pdeName, "-")) {
      // go to direct-coupling List
      // iterate over all nodes and pick the one
      // which has a "nonLinear" element
      // get current multiSequenceStep
      UInt actMsStep = 
        domain->GetSingleDriver()->GetActSequenceStep();
      StdVector<ParamNode *> pairCouplings;
        ParamNode * cplList = 
        param->Get("sequenceStep","index", GenStr(actMsStep ) )
        ->Get("couplingList", false);
      if( cplList ) {
        StdVector<ParamNode *> pairCouplings = 
          cplList->Get("direct")->GetChildren();
        // look for each direct coupling, if it has a "nonLin"-node
        for( UInt iCpl = 0; iCpl < pairCouplings.GetSize(); iCpl++ ) {
          nonLinNode = pairCouplings[iCpl]->Get("nonLinear", false );
          if( nonLinNode ) {
            break;
          }
        }
      }
      
    } else {
      nonLinNode = PDE_.GetParamNode()->Get("nonLinear", false );
    }

    // Read data, if "nonLinear" element was found
    if( nonLinNode ) {

      // solution method
      nonLinNode->Get( "method", nonLinMethod_, false );

      // perform logging?
      nonLinNode->Get( "logging", nonLinLogging_, true );

      // type of line search
      nonLinNode->Get( "lineSearch")->Get( "type", lineSearch_ );

      // incremental stopping criterion
      nonLinNode->Get( "incStopCrit", incStopCrit_ );

      // residual stopping criterion
      nonLinNode->Get( "resStopCrit", residualStopCrit_ );

      // maximal number of NL-iterations
      nonLinNode->Get( "maxNumIters", nonLinMaxIter_ );
    }
    
  }
  

  void StdSolveStep::WriteClaNlNorms(const UInt iterationCounter, 
                                     const Double residualL2Norm, 
                                     const Double extForcesL2Norm,
                                     const Double residualErr, 
                                     const Double solIncrL2Norm, 
                                     const Double actSolL2Norm, 
                                     const Double incrementalErr)
  {
    ENTER_FCN( "StdSolveStep::WriteClaNlNorms" );
  
    *cla << std::endl << " ======================================================= " 
         << std::endl
         << " NONLINEAR ITERATION " << iterationCounter << std::endl
         << " ======================================================= " << std::endl;
    *cla << " === Residual norm:          " << residualL2Norm << std::endl;
    *cla << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
    *cla << "     Residual error          " << residualErr << std::endl;
  
    *cla << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
    *cla << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
    *cla << "     Incremental error       " << incrementalErr << std::endl;      
  }



} // end of namespace

