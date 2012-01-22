#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

#include "StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultCache.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/EvalIntegrals/BiotSavart.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"

#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField {

  StdSolveStep::StdSolveStep(StdPDE & apde)
    :BaseSolveStep(),
     PDE_(apde)
  {



    pdename_      = PDE_.GetName();
    isaxi_        = PDE_.GetIsaxi();
    subdoms_      = PDE_.getPDE_subdoms();
    materialData_ = PDE_.getPDEMaterialData();
    ptgrid_       = PDE_.getPDE_grid();
    algsys_       = PDE_.getPDE_algsys();
    assemble_     = PDE_.getPDE_assemble();

    results_      = PDE_.GetResultInfos();

    // copy FE functions of PDE
    feFunctions_ = PDE_.GetFeFunctions();
    rhsFeFunctions_ = PDE_.GetRhsFeFunctions();
    
    // Copy vectors FE functions in SBM-vector for communication
    // with OLAS and time stepping
    solVec_.SetSize( feFunctions_.size() );
    rhsVec_.SetSize( feFunctions_.size() );
    


    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it;
    it = feFunctions_.begin();
    UInt pos = 0;
    for( ; it != feFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      solVec_.SetSubVector(ptFct->GetSingleVector(), pos++);
    }
    pos = 0;
    it = rhsFeFunctions_.begin();
    for( ; it != rhsFeFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      rhsVec_.SetSubVector(ptFct->GetSingleVector(), pos++);
    }
    // Make sure to have both vectors as "weak" vectors,
    // as the feFunctions themselves are responsible for
    // creation and destruction.
    solVec_.SetOwnership(false);
    rhsVec_.SetOwnership(false);
    
    // nonlinear parameters
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;

    nonLin_            = PDE_.IsNonLin();
    nonLinMaterial_    = PDE_.IsNonLinMaterial();
    isHyst_            = PDE_.IsHysteresis();
    totalFormulation_  = PDE_.IsTotaFormulation();
    regionNonLinTypes_ = PDE_.GetNonLinRegionTypes();

    startStep_ = 1;
    
    // set entry type of SBM vector
    oldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    tmpOldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    DeltaRhsLinVal_ =  SBM_Vector(BaseMatrix::DOUBLE);
    RhsLinVal_ =  SBM_Vector(BaseMatrix::DOUBLE);

    // In the end, read nonlinear data from xml-file
    if( nonLin_ || nonLinMaterial_ ) {
      ReadNonLinData();
    }
    
    //logFile_.open("nonlin.txt");

    mHandle_ = domain->GetMathParser()->GetNewHandle();
    mParser_ = domain->GetMathParser();
    mParser_->SetExpr(mHandle_,"step");
  }


  //! Destructor
  StdSolveStep::~StdSolveStep() {
    //logFile_.close();
    mParser_->ReleaseHandle(mHandle_);
  }


  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================

  void StdSolveStep::PreStepStatic( ) {

    // init RHS at this place, because e.g. forces of other PDEs are added
    // to RHS afterwards
    algsys_->InitRHS();

  }


  void StdSolveStep::PostStepStatic() {

    // increment coupling counter
    if ( PDE_.IsIterCoupled() ) {
      UInt& counter = PDE_.GetIterCoupledCounter();
      counter++;
    }
    
    // check for Biot Savart
    if ( PDE_.IsBiotSavart() ) {
      WARN("BIot Savart not yet adpated to new structure");
//      Vector<Double> & sol = 
//          dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
//      Vector<Double>& magVecBiotSavart = 
//          PDE_.GetBiotSavart()->CalcFieldAllEqns(false); 
//      sol += magVecBiotSavart;
    }

  }


  void StdSolveStep::SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams) {

    if (nonLin_) {
      if ( totalFormulation_ )
        StepStaticNonLinTotal(analysis_id);
      else
        StepStaticNonLin(analysis_id);
    }
    else {
      StepStaticLin(analysis_id, adjointParams);
    }
  }


  void StdSolveStep::StepStaticLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) {

    assemble_->AssembleMatrices();

    // The RHS-sources and boundary conditions
    // have to be reassembled each time
    assemble_->AssembleLinRHS(adjointParams);
    PDE_.SetBCs();

    // store rhs vector back to PDE
    algsys_->GetRHSVal(rhsVec_);

    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }

    // Incorporate Boundary conitions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond(analysis_id);

      algsys_->SetupSolver(analysis_id);
    }

    // Solve problem
    algsys_->Solve(analysis_id);

    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }


  void StdSolveStep::StepStaticNonLin(PtrParamNode analysis_id)
  {
    REFACTOR;

  }


  void StdSolveStep::StepStaticNonLinTotal(PtrParamNode analysis_id)
  {
    REFACTOR;

  }


  // ======================================================
  // Solve Step Transient SECTION
  // ======================================================

  void StdSolveStep::InitTimeStepping(){
    //also initialize vectors for the time stepping scheme

    stageRHS_.Resize(feFunctions_.size());

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    UInt pos = 0;
    UInt rhsSize = 0;

    //reserve memory for the rhs
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++,pos++){
      rhsSize = fncIt->second->GetSingleVector()->GetSize();
      stageRHS_.SetSubVector(new Vector<Double>(),pos);
      stageRHS_.GetPointer(pos)->Resize(rhsSize);
    }
    stageRHS_.Init();

  }

  void StdSolveStep::PreStepTrans() {
    ResultCache::SetStepValue( actTime_ );

    // due to coupling-pdes, the RHS has to be initialized BEFORE
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();

    UInt step = (UInt) mParser_->Eval(mHandle_);

    PDE_.ReadDisplacementAndUpdateGrid( step );



  }


  void StdSolveStep::SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams) {


    // do a nonlinear material time step
    if ( nonLin_ && nonLinMaterial_ ) {
      StepTransNonLinMaterial(analysis_id);
    }
    else if ( isHyst_ ) {
      StepTransNonLinHysteresis(analysis_id);
    }
    // do a nonlinear time step
    else if (nonLin_){
      if ( totalFormulation_ )
        StepTransNonLinTotal(analysis_id);
      else
        StepTransNonLin(analysis_id);
    }
    // do a linear time step
    else {
      StepTransLin(analysis_id, adjointParams);
    }
  }


  void StdSolveStep::StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) 
  {
    //TODO: add consistency check here
    //basically loop over all functions and check if the solution order is the same...

    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();

    std::map<FEMatrixType,Integer>::iterator matIt;

    UInt pos = 0;
    bool effectiveMatrixUpdated = false;

    for(UInt i=0;i<numStages;i++){
      effectiveMatrixUpdated = false;

      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++,pos++){
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
      }
      stageSol.SetOwnership(false);
      stageRHS_.Init();

      algsys_->InitRHS();

      //account for RHS
      assemble_->AssembleLinRHS(adjointParams);
      PDE_.ComputeRHS( actTime_ );

      // store rhs vector back to PDE for e.g. postProc
      algsys_->GetRHSVal(rhsVec_);

      assemble_->AssembleMatrices();
      if(assemble_->IsMatrixUpdated()){
        matrix_factor_.clear();
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_);
        }
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        effectiveMatrixUpdated = true;
      }

      //now compute the effective right hand side
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++,pos++){
          fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_);
      }

      // set boundary conditions
      PDE_.SetBCs();
      algsys_->BuildInDirichlet();

      if( effectiveMatrixUpdated ){
        algsys_->SetupPrecond( analysis_id);
        algsys_->SetupSolver(analysis_id);
      }

      algsys_->Solve(analysis_id);

      algsys_->GetSolutionVal(stageSol);
    }

    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
      fncIt->second->GetTimeScheme()->FinishStep();
    }

  }



  void StdSolveStep::StepTransNonLin(PtrParamNode base_analysis_id) {

    REFACTOR;
/*
    bool performOneMoreStep;
    SBM_Vector solInc(BaseMatrix::DOUBLE);

    //get actual solution
    algsys_->GetSolutionVal( actSol_ );

    // perform predictor step
    if ( TS_alg_== NULL ) {
      EXCEPTION ( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin" );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol_);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs();

    // currently just for testing!!
    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS();  

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(actSol_);

      // set iteration counter
      UInt iterationCounter=0;


      do {
        iterationCounter++;

        // RHS is already set up!!
        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(base_analysis_id, "load", iload, "nonLin", iterationCounter);

        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();

        algsys_->SetupPrecond(child_id);
        algsys_->SetupSolver(child_id);
        algsys_->Solve(child_id);

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solInc );

				// initialize residualL2Norm or receive compiler warning
        Double residualL2Norm = 0.0;
        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none" ) {
          actSol_.Add(1.0, solInc);
        }
        else {
          residualL2Norm = LineSearch(solInc, actSol_, etaLineSearch, true);
        }

        //store A_(n+1) in the solution-object sol_
        solVec_ = actSol_;

        if ( lineSearch_ == "none" ) {
          // calculation of error norms
          // recalculate RHS with new values to get new residual (f^(k+1))=======
          algsys_->InitRHS(RhsLinVal_);

          //Update RHS (mass matrix on right hand side)
          TS_alg_->UpdateRHS(actSol_);

          // inner forces due to nonlin formulation
          assemble_->AssembleNonLinRHS();

          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
          residualL2Norm = actRHS.NormL2();
        }

        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr    = residualL2Norm /  RhsLinL2Norm;
        else
          residualErr    = residualL2Norm;

        // calculate incremental error
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm = actSol_.NormL2();
        Double incrementalErr;

        if ( actSolL2Norm > 1.0)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        // --------------------------------------------------------------------
        // output of norms and data
        // --------------------------------------------------------------------
        if ( nonLinLogging_ == true )
          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);

        // boolean variable, holds condition if another iteration step
        // is necessary
        performOneMoreStep =
          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } // load step loop

    //perform corrector step
    TS_alg_->Corrector(actSol_);
    */
  }


  void StdSolveStep::StepTransNonLinTotal(PtrParamNode base_analysis_id) {
    REFACTOR;
/*
    std::cout << "\n In :StepTransNonLinTotal  \n " << std::endl;
 
    bool performOneMoreStep;
    UInt iterationCounter=0;
    Double etaLineSearch;
    Double incrementalErr, residualErr;

    SBM_Vector newSol(BaseMatrix::DOUBLE), 
               oldSol(BaseMatrix::DOUBLE), solPrev(BaseMatrix::DOUBLE);
  
    // remember previous solution  
    solPrev = solVec_;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Exception( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin");
    }
    else {
      //compute predictors
      TS_alg_->Predictor(solPrev);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( );    

    Double loadFactor = 0.0;
    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
//      Info->PrintF(pdename_, "\n");
//      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS();  

      do {
        iterationCounter++;
        
        // set solution of previous iteration
        if (iterationCounter == 1 ) {
          oldSol = solPrev;
        }
        else {
          oldSol = newSol;
        }
        
        //init RHS with linear part!
        algsys_->InitRHS(RhsLinVal_);
        
        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(base_analysis_id, "load", iload, "nonLin", iterationCounter);
        
        //perform new assembly
        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        
        
        //Update RHS (mass matrix on right hand side)
        TS_alg_->UpdateRHS();
        
        // build in the Dirichlet vales in system matrix and rhs
        algsys_->BuildInDirichlet();
        
        algsys_->SetupSolver(child_id);
        algsys_->SetupPrecond(child_id);
        
        algsys_->Solve(child_id);
        algsys_->GetSolutionVal(newSol); 
        
        //store solution for (n+1)
        solVec_ = newSol;
        
        // compute L2-Norm of error between last incremental solution and
        // actual incremental solution
        Double solIncrL2Norm=0;
        SBM_Vector temp(BaseMatrix::DOUBLE);
        temp.Add(1.0, newSol, -1.0, oldSol );
//        for (UInt i=0; i<newSol.GetSize(); i++)
//          solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
//        
//        solIncrL2Norm = sqrt(solIncrL2Norm);
        solIncrL2Norm = temp.NormL2();
        Double actSolL2Norm = newSol.NormL2();
        
        if (actSolL2Norm > 1)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        //just dummy things
        etaLineSearch = RhsLinL2Norm;
        etaLineSearch = 1.0;
        residualErr   = incrementalErr;
        
        // output of norms and data
        nonLinLogging_ = true;
        if ( nonLinLogging_ == true )
          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        
        performOneMoreStep = (incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
        
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
    } // load step loop

    if ( incrementalErr > incStopCrit_ ) 
      std::cout << "Not converged, norm is: " <<incrementalErr << std::endl; 

    //perform corrector step  
    TS_alg_->Corrector(newSol);
    */
  }

  
  void StdSolveStep::StepTransNonLinMaterial(PtrParamNode analysis_id) {

    REFACTOR;
//
//    bool performOneMoreStep;
//
//    SBM_Vector solInc, actSol;
//    actSol.Init();
//
//    // set iteration counter
//    UInt iterationCounter=0;
//    Double RhsLinL2Norm;
//    SBM_Vector uOld, actRHS;
//
//    StepTransLin(analysis_id);
//    algsys_->GetSolutionVal( actSol );
//    PDE_.SaveSolution( actSol );
//
//    // to incorporate loads
//    Double loadFactor = 1.0;
//    RhsLinL2Norm = SetLinRHS(loadFactor);
//
//
//    do {
//      uOld=actSol;
//      // compute u_{n+1}^k+1
//      iterationCounter++;
//
//      PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "nonLin", iterationCounter);
//      
//      // re initialize RHS and system matrix
//      algsys_->InitRHS();
//
//      assemble_->AssembleLinRHS();
//
//      assemble_->AssembleMatrices();
//
//      // account for Dirichlet BCs
//      PDE_.SetBCs();
//
//      algsys_->ConstructEffectiveMatrix(matrix_factor_);
//
//      algsys_->BuildInDirichlet();
//
//      // put mass and damping on RHS
//      TS_alg_->UpdateRHS(actSol);
//
//      algsys_->RemoveIDBCInfoFromMatrix();
//
//      // substract K^* u^k from RHS
//      TS_alg_->SubstractStiffnessFromRHS(actSol);
//
//      algsys_->SetupPrecond(analysis_id);
//      algsys_->SetupSolver(analysis_id);
//      algsys_->Solve(analysis_id);
//
//      // new solution is only an increment of the full solution =============
//      algsys_->GetSolutionVal( solInc );
//      Double residualL2Norm;
//      Double etaLineSearch = 1.0;
//
//      residualL2Norm = solInc.NormL2();
//
//      if ( lineSearch_ == "none" ) {
//        actSol.Add( 1.0, solInc );
//      }
//      else {
//        residualL2Norm = LineSearchMaterial(solInc, actSol, etaLineSearch, RhsLinL2Norm);
//      }
//
//      residualL2Norm = solInc.NormL2();
//
//      PDE_.SaveSolution( actSol );
//
//      SBM_Vector actRHS;
//      algsys_->GetRHSVal( actRHS );
//
//      Vector<Double> u_uOld(uOld.GetSize());
//      u_uOld.Init();
//      for (UInt ii=0;ii<uOld.GetSize();ii++){
//        if(uOld[ii]!=0)
//          u_uOld[ii]=(actSol[ii]-uOld[ii])/uOld[ii];
//
//      }
//      Double incrementL2Norm = u_uOld.NormL2();
//      std::cout<<"-- residual2Norm = " << residualL2Norm
//               <<", incrementL2Norm = "<<incrementL2Norm<< std::endl;
//
//      Double residualErr;
//      if ( RhsLinL2Norm > 1.0 )
//        residualErr    = residualL2Norm /  RhsLinL2Norm;
//      else
//        residualErr    = residualL2Norm;
//
//      // calculate incremental error
//      Double solIncrL2Norm = solInc.NormL2();
//      Double actSolL2Norm = actSol.NormL2();
//      Double incrementalErr;
//
//      if ( actSolL2Norm > 1.0)
//        incrementalErr = solIncrL2Norm / actSolL2Norm;
//      else
//        incrementalErr = solIncrL2Norm;
//
//
//      // --------------------------------------------------------------------
//      // output of norms and data
//      // --------------------------------------------------------------------
//      if ( nonLinLogging_ == true )
//        WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
//
//      // boolean variable, holds condition if another iteration step
//      // is necessary
//      performOneMoreStep =
//        (incrementL2Norm > incStopCrit_)||(residualErr > residualStopCrit_);
//
//    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
//
//    // perform corrector step
//    TS_alg_->Corrector(actSol);

  }


  void StdSolveStep::StepTransNonLinHysteresis(PtrParamNode analysis_id) {

    REFACTOR;
    //    bool performOneMoreStep;
//
//    Vector<Double> solInc( numEqns_ );
//    Vector<Double> actSol( numEqns_ );
//    Vector<Double> oldSol( numEqns_ );
//    Vector<Double> prevSol( numEqns_ );
//    actSol.Init();
//
//    // get actual solution
//    algsys_->GetSolutionVal( actSol );
//
//    //solution from previous time step
//    prevSol = actSol;
//
//    // perform predictor step
//    if ( TS_alg_== NULL ) {
//      EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin" );
//    }
//    else {
//      //compute predictors
//      TS_alg_->Predictor(actSol);
//    }
//
//    //! account for Dirichlet BCs
//    PDE_.SetBCs();
//
//    // currently just for testing!!
//    // loop over load factor
//    Double loadFactor = 0.0;
//
//    for ( UInt iload=0; iload<1; iload++ ) {
//      loadFactor += 1.0;
//      info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);
//
//      // setup right hand side
//      Double RhsLinL2Norm = SetLinRHS(loadFactor);
//
//      // inner forces due to nonlin formulation
//      assemble_->AssembleNonLinRHS();  
//
//      //Update RHS (mass matrix on right hand side)
//      TS_alg_->UpdateRHS();
//
//      // set iteration counter
//      UInt iterationCounter=0;
//
//      do {
//
//        iterationCounter++;
//        oldSol = actSol;
//
//        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "load", iload, "nonLin", iterationCounter);
//
//        //        RHS is already set up!!
//        if ( iterationCounter > 0 ) {
//          // setup linear part of right hand side
//          algsys_->InitRHS(RhsLinVal_);
//
//          // inner forces due to nonlin formulation
//          assemble_->AssembleNonLinRHS();  
//
//          //Update RHS (mass matrix on right hand side)
//          TS_alg_->UpdateRHS();
//        }
//
//        assemble_->AssembleMatrices();
//        algsys_->ConstructEffectiveMatrix(matrix_factor_);
//        algsys_->BuildInDirichlet();
//
//        algsys_->SetupPrecond(analysis_id);
//        algsys_->SetupSolver(analysis_id);
//        algsys_->Solve(analysis_id);
//
//        // new solution is only an increment of the full solution =============
//        algsys_->GetSolutionVal( solInc );
//
//        Double residualL2Norm;
//        Double etaLineSearch = 1.0;
//
//        if ( lineSearch_ == "none" ) {
//          actSol = solInc;
//        }
//        else {
//          EXCEPTION("Currently lineSreach not supported" );
//          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
//        }
//
//        //store actual solution to the solution-object sol_
//        PDE_.SaveSolution(actSol.GetPointer(), actSol.GetSize() );
//
//        if ( lineSearch_ == "none" ) {
//          // calculation of error norms
//          // recalculate RHS with new values to get new residual (f^(k+1))=======
//          algsys_->InitRHS(RhsLinVal_);
//
//          //substract from RHS: intFactro*MASS*acc - intFactor*DAMP*vel
//          TS_alg_->UpdateRHS(actSol);
//
//          // substracte stiff-matrix from RHS
//          TS_alg_->SubstractStiffnessFromRHS(actSol);
//
//          // inner forces due to nonlin formulation
//          //assemble_->AssembleNonLinRHS( actTime_ );
//
//          Vector<Double> actRHS;
//          algsys_->GetRHSVal( actRHS );
//
//          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
//          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
//        }
//
//        Double residualErr;
//        if ( RhsLinL2Norm > 1.0 )
//          residualErr    = residualL2Norm /  RhsLinL2Norm;
//        else
//          residualErr    = residualL2Norm;
//
//
//        residualErr = 0.0;
//
//        // calculate incremental error
//
//        // compute L2-Norm of error between last incremental solution and
//        // actual incremental solution
//        Double solIncrL2Norm=0;
//        for (UInt i=0; i<actSol.GetSize(); i++)
//          solIncrL2Norm += (actSol[i]-oldSol[i])*(actSol[i]-oldSol[i]);
//
//        solIncrL2Norm = sqrt(solIncrL2Norm);
//        Double actSolL2Norm = actSol.NormL2();
//
//        Double incrementalErr;
//        if (actSolL2Norm > 1)
//          incrementalErr = solIncrL2Norm / actSolL2Norm;
//        else
//          incrementalErr = solIncrL2Norm;
//
//        // --------------------------------------------------------------------
//        // output of norms and data
//        // --------------------------------------------------------------------
//        if ( nonLinLogging_ == true )
//          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
//
//        // boolean variable, holds condition if another iteration step
//        // is necessary
//        performOneMoreStep =
//          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
//
//      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
//
//    } // load step loop
//
//    //perform corrector step
//    TS_alg_->Corrector(actSol);
  }


  void StdSolveStep::PostStepTrans( ) {

//    WARN("Biot-Savart not yet included";)
//    Vector<Double> & solHelp =
//      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
//
//    // Following method is essential for fractional damping model
//    TS_alg_->AdvanceTimestep(solHelp);
//    
//    // check for Biot Savart
//    if ( PDE_.IsBiotSavart() ) {
//      Vector<Double> & sol = 
//          dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
//      Vector<Double>& magVecBiotSavart = 
//          PDE_.GetBiotSavart()->CalcFieldAllEqns(false);
//      sol += magVecBiotSavart;
//    }
  }



  // ======================================================
  // Solve Step Harmonic  SECTION
  // ======================================================w

  void StdSolveStep::PreStepHarmonic() {
    ResultCache::SetStepValue( actFreq_ );
    algsys_->InitRHS();
  }


  void StdSolveStep::SolveStepHarmonic(PtrParamNode analysis_id) {
    if ( nonLin_ ) {
      StepHarmonicNonLin(analysis_id);
    }
    else {
      StepHarmonicLin(analysis_id);
    }
  }


  void StdSolveStep::StepHarmonicLin(PtrParamNode analysis_id) {


    //this has to be done each frequency!
    assemble_->AssembleLinRHS();

    if ( PDE_.IsComputeRHS4HarmSet() ) {
      // Evaluating RHS with nodal srcs for harmonic flownoise problems
      PDE_.ComputeRHS(actFreq_);
    }

    assemble_->AssembleMatrices( );
    PDE_.SetBCs();

    // store rhs vector back to PDE
    algsys_->GetRHSVal( rhsVec_ );

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }

    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond(analysis_id);
      algsys_->SetupSolver(analysis_id);
    }

    algsys_->Solve(analysis_id);
    algsys_->GetSolutionVal(solVec_);
  }


  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, false );

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds);

    return frequencies.GetSize();
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, true );

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds );

    return frequencies.GetSize();
  }

  void StdSolveStep::CalcEigenMode( UInt numMode ) {


    algsys_->CalcEigenMode( numMode );

    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }




  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor)
  {

    Double RhsLinL2Norm;


    // to incorporate loads
    assemble_->AssembleLinRHS(); 

    // Stores rhs vector into extForces and returns that L2-norm
    algsys_->GetRHSVal( RhsLinVal_ );
    RhsLinVal_.ScalarMult(loadFactor);

    RhsLinL2Norm = RhsLinVal_.NormL2();

    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      WARN("Zero external force vector!! ");
    }

    return RhsLinL2Norm;
  }

  UInt StdSolveStep::SetDeltaLinRHS()
  {

    // to incorporate loads
    assemble_->AssembleLinRHS(); 
    
    SBM_Vector newRhsLinVal(BaseMatrix::DOUBLE); //!< external forces (for nonlin simulations)
    algsys_->GetRHSVal( newRhsLinVal );
    DeltaRhsLinVal_.Add( 1.0, newRhsLinVal,
                        -1.0, tmpOldRhsLinVal_ );
    
    RhsLinVal_=tmpOldRhsLinVal_;

    Double DeltaNorm=DeltaRhsLinVal_.NormL2();
    Double oldNorm  =tmpOldRhsLinVal_.NormL2();
    Double aux;
    UInt nrLoadSteps;

    UInt minNrLoadSteps = 1;
    UInt maxNrLoadSteps = 1;

    if(oldNorm<1e-13)
      nrLoadSteps=maxNrLoadSteps;
    else{
      aux= Double(maxNrLoadSteps)*(DeltaNorm/oldNorm);
      nrLoadSteps=UInt(aux);
    }
    if(nrLoadSteps<minNrLoadSteps)
      nrLoadSteps=minNrLoadSteps;
    else if (nrLoadSteps>maxNrLoadSteps)
      nrLoadSteps=maxNrLoadSteps;

    oldRhsLinVal_=tmpOldRhsLinVal_;
    tmpOldRhsLinVal_=newRhsLinVal;

    return nrLoadSteps;
  }



  Double StdSolveStep::LineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol,
                                  Double& etaLineSearch, bool trans)
  {
    REFACTOR;
/*
    SBM_Vector solOld;
    solOld = actSol;
    const UInt nrEtas = 5;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125, 0.1};
		// initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      actSol.Add( 1.0, solOld, eta[i], solIncrement);
//      actSol = solIncrement * eta[i];
//      actSol += solOld;

      //store new solution
      solVec_ = actSol;
      
      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS(RhsLinVal_ );

      if( trans ) {
        assemble_->AssembleNonLinRHS();
        TS_alg_->UpdateRHS(actSol);
      }
      else {
        assemble_->AssembleNonLinRHS();
      }


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2(); // L2Norm of  ( f_i^(k+1) - f_a )

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;
    
    // Careful: in the end, we have to re-assemble the RHS with the correct
    // value i.e. use the "optimal" solution
    actSol.Add( 1.0, solOld, etaOpt, solIncrement );

    return residualL2NormOpt;
    */
    return 0.0;
  }


  Double StdSolveStep::LineSearchMaterial(SBM_Vector& solIncrement, SBM_Vector& actSol,
                                          Double& etaLineSearch, Double& RhsLinL2Norm, bool trans)
  {
    REFACTOR;
/*
    SBM_Vector solOld;
    solOld = actSol;
    const UInt nrEtas = 3;
    const Double eta[nrEtas] = {0.9, 0.5, 0.3};
    //    const Double eta[nrEtas] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0};
		// initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    SBM_Vector tmpSol(BaseMatrix::DOUBLE);
    algsys_->GetSolutionVal( tmpSol );

    for( UInt i=0; i<nrEtas; i++) {
      //       if (i>0)
      //         actSol=solOld;
      actSol.Add( 1.0, solOld,
                  eta[i], solIncrement );
      //actSol = solOld + solIncrement * eta[i];
      //      actSol -= solOld;


      //store new solution
      tmpSol = actSol;
      solVec_ = tmpSol;

      // Recalculate residual, f-Cu-Mu-K*u
      algsys_->InitRHS();

      assemble_->AssembleLinRHS();

      // assemble!
      assemble_->AssembleMatrices();

      // account for Dirichlet BCs
      PDE_.SetBCs();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      algsys_->BuildInDirichlet();

      TS_alg_->UpdateRHS(actSol);
      // substract K^* u^k from RHS

      algsys_->RemoveIDBCInfoFromMatrix();

      TS_alg_->SubstractStiffnessFromRHS(actSol);


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2(); // L2Norm of  (f-Ku )

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

    actSol.Add( 1.0, solOld, etaOpt, solIncrement );

    return residualL2NormOpt;
    */
    return 0.0;
  }



  // read nonlinear parameters from xml file
  void StdSolveStep::ReadNonLinData() {

    // Get ParamNode of pde
    PtrParamNode nonLinNode;
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
      StdVector<PtrParamNode> pairCouplings;
      PtrParamNode cplList =
        param->GetByVal("sequenceStep",std::string("index"), actMsStep)
        ->Get("couplingList", ParamNode::PASS);
      if( cplList ) {
        StdVector<PtrParamNode> pairCouplings =
          cplList->Get("direct")->GetChildren();
        // look for each direct coupling, if it has a "nonLin"-node
        for( UInt iCpl = 0; iCpl < pairCouplings.GetSize(); iCpl++ ) {
          nonLinNode = pairCouplings[iCpl]->Get("nonLinear", ParamNode::PASS );
          if( nonLinNode ) {
            break;
          }
        }
      }
      if( !nonLinNode ) {
        // in this case we iterate over all single PDEs and try to
        // find the related entry
        ParamNodeList pdes =
          param->GetByVal("sequenceStep",std::string("index"), actMsStep)
          ->Get("pdeList")->GetChildren();
        for (UInt iPde = 0; iPde < pdes.GetSize(); iPde++ ) {
          if( boost::find_first(pdeName, pdes[iPde]->GetName() ) ) {
            nonLinNode = pdes[iPde]->Get("nonLinear", ParamNode::PASS );
            if( nonLinNode ) {
              break;
            }
          }
        }
      }
    } else {
      nonLinNode = PDE_.GetParamNode()->Get("nonLinear", ParamNode::PASS );
    }

    // Check, if any nonlinear node was found
    if( !nonLinNode ) {
      WARN("Taking default parameters for nonlinear data" );
    }

    // Read data, if "nonLinear" element was found
    if( nonLinNode ) {

      // solution method
      nonLinNode->GetValue( "method", nonLinMethod_, ParamNode::PASS );

      // perform logging?
      nonLinNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );

      // type of line search
      nonLinNode->Get( "lineSearch")->GetValue( "type", lineSearch_ );

      // incremental stopping criterion
      nonLinNode->GetValue( "incStopCrit", incStopCrit_ );

      // residual stopping criterion
      nonLinNode->GetValue( "resStopCrit", residualStopCrit_ );

      // maximal number of NL-iterations
      nonLinNode->GetValue( "maxNumIters", nonLinMaxIter_ );
    }

  }

  void StdSolveStep::ReInit(){
    EXCEPTION("Adjust to new implementation" );
//    dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector()).Init();
//    PDE_.getTimeStepping()->ReInit();
  }
  
  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName, const UInt iterationCounter,
          const Double residualErr, const Double incrementalErr, double etaLineSearch)
  {
    PtrParamNode iter = info->Get("PDE")->Get(pdeName)->Get("nonlinear_iteration", ParamNode::APPEND); 
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
  }
  
} // end of namespace
