#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

#include "StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/EvalIntegrals/BiotSavart.hh"
#include "Utils/Timer.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField {

// declare logging stream
DECLARE_LOG(stdsolvestep)
DEFINE_LOG(stdsolvestep, "stdsolvestep")

  StdSolveStep::StdSolveStep(StdPDE & apde)
    :BaseSolveStep(),
     PDE_(apde)
  {

    pdename_      = PDE_.GetName();
    isaxi_        = PDE_.GetIsaxi();
    subdoms_      = PDE_.GetRegions();
    materialData_ = PDE_.GetMaterialData();
    ptgrid_       = PDE_.GetGrid();
    algsys_       = PDE_.GetAlgSys();
    assemble_     = PDE_.GetAssemble();
    solStrat_     = algsys_->GetSolStrategy();

    results_      = PDE_.GetResultInfos();
    couplingIter_ = 0;
    solutionLimit_ = NO_SOLUTION_TYPE;

    // copy FE functions of PDE
    feFunctions_ = PDE_.GetFeFunctions();
    rhsFeFunctions_ = PDE_.GetRhsFeFunctions();
    
    // Copy vectors FE functions in SBM-vector for communication
    // with OLAS and time stepping
    solVec_.SetSize( feFunctions_.size() );
    rhsVec_.SetSize( feFunctions_.size() );
    


    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator it;
    it = feFunctions_.begin();
    // UInt pos = 0;
    for( ; it != feFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
      solVec_.SetSubVector(ptFct->GetSingleVector(), id);
    }
    //pos = 0;
    it = rhsFeFunctions_.begin();
    for( ; it != rhsFeFunctions_.end(); ++it ){
      shared_ptr<BaseFeFunction> & ptFct = it->second;
      FeFctIdType id = ptFct->GetFctId();
      rhsVec_.SetSubVector(ptFct->GetSingleVector(), id);
    }
    // Make sure to have both vectors as "weak" vectors,
    // as the feFunctions themselves are responsible for
    // creation and destruction.
    solVec_.SetOwnership(false);
    rhsVec_.SetOwnership(false);
    
    // nonlinear parameters
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;
    nonLinMaxIter_ = 10;
    minValidValue_ = -std::numeric_limits<double>::max();// = DBL_MAX;
    maxValidValue_ = std::numeric_limits<double>::max();// = DBL_MAX;

    nonLin_                 = PDE_.IsNonLin();
    nonLinMaterial_         = PDE_.IsNonLinMaterial();
    nonLinTotalFormulation_ = PDE_.IsTotalNonLinFormulation();
    isHyst_                 = PDE_.IsHysteresis();
    regionNonLinTypes_      = PDE_.GetNonLinRegionTypes();

    startStep_ = 1;
    
    // set entry type of SBM vector
    //oldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    //tmpOldRhsLinVal_ = SBM_Vector(BaseMatrix::DOUBLE);
    //DeltaRhsLinVal_ =  SBM_Vector(BaseMatrix::DOUBLE);
    RhsLinVal_ = rhsVec_;

    // In the end, read nonlinear data from xml-file
    if( nonLin_ || nonLinMaterial_ ) {
      ReadNonLinData();
    }
    
    //logFile_.open("nonlin.txt");

    mHandle_ = PDE_.GetDomain()->GetMathParser()->GetNewHandle();
    mParser_ = PDE_.GetDomain()->GetMathParser();
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

  }


  void StdSolveStep::SolveStepStatic() {

    if (nonLin_) {
      StepStaticNonLin();
    }
    else {
      StepStaticLin();
    }
  }


  void StdSolveStep::StepStaticLin() {

    assemble_->AssembleMatrices();

    // The RHS-sources and boundary conditions
    // have to be reassembled each time
    assemble_->AssembleLinRHS();
    //Set special RHS Values
    PDE_.SetRhsValues();

    PDE_.SetBCs();

    // store rhs vector back to PDE
    algsys_->GetRHSVal(rhsVec_);

    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( NO_FCT_ID, 
                                         matrix_factor_[NO_FCT_ID] );
    }

    // Incorporate Boundary conitions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();

      algsys_->SetupSolver();
    }

    // Solve problem
    algsys_->Solve();

    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }


  void StdSolveStep::StepStaticNonLin() {

    bool performOneMoreStep;
    bool isNewton = false;

    SBM_Vector solInc(BaseMatrix::DOUBLE);

    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;

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

      // set the boundary conditions
      PDE_.SetBCs();

      //perform the load-steps
      Double loadFactor = 1.0;
      PDE_.GetInfoNode()->Get("PDE")->Get(pdename_)->
          Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      algsys_->InitRHS();
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // set iteration counter
      UInt iterationCounter=0;

      // =================================
      //  Inner nonlinear loop 
      // =================================
      do {
        iterationCounter++;

        if ( lineSearch_ != "none" || iterationCounter == 1) {
          // set linear part of RHS
          algsys_->InitRHS(RhsLinVal_);

          // setup the matrices
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(SYSTEM,solVec_,true);
          solVec_.ScalarMult(-1.0);
        }


        // assemble Newton bilinear forms
        isNewton = true;
        assemble_->AssembleMatrices(isNewton);

        //compute effective matrix
        algsys_->ConstructEffectiveMatrix( NO_FCT_ID,
                                           matrix_factor_[NO_FCT_ID] );

        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();

        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;

        algsys_->Solve(setIDBC);

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solInc, setIDBC );


        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;
        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration

          actSol.Add(1.0, solInc);
          // store the new solution
          solVec_ = actSol;

          //=================compute residual norm
          // set linear part of RHS
          algsys_->InitRHS(RhsLinVal_);

          // setup the matrices
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(SYSTEM,solVec_,true);
          solVec_.ScalarMult(-1.0);

          //get RHS vector
          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error =======================================
          residualL2Norm = actRHS.NormL2();
          //std::cout << "ResAbsolut: " << residualL2Norm << std::endl;
        }
        else {
          // do line search
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);

          // store the new solution
          solVec_ = actSol;
        }

        // calculation relative residual error ====================================
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

        // output of norms and data
        if ( nonLinLogging_ == true ) {
          //UInt actStep = PDE_.GetSolveStep()->GetActStep();

        if (PDE_.IsIterCoupled()) {
          WriteNonLinIterToInfoXML(pdename_, couplingIter_, iLevel+1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        } else {
            WriteNonLinIterToInfoXML(pdename_, iLevel+1, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        }

        // write norm to file
        logFile_ <<  iterationCounter << "\t"
            << residualErr << "\t"
            << incrementalErr << "\t"
            << etaLineSearch << std::endl;
        }

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
            (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

        
        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' in step no '" << iLevel+1 
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

      // stop timer
      timer->Stop();
    } // loop over levels
  }



  // ======================================================
  // Solve Step Transient SECTION
  // ======================================================

  void StdSolveStep::InitTimeStepping(){
    //also initialize vectors for the time stepping scheme

    stageRHS_.Resize(feFunctions_.size());

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    UInt rhsSize = 0;

    //reserve memory for the rhs
    for(fncIt = feFunctions_.begin(); fncIt != feFunctions_.end();++fncIt ){
      rhsSize = fncIt->second->GetSingleVector()->GetSize();
      FeFctIdType id = fncIt->second->GetFctId();
      stageRHS_.SetSubVector(new Vector<Double>(),id);
      stageRHS_.GetPointer(id)->Resize(rhsSize);
    }
    stageRHS_.Init();

  }

  void StdSolveStep::PreStepTrans() {
    // Update moving ncInterfaces as needed
    ptgrid_->MoveNcInterfaces();
    
    // due to coupling-pdes, the RHS has to be initialized BEFORE
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();
    PDE_.FinilizeBeforTimeStep();
  }


  void StdSolveStep::SolveStepTrans() {
	if ( isHyst_ ) {//&& PDE_.IsHysteresis_Fixpoint() == false ) {
		StepTransNonLinHysteresis();
	}
	//currently not supported
//	else if ( nonLin_ && nonLinMaterial_ ) {
//      StepTransNonLinMaterial();
//    }
    // do a nonlinear time step
    else if (nonLin_ ){//|| PDE_.IsHysteresis_Fixpoint() == true){
      if ( nonLinTotalFormulation_ ){ //|| PDE_.IsHysteresis_Fixpoint() == true){
        //std::cout << "TOTAL" << std::endl;
        StepTransNonLinTotal();
      }
      else
        StepTransNonLin();
    }
    // do a linear time step
    else {
      StepTransLin();
    }
  }


  void StdSolveStep::StepTransLin()
  {
    //TODO: add consistency check here
    //basically loop over all functions and check if the solution order is the same...

    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();

    std::map<FEMatrixType,Integer>::iterator matIt;

    bool effectiveMatrixUpdated = false;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 ); 
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    for(UInt i=0;i<numStages;i++){
      effectiveMatrixUpdated = false;
      rhsVec_.Init();
      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
        FeFctIdType fctId = fncIt->second->GetFctId();
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),fctId);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);


      algsys_->InitRHS();

      //account for RHS
      assemble_->AssembleLinRHS();
      //Set special RHS Values
      PDE_.SetRhsValues();

      // store rhs vector back to PDE 
      algsys_->GetRHSVal(rhsVec_);
      

      // if we want to use static condensation we have to perform the timestepping on element level
      if(algsys_->UseStaticCondensation()){

        matrix_factor_.clear();
        // get timeintegration factors and store them in a map
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
                                  ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
        }
        // assemble all matrix parts (needed to update rhs) and also calculate complete
        // system matrix
        assemble_->AssembleMatrices_CondTrans(false,i,matrix_factor_);

        if(assemble_->IsMatrixUpdated()){

          //std::cout << "in SolveStepLin: new matrices computed" << std::endl;

          // the system matrix was already created so do not init it!
          //        algsys_->InitMatrix(SYSTEM);

          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();

            // we need to call this function
            // but as we already have summed up all entries on element level we just have to give
            // this function a matrix_factor_ list which contains only 0 else
            std::map<FEMatrixType,Double> zero_factors;
            zero_factors[MASS] = 0.0;
            zero_factors[MASS_UPDATE] = 0.0;
            zero_factors[STIFFNESS] = 0.0;
            zero_factors[STIFFNESS_UPDATE] = 0.0;
            zero_factors[DAMPING] = 0.0;
            zero_factors[DAMPING_UPDATE] = 0.0;

            algsys_->ConstructEffectiveMatrix(fctId, zero_factors); //matrix_factor_[fctId]);
          }
          //        matrix_factor_[fctId]
          effectiveMatrixUpdated = true;
        }

      } else {


        assemble_->AssembleMatrices();
        if(assemble_->IsMatrixUpdated()){
          //std::cout << "in SolveStepLin: new matrices computed" << std::endl;

          // set system matrix to zero initially, as ConstructEffectiveMatrix only
          // sums up the contributions
          algsys_->InitMatrix(SYSTEM);
          matrix_factor_.clear();
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()
                  ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
            algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
          }
          effectiveMatrixUpdated = true;
        }
      }

      //now compute the effective right hand side
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        
      
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(fctId));
          }
        algsys_->UpdateRHS(matIt->first,stageRHS_,effectiveMatrixUpdated);
      }

      // set boundary conditions
      PDE_.SetBCs();
      algsys_->BuildInDirichlet();

      if( effectiveMatrixUpdated ){
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
      }

      algsys_->Solve();

      algsys_->GetSolutionVal(stageSol);
    }

    //update stage
    
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end(); ++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep(  );
    }

  }


  void StdSolveStep::StepTransNonLin() {

    //std::cout << "in StepTransNonLin" << std::endl;
    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransNonLin";
    bool performOneMoreStep;
    bool isNewton = false;

    SBM_Vector solInc(BaseMatrix::DOUBLE);

    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;

    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;

    UInt pos = 0;
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 ); 
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }

    for(UInt i=0;i<numStages;i++){
      //do initialization 
      rhsVec_.Init();
      LOG_DBG(stdsolvestep) << "StepTransNonLin: Stage: " << i ;

      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);


      //initialize solution vector for each stage
      if ( i > 0 )
        actSol = stageSol;
      else{
        //special case of incremental non-linearity, we set the stage vector to the solution vector
        stageSol = actSol;
      }

      solVec_  = actSol;

      // setup right hand side
      Double loadFactor = 1.0;
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // set iteration counter
      UInt iterationCounter=0;

      do {
        iterationCounter++;

        if ( lineSearch_ != "none" || iterationCounter == 1) {

          // set linear part of RHS
          algsys_->InitRHS(RhsLinVal_);

          // setup the matrices
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

          //now update RHS according to time stepping
          for(matIt = matrices.begin();matIt != matrices.end();matIt++){

            if(matIt->second < 0)
              continue;
            for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
            }
            algsys_->UpdateRHS(matIt->first,stageRHS_,true);
          }
          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(STIFFNESS,solVec_,true); // we also or only need the updated version
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
          solVec_.ScalarMult(-1.0);
        }

        //now assemble the Newton bilinear forms
        isNewton = true;
        assemble_->AssembleMatrices(isNewton);

        matrix_factor_.clear();
        
        // set system matrix to zero initially, as ConstructEffectiveMatrix only
        // sums up the contributions
        algsys_->InitMatrix(SYSTEM);
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
            ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
        }
        // setup the matrices to compute correct error norms
        
        PDE_.SetBCs();
        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();

        // just set inh. Dirichlet BCs for the first iteration
        bool setIDBC = false;
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          setIDBC = true;

        algsys_->Solve(setIDBC);
        algsys_->GetSolutionVal(solInc, setIDBC );

        //std::cout << solInc.ToString(1,',') << std::endl;

        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;

        //necessary due to inh. Dirichlet BCs!!
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          stageSol.Init();

        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration

          stageSol.Add(1.0, solInc);
          solVec_  = stageSol;

          // setup the matrices with new solution
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

          // set linear part of RHS
          algsys_->InitRHS(RhsLinVal_);
          //assemble_->AssembleNonLinRHS();

          //now update RHS according to time stepping
          for(matIt = matrices.begin();matIt != matrices.end();matIt++){
            //std::cout << "Matrix: " << matIt->first << std::endl;
            if(matIt->second < 0)
              continue;
            for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
            }
            algsys_->UpdateRHS(matIt->first,stageRHS_,true);
          }

          //substract from RHS the term K*sol
          solVec_.ScalarMult(-1.0);
          algsys_->UpdateRHS(STIFFNESS,solVec_,true);
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
          solVec_.ScalarMult(-1.0);

          //get RHS vector
          SBM_Vector actRHS(BaseMatrix::DOUBLE);
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error =======================================
          residualL2Norm = actRHS.NormL2();
        }
        else {
           //solVec_  = stageSol;
          residualL2Norm = LineSearch(solInc, stageSol, etaLineSearch,true);
          solVec_  = stageSol;
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
        Double actSolL2Norm  = stageSol.NormL2();

        if ( actSolL2Norm )
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else {
          incrementalErr = solIncrL2Norm;
          //WARN("Zero solution vector!! ");
        }

        // output of norms and data
        if ( nonLinLogging_ == true ) {
          // get current step 
          UInt actStep = PDE_.GetSolveStep()->GetActStep();

          if (PDE_.IsIterCoupled()) {
            WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
      } else {
            WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
      }
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
              << residualErr << "\t"
              << incrementalErr << "\t"
              << etaLineSearch << std::endl;
        }

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
          (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } //stages
      
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
    limitFeFctIt = feFunctions_.find(solutionLimit_);
    if (limitFeFctIt != feFunctions_.end() ) {
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
          for (UInt j=0; j < dsubVec.GetSize(); j++) {
            if (dsubVec[j] >= maxValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ << 
                      "' is larger than the allowed maximum limit set in the XML: "
                      << maxValidValue_); 
            }
            if (dsubVec[j] <= minValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ << 
                      "' is smaller than the allowed minimum limit set in the XML: "
                      << minValidValue_); 
            }
          }
        }
      }
    }

    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
  }


  void StdSolveStep::StepTransNonLinTotal() {

    bool performOneMoreStep;
    bool isNewton;
    Double incrementalErr;

    SBM_Vector solNew(BaseMatrix::DOUBLE);
    SBM_Vector diffSol(BaseMatrix::DOUBLE);

    //get actual solution
    SBM_Vector  actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;

    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;

    UInt pos = 0;

    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }

    for(UInt i=0;i<numStages;i++){
      //do initialization
      rhsVec_.Init();

      //we obtain a reference to the stage vectors of the scheme
      SBM_Vector stageSol;
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);

      //initialize solution vector for each stage
      if ( i > 0 )
        actSol = stageSol;
      else{
        //special case of incremental non-linearity, we set the stage vector to the solution vector
        stageSol = actSol;
      }

      solVec_  = actSol;

      // setup right hand side
      Double loadFactor = 1.0;
      incrementalErr = SetLinRHS(loadFactor);

      // set iteration counter
      UInt iterationCounter=0;

      do {
        iterationCounter++;

        // do matrices: Newton is not working for total formulation!!
        isNewton = false;
        assemble_->AssembleMatrices(isNewton);

        // set RHS
        algsys_->InitRHS(RhsLinVal_);

        /*
         * set nonlinrhs via assmble and not via setlinrhs(...,true)
         */
        assemble_->AssembleNonLinRHS();

        //now update RHS according to time stepping
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }

        // set system matrix to zero initially, as ConstructEffectiveMatrix only
        // sums up the contributions
        matrix_factor_.clear();
        algsys_->InitMatrix(SYSTEM);
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()
            ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
        }

        PDE_.SetBCs();
        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();

        //always set inhomog. Dirichlet BCs
        bool setIDBC = true;

        algsys_->Solve(setIDBC);
        algsys_->GetSolutionVal(solNew, setIDBC );

        // calculate incremental error ========================================
        diffSol = solNew;
        diffSol.Add( -1.0, actSol);
        Double solIncrL2Norm = diffSol.NormL2();
        Double solNewL2Norm = solNew.NormL2();

        if (solNewL2Norm > 1)
          incrementalErr = solIncrL2Norm / solNewL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        //just dummy things
        Double etaLineSearch = 1.0;
        Double residualErr = incrementalErr;

        // output of norms and data
        if ( nonLinLogging_ == true ) {
          WriteNonLinIterToInfoXML(pdename_, 1,iterationCounter, residualErr, incrementalErr, etaLineSearch);
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
              << residualErr << "\t"
              << incrementalErr << "\t"
              << etaLineSearch << std::endl;
        }

        stageSol = solNew;
        solVec_  = stageSol;

        //store new solution
        actSol = solNew;

        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep =
          (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                  << "' at iteration '" << iterationCounter 
                  << "'.\n ==> incremental error: " << incrementalErr
                  << "\n ==> residual error: " << residualErr);
        }

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } //stages

    //update stage
    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }

  }

  
  void StdSolveStep::StepTransNonLinHysteresis() {

   bool performOneMoreStep;
   bool isNewton;
   Double incrementalErr;

   SBM_Vector newSol(BaseMatrix::DOUBLE);
   SBM_Vector oldSol(BaseMatrix::DOUBLE);
   SBM_Vector diffSol(BaseMatrix::DOUBLE);

   //get solution of previous time step
   SBM_Vector  prevSol(BaseMatrix::DOUBLE);
   prevSol = solVec_;
   oldSol = solVec_;

   std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
   std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
   std::map<FEMatrixType,Integer>::iterator matIt;

   UInt pos = 0;
   bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
   for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
     fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
   }

   //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 )
  	   EXCEPTION("StepTransNonLinHysteresis: just one stage time-stepping allowed");

   for(UInt i=0;i<numStages;i++){
     //do initialization
     rhsVec_.Init();

     //we obtain a reference to the stage vectors of the scheme
     SBM_Vector stageSol;
     stageSol.Resize(feFunctions_.size());
     for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
       stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
       fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
     }
     stageSol.SetOwnership(false);

     //initialize solution vector for each stage
     if ( i > 0 )
       oldSol = stageSol;
     else{
       //special case of incremental non-linearity, we set the stage vector to the solution vector
       stageSol = oldSol;
     }

     solVec_  = oldSol;

     // setup right hand side
     Double loadFactor = 1.0;

     incrementalErr = SetLinRHS(loadFactor);

     // set iteration counter
     UInt iterationCounter=0;

     do {
       iterationCounter++;

       //std::cout << "Current iteration: " << iterationCounter << std::endl;

       // reset rhs
     //  RhsLinVal_.Init();
     //  algsys_->InitRHS(RhsLinVal_);
//
//       if(PDE_.IsHysteresis_Fixpoint() == true){
//
//         // nach dem löschen muss er eigentlich sowohl die linearen als auch die nichtlinearen teile nochmals hinzufügen
//         // die unterscheidung zwischen linear und nichtlinear ist also obsolet
//
//         RhsLinVal_.Init();
//         algsys_->InitRHS(RhsLinVal_);
//
//         //incrementalErr = SetLinRHS(loadFactor,true);
//         /*
//          * Fixpoint iteration used to work, however, P was just postprocessed upon E
//          * as we RhsLinVal_.Init() followed by SetLinRHS( ,true) we just had nonlinear terms on the rhs
//          * although P depends on E it was not marked as solDependent and such we solved div(eps0*E) = 0
//          * this worked of course
//          * Doing a real updating on the rhs will however not work
//          */
//         //incrementalErr = SetLinRHS(loadFactor,false);
//       } else {
//         //incrementalErr = SetLinRHS(loadFactor,false);
//       }

       // do matrices: Newton is not working for total formulation!!
       isNewton = false;
       assemble_->AssembleMatrices(isNewton);

       // set RHS
       algsys_->InitRHS(RhsLinVal_);

       //now update RHS according to time stepping
       for(matIt = matrices.begin();matIt != matrices.end();matIt++){
         if(matIt->second < 0)
           continue;
         for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
           fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
         }
         algsys_->UpdateRHS(matIt->first,stageRHS_,true);
       }

       if(PDE_.IsHysteresis_Fixpoint() == false){
         //due to incremental material formulation
         //std::cout << "RHS updated with old solution" << std::endl;
         algsys_->UpdateRHS(STIFFNESS,prevSol,true);
       }
       // set system matrix to zero initially, as ConstructEffectiveMatrix only
       // sums up the contributions
       matrix_factor_.clear();
       algsys_->InitMatrix(SYSTEM);

       for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
         FeFctIdType fctId = fncIt->second->GetFctId();
         fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
         algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
       }

       PDE_.SetBCs();
       algsys_->BuildInDirichlet();
       algsys_->SetupPrecond();
       algsys_->SetupSolver();

       //always set inhomog. Dirichlet BCs
       bool setIDBC = true;

       algsys_->Solve(setIDBC);
       algsys_->GetSolutionVal(newSol, setIDBC );

       // calculate incremental error ========================================
       diffSol = newSol;
       diffSol.Add( -1.0, oldSol);
       Double solIncrL2Norm = diffSol.NormL2();
       Double solNewL2Norm = newSol.NormL2();

       if (solNewL2Norm > 1)
         incrementalErr = solIncrL2Norm / solNewL2Norm;
       else
         incrementalErr = solIncrL2Norm;

       //std::cout << "rel error: " << incrementalErr << std::endl;
       //std::cout << "abs error: " << solIncrL2Norm << std::endl;

       //just dummy things
       Double etaLineSearch = 1.0;

       // Why is residualErr = incrementalErr?
       Double residualErr = incrementalErr;

       // output of norms and data
       if ( nonLinLogging_ == true ) {
         WriteNonLinIterToInfoXML(pdename_, 1,iterationCounter, residualErr, incrementalErr, etaLineSearch);
         // write norm to file
         logFile_ <<  iterationCounter << "\t"
             << residualErr << "\t"
             << incrementalErr << "\t"
             << etaLineSearch << std::endl;
       }

       //use relaxated form
       stageSol = newSol;
       //Double relaxfac = 0.1;
       //stageSol.Add( -relaxfac, newSol);
       //stageSol.Add( relaxfac, oldSol);
       //stageSol = newSol;
       solVec_  = stageSol;

       //store new solution
       oldSol = stageSol; //newSol;

       // boolean variable, holds condition if another iteration step is necessary
       performOneMoreStep =
         (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

      // std::cout << "IncError: " << incrementalErr << std::endl;
      // std::cout << "ResError: " << residualErr << std::endl;

       if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
         EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
                 << "' at iteration '" << iterationCounter
                 << "'.\n ==> incremental error: " << incrementalErr
                 << "\n ==> residual error: " << residualErr);
       }

     } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

   } //stages

   //update stage
   for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
     fncIt->second->GetTimeScheme()->FinishStep();
   }

   PDE_.FinalizeAfterTimeStep();
  }



  void StdSolveStep::StepTransNonLinMaterial() {

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
//    StepTransLin;
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
//      algsys_->SetupPrecond;
//      algsys_->SetupSolver;
//      algsys_->Solve;
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
    algsys_->InitRHS();
  }


  void StdSolveStep::SolveStepHarmonic() {
    if ( nonLin_ ) {
      StepHarmonicNonLin();
    }
    else {
      StepHarmonicLin();
    }
  }


  void StdSolveStep::StepHarmonicLin() {

    //JUST A HACK!!!!
    //matrix_factor_Complex_[NO_FCT_ID][STIFFNESS] = Complex(1.0,0);
    //matrix_factor_Complex_[NO_FCT_ID][DAMPING] = CompSetLinRHSlex(0.0,actFreq_*2*M_PI);
    //matrix_factor_Complex_[NO_FCT_ID][MASS] = Complex(-1.0 * actFreq_*actFreq_*4*M_PI*M_PI,0);

    //matrix_factor_Complex_[NO_FCT_ID][STIFFNESS] = Complex(1.0,0);
    //matrix_factor_Complex_[NO_FCT_ID][DAMPING] = Complex(1.0,0.0);
    //matrix_factor_Complex_[NO_FCT_ID][MASS] = Complex(1.0,0.0);


    //this has to be done each frequency!
    assemble_->AssembleLinRHS();

    //Set special RHS Values
    PDE_.SetRhsValues();

    assemble_->AssembleMatrices( );
    PDE_.SetBCs();

    // store rhs vector back to PDE
    algsys_->GetRHSVal( rhsVec_ );

    // Where should we get the matrix factors from in a harmonic case?
    // In my opinion this method
    //if( assemble_->IsMatrixUpdated() ) {
    std::map<FEMatrixType,Double> empty;
    algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty );

    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }

    algsys_->Solve();
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
    algsys_->SetupEigenSolver( numFreq, shift, false, false);

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies(frequencies, errBounds);

    return frequencies.GetSize();
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift, bool bloch) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver  - we cannot be quadratic and bloch concurrently!
   algsys_->SetupEigenSolver( numFreq, shift, !bloch, bloch);

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds );
    algsys_->ExportLinSys(false, false, true); // setup

    return frequencies.GetSize();
  }

  void StdSolveStep::GetEigenMode( UInt numMode ) {


    algsys_->GetEigenMode( numMode );

    // Get the solution and store it
    algsys_->GetSolutionVal(solVec_);
  }




  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor, bool nonlin)
  {

    Double RhsLinL2Norm;

    // to incorporate loads+
    if(nonlin){
      assemble_->AssembleNonLinRHS();
    } else {
      assemble_->AssembleLinRHS();
    }

    //Set special RHS Values
    PDE_.SetRhsValues();

    // Stores rhs vector into extForces and returns that L2-norm
    algsys_->GetRHSVal( RhsLinVal_ );

    RhsLinVal_.ScalarMult(loadFactor);

    //std::cout << "RHS: " << RhsLinVal_.ToString() << std::endl;

    RhsLinL2Norm = RhsLinVal_.NormL2();

    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
        // Note: there are PDEs, such as elecconduction, which always have rhs=0. Those should not emit this warning. SE.
      if (pdename_ != "elecConduction")
        WARN("Zero external force vector!! ");
    }

    return RhsLinL2Norm;
  }

  UInt StdSolveStep::SetDeltaLinRHS()
  {

    // to incorporate loads
    assemble_->AssembleLinRHS(); 
    //Set special RHS Values
    PDE_.SetRhsValues();
    
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
                                  Double& etaLineSearch, bool trans)  {

    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {0.1, 0.25, 0.5, 1.0}; //, 0.5, 0.25, 0.125, 0.1};
        // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      actSol.Add( 1.0, solOld, eta[i], solIncrement);

      //store new solution
      solVec_ = actSol;
      
      // set RHS: linear part
      algsys_->InitRHS(RhsLinVal_ );
      // and nonlinpart if any
      assemble_->AssembleNonLinRHS();

      // setup the matrices
      bool isNewton = false;
      assemble_->AssembleMatrices(isNewton);

      if( trans ) {
        //now update RHS according to time stepping
        std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
        std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
        std::map<FEMatrixType,Integer>::iterator matIt;
        UInt pos = 0;
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(0,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }

        //substract from RHS the term K*sol
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(STIFFNESS,solVec_,true);
        algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
        solVec_.ScalarMult(-1.0);
      }
      else {
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(SYSTEM,solVec_,true);
        solVec_.ScalarMult(-1.0);
      }


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = actRHS.NormL2();

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;
    
    // Set new solution
    actSol.Add( 1.0, solOld, etaOpt, solIncrement );

    return residualL2NormOpt;
  }


  Double StdSolveStep::LineSearchMag(SBM_Vector& solIncrement, SBM_Vector& actSol,
                                          Double& etaLineSearch, bool trans)
  {

    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.1};
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
        //EXCEPTION("Line Search for nonlinear transient problems yet verified")
//        assemble_->AssembleNonLinRHS();
//        TS_alg_->UpdateRHS(actSol);
        assemble_->AssembleNonLinRHS();
        //PDE_.SetRhsValues();
        //now update RHS according to time stepping
        std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
        std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
        std::map<FEMatrixType,Integer>::iterator matIt;
        UInt pos = 0;
        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
          if(matIt->second < 0)
            continue;
          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
            fncIt->second->GetTimeScheme()->ComputeStageRHS(0,matIt->second,stageRHS_.GetPointer(pos));
          }
          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
        }
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
    actSol.Add(1.0, solOld, etaOpt, solIncrement );

    return residualL2NormOpt;
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
    PtrParamNode nonLinNode = solStrat_->GetNonLinNode();

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
      if( nonLinNode->Has("lineSearch") ) {
        nonLinNode->Get( "lineSearch")
            ->GetValue( "type", lineSearch_,ParamNode::PASS );
      }

      // incremental stopping criterion
      nonLinNode->GetValue( "incStopCrit", incStopCrit_, ParamNode::PASS );

      // residual stopping criterion
      nonLinNode->GetValue( "resStopCrit", residualStopCrit_, 
                            ParamNode::PASS );

      // maximal number of NL-iterations
      nonLinNode->GetValue( "maxNumIters", nonLinMaxIter_, ParamNode::PASS );

      // abort if max number of iterations is reached?
      nonLinNode->GetValue("abortOnMaxNumIters",abortOnMaxIter_,ParamNode::INSERT);

      LOG_TRACE(stdsolvestep) << "Nonlinear convergence criteria were read:";
      LOG_DBG3(stdsolvestep) << "\tincremental Stopping Criterion: " << incStopCrit_;
      LOG_DBG3(stdsolvestep) << "\tresidual Stopping Criterion: " << residualStopCrit_;
      LOG_DBG3(stdsolvestep) << "\tmaxNumIters: " << nonLinMaxIter_;
      if( nonLinNode->Has("stopOnLimit") ) {
        std::string solutionString;
        // quantity
        solutionString = nonLinNode->Get( "stopOnLimit")
            ->Get( "quantity" )->As<std::string>();
        solutionLimit_ = SolutionTypeEnum.Parse(solutionString);
        if (feFunctions_.find(solutionLimit_) == feFunctions_.end() ) 
          EXCEPTION("ERROR: Solution type '" << solutionString << "' is not part of PDE '" << pdename_ << 
                      "' and cannot serve as stopping limit criterion");
        // minimum value
        nonLinNode->Get( "stopOnLimit")
            ->GetValue( "min", minValidValue_, ParamNode::PASS );
        // maximum value
        nonLinNode->Get( "stopOnLimit")
            ->GetValue( "max", maxValidValue_, ParamNode::PASS );
        // region
        std::string solutionRegion;
        solutionRegion = nonLinNode->Get( "stopOnLimit")
            ->Get( "region" )->As<std::string>();
        solutionLimitReg_ = ptgrid_->GetRegion().Parse(solutionRegion);
        if (solutionLimitReg_ != ALL_REGIONS) {
          if( subdoms_.Find(solutionLimitReg_) == -1 ) 
              EXCEPTION("ERROR: Region '" << solutionRegion <<"' is not part of PDE '" << pdename_ << 
                      "' and cannot serve as stopping limit region");

          // don't know how to implement region checking as we only have the solution as vector (glmvector_[0])
          EXCEPTION("checking limits in a region is not implemented (and I don't know how) (SE). " 
                  << "Just remove the 'region=XXX' tag");
        }
        LOG_DBG3(stdsolvestep) << "\tStop on Limit: " << solutionString;
        LOG_DBG3(stdsolvestep) << "\t\tmin:     " << minValidValue_;
        LOG_DBG3(stdsolvestep) << "\t\tmax:     " << maxValidValue_;
        LOG_DBG3(stdsolvestep) << "\t\tin Region:     " << solutionRegion;
        LOG_DBG3(stdsolvestep) << "\t\tRegionID:     " << solutionLimitReg_;
      }
    }
  }

  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName,
                                              const UInt solStep,
                                              const UInt iterationCounter,
                                              const Double residualErr, 
                                              const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT)
        ->Get("iteration",ParamNode::APPEND);
    iter->Get("pdeName")->SetValue(pdeName);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
  }

  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName,
                                              const UInt coupledIterStep,
                                              const UInt solStep,
                                              const UInt iterationCounter,
                                              const Double residualErr, 
                                              const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT);
    iter = iter->GetByVal("couplingStep", "value", coupledIterStep, ParamNode::INSERT)
        ->Get("iteration",ParamNode::APPEND);
    iter->Get("pdeName")->SetValue(pdeName);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    // SE: include PDE name and time step
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
  }
  
} // end of namespace
