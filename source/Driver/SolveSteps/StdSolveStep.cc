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
//#include "MatVec/SingleVector.hh"


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
    
//    resVec_.SetSize( feFunctions_.size() );
//    nonLinRHS_.SetSize( feFunctions_.size() );

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
    /*
     * also set residual and nonLinRHS
     */
    resVec_ = rhsVec_;
    nonLinRHS_ = rhsVec_;

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
    // Set special RHS Values
    PDE_.SetRhsValues();

    PDE_.SetBCs();

    // store rhs vector back to algsys
    algsys_->GetRHSVal(rhsVec_);

    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( NO_FCT_ID, 
                                         matrix_factor_[NO_FCT_ID] );
    }

    // Check if the AMG-framework is used (if so, we have
    // to gather some geometry information at this point)
    if(algsys_->UseAMG() ){
      // only works for elimination
      if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
      PDE_.SetGeomInfo();
      algsys_->BuildAMGAuxMatrix();
    }


    // Incorporate Boundary conditions and
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
          //add linear right hand side
          algsys_->InitRHS(RhsLinVal_);

          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
            assemble_->AssembleNonLinRHS();
          }

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

        //compute norms (residual and incremenal ones)
        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;
        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration

          actSol.Add(1.0, solInc);
          // store the new solution
          solVec_ = actSol;

          //=================compute residual norm
          algsys_->InitRHS(RhsLinVal_);
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
            assemble_->AssembleNonLinRHS();
          }

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
    if ( isHyst_ ){
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
          //if AMG is used, rebuild auxiliary matrix
          auxSet_ = false;
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

      // Check if the AMG-framework is used (if so, we have
      // to gather some geometry information at this point)
      // needs only be built once, doesn't change over frequency
      if( (algsys_->UseAMG()) && (auxSet_ == false) ){
        // only works for elimination
        if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
        PDE_.SetGeomInfo();
        algsys_->BuildAMGAuxMatrix();
        auxSet_ = true;
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
    /*!
     * Comments added to better understand what's going on. If you find any errors, please correct.
     */
    /*! Solve step for transient, nonlinear simulations using material non-linearities
     *  Works for std fixpoint iteration as well as for Newton iteration.
     *  Instead of computing the total value u_n+1^k+1 at time step n+1 and iteration k+1,
     *  only the difference deltaU^k+1 = u_n+1^k+1 - u_n+1^k is computed.
     *
     *  This leads to the following iterative scheme (only fixpoint case described here)
     *    n = time step counter; k = iteration counter
     *    K(u) = nonlinear system matrix; f = rhs load
     *    Known from last iteration: u_n+1^k
     *    Aim: compute u_n+1^k+1 such that
     *      a) u_n+1^k+1 - u_n+1^k = deltaU^k+1 < tol        -> small incremental error
     *      b) K( u_n+1^k+1 )*u_n+1^k+1 - f     < tol        -> small residual error
     *
     *    Approach:
     *    1. Build up system:
     *      K( u_n+1^k )*u_n+1^k+1 = f
     *    > K( u_n+1^k )*(u_n+1^k + deltaU^k+1) = f
     *    > K( u_n+1^k )*deltaU^k+1 = f - K( u_n+1^k )*u_n+1^k
     *
     *    2. solve for solution increment deltaU^k+1
     *
     *    3. calculate new solution by updating the old one
     *      u_n+1^k+1 = u_n+1^k + deltaU^k+1
     *
     *    4. check for incremental and residual error
     *
     *    5. continue if one of the errors is too large and the maximal number of iterations
     *        has not been reached yet
     *
     *    Three possible implementations (I and II work well at least for the NLtrans magnetic example from
     *      testsuite, III was not tested in that context)
     *      I. u_n+1^0 = u_n
     *         u_n+1^1 = deltaU^1
     *         u_n+1^2 = deltaU^1 + deltaU^2
     *
     *         This is how it is implemented below. During the first iteration, we calculate
     *         deltaU^1 between u_n and u_n+1^1, but we do NOT add deltaU^1 to u_n
     *         to approximate u_n+1^1.
     *         Instead we set u_n+1^1 = deltaU^1.
     *         This leads to (completely) different evaluation points for the material curve and seems
     *         a little bit odd at first. However, it converges to the same solution as II.
     *
     *      II. u_n+1^0 = u_n, but rhs will NOT be updated during first iteration
     *         -> we solve not for deltaU^1 but for the whole vector u_n+1^1 during the first iteration
     *         -> i.e. we perform one "total" step (like it is done in StepTransNonLinTotal())
     *
     *         K( u_n+1^0 )*u_n+1^1 = f
     *         -> u_n+1^1 will thus contain u_n
     *
     *         u_n+1^2 = u_n + deltaU^1 + deltaU^2
     *         -> all further iterations are standard delta-steps
     *         -> compared to I. we do not evaluate the material tensor only with the summed up increments,
     *            but with the increments added to the previous solution
     *         -> tests showed, that I. and II. will converge to the same result, although the calculated steps
     *            to come there are different
     *
     *      III. u_n+1^0 = u_n
     *         u_n+1^1 = u_n + deltaU^1
     *         u_n+1^2 = u_n + deltaU^1 + deltaU^2
     *         ...
     *         -> solution from old time step is the basis for the new solution and will be kept
     *         -> seems to me the most natural way but it has the problem when the system is solved
     *            for the first time during each time step. To incorporate changing inhomogeneous Dirichlet
     *            Boundary Conditions (IDBC) into the result vector, we would have to replace the IDBC values
     *            from the last time step with the ones from the new time step. A simple addition would lead
     *            to wrong results here
     *         -> Added function ClearIDBCFromSolutionVal to set all IDBC nodes to 0
     *            -> call this function with the OLD solution before solving the system
     *            -> the old solution will have 0 at IDBC nodes so that the new solution can simply be added
     *
     *
     *
     *
     */

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
        	//add linear right hand side
        	algsys_->InitRHS(RhsLinVal_);

        	// if the RHS depends on the nonlinearity, we have to re-assemble it
        	if( assemble_->IsRhsSolDependent() ) {
        		assemble_->AssembleNonLinRHS();
        	}

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
        // if setIDBC is true, solInc will contain the inhom. Dirichlet values
        algsys_->GetSolutionVal(solInc, setIDBC );

        Double residualL2Norm = 0.0;
        Double etaLineSearch  = 1.0;

        //necessary due to inh. Dirichlet BCs!!
        if ( iterationCounter == 1 && couplingIter_ == 0 )
          stageSol.Init();

        if ( lineSearch_ == "none" || iterationCounter == 1) {
          //to incooperate the inhomog. Dirichlet BCs we need a full
          //step for the first iteration

          /*
           * note: stageSol was initialized to 0.0, so
           *      during the iterations, stageSol will be the sum of all solultion increments but will NOT
           *      contain the old solution itself
           */
          stageSol.Add(1.0, solInc);

          solVec_  = stageSol;

          //=================compute residual norm
          algsys_->InitRHS(RhsLinVal_);
          // if the RHS depends on the nonlinearity, we have to re-assemble it
          if( assemble_->IsRhsSolDependent()) {
        	  assemble_->AssembleNonLinRHS();
          }

          // setup the matrices with new solution
          isNewton = false;
          assemble_->AssembleMatrices(isNewton);

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
      /*
       * here we finally compute the new solution vector
       *  solution_new = solution_old + stage_solutions
       */
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

  void StdSolveStep::CalcResidualAndConfigSystemForHysteresis(SBM_Vector& oldSolution,SBM_Vector& solIncrement,Double usedEta, UInt stage, UInt callingCnt, UInt evalVersion, bool trans){
    /*
     * GENERAL NOTE 1
     *  this function does the following
     *    1. evaluate residual using newSolution = oldSolution+usedEta*solIncrement
     *      a) overwrite resVec_ with new residual
     *      b) overwrite nonLinRHS_ with current value; note: nonLinRHS_ does NOT
     *          contain hysteresis values
     *      c) overwrite solVec_ with oldSolution + usedEta*solIncrement
     *    2. setup system to be solved for the next iteration
     *
     *  if this function is called during lineSearch, i.e. multiple times with
     *    different values of usedEta, make sure to reset nonLinRHS_, resVec_
     *    and solVec_ after each call!
     *
     * GENERAL NOTE 2
     *  this function is called twice during the first iteration
     *  a) at the beginning of the iteration, before solving the system to
     *      simply setup the equation system to be solved;
     *      usedEta = 0 in that case and oldSolution is the final solution of
     *      the previous time step
     *      -> callingCnt = 0
     *  b) after the solution of the system with usedEta = 1
     *      -> callingCnt = iterationCnt
     *
     *  afterwards this function is called only once during each iteration
     *      -> callingCnt = iterationCnt
     */

    std::cout << "Residual at start of CalcResidual" << std::endl;
    std::cout << resVec_.ToString() << std::endl;
    std::cout << "callingCnt: " << callingCnt << std::endl;



    /*
     * set current solution by overwriting solVec_
     */
    if(usedEta != 0){
      solVec_.Add( 1.0, oldSolution, usedEta, solIncrement);
    }

    /*
     * get empty SBM_Vector for storing temporalRHS
     */
    SBM_Vector tmpRHS(BaseMatrix::DOUBLE);
    algsys_->GetRHSVal( tmpRHS );
    tmpRHS.Init();

    /*
     * assemble nonlinear rhs but without hysteresis (if any)
     * use current/new solution for evaluation
     */
    algsys_->InitRHS(tmpRHS); //rhs = 0

    /*
     * set returnZeroValues to true -> coefFncHyst will return 0 when getVector /getScalar is called
     */
    PDE_.SetFlagInCoefFncHyst("returnZeroValues",true);
    assemble_->AssembleNonLinRHS(); //rhs = f_nonlin without hysteresis

    /*
     * store nonlinear rhs without hysteresis in tmpRHS;
     * at end of function, use tmpRHS to overwrite nonLinearRHS_;
     */
    algsys_->GetRHSVal( tmpRHS );

    /*
     *  To make the following steps easier to follow, consider an Electrostatic
     *  problem of the kind
     *    div( D ) = div( eps0*E + P(E) ) = f_lin + f_nonlin(E)
     *
     *  Nomenclature:
     *    E = field intensity
     *    P(.) = hysteresis operator
     *    P_n+1^k = evaluated value of hysteresis operator
     *    f_lin = linear rhs
     *    f_nonlin(E) = nonlinear rhs EXCLUDING hysteresis
     *    n = time step counter, starting at 1
     *    k = iteration counter during time step, starting at 1
     *    c = calling counter during time step, starting at 0
     *          basically the same as k, only for an additional 0-th entry
     *          as function is called twice during first iteration
     *    final = last iteration during a time step
     *
     *  Starting values at time step n+1
     *    E_n+1^0 = E_n^final
     *    P_n+1^0 = P( E_n^final )
     *    deltaE^k = E_n+1^k - E_n+1^k-1
     *
     *  Approach for k = 0
     *    E_n+1^-1 = E_n^0 = E_n-1^final
     *    deltaE^0 = E_n+1^0 - E_n^0
     *             = E_n^final - E_n-1^final
     *
     *
     *  Updating scheme
     *    E_n+1^k = E_n+1^k-1 + usedEta * deltaE^k
     *
     *  DeltaMaterial-Approach
     *    1. Assume Es = E_n^final solves system
     *      div( eps0*Es + P(Es) ) = f_lin + f_nonlin(Es)
     *
     *    2. Ansatz E_n+1^k+1 = E_n+1^k + deltaE^k+1
     *    3. System to be solved
     *      div( eps0*(E_n+1^k + deltaE^k+1) + P(E_n+1^k) + deltaP^k+1 )
     *        = f_lin + f_nonlin(E_n+1^k+1)
     *
     *      with
     *        deltaP^k+1 = P(E_n+1^k+1) - P(E_n+1^k)
     *
     *    3.0 hysteresis free approach -> only for testing of hysteresis operator
     *      ( Idea: solve standard (non-)linear system and evaluate hysteresis
     *              operator only during calculation of output )
     *
     *       Actual system that can be solved
     *
     *       div( [eps' * deltaE^k+1 )
     *            = f_lin + f_nonlin(E_n+1^k) - div( eps0*E_n+1^k )
     *            = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
     *               + sum_jj_from_1_to_k ( f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1) )
     *               - sum_jj_from_1_to_k ( eps0*deltaE^jj )
     *
     *       with res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
     *       and deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                      - usedEta^jj*eps0*deltaE^jj
     *
     *       div( eps' * deltaE^k+1 )
     *            = res0 + sum_jj_from_1_to_k ( deltaRes^jj )
     *
     *      Residuals
     *           res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
     *           res^k+1 = res0 + sum_jj_from_2_to_k+1 ( deltaRes^jj )
     *
     *           deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                           - usedEta^jj*eps0*deltaE^jj
     *
     *    3.1 deltaP/deltaE approach
     *      -> put P(E_n+1^k) on rhs
     *
     *      -> use deltaMatrix of form
     *            eps' + deltaP^k/deltaE^k (eps' basically is eps0, but during
     *                                    the first iteration a larger value
     *                                    might be reasonable for better convergence)
     *
     *      Actual system that can be solved
     *      ( P(E_n+1^k+1) and f_nonlin(E_n+1^k+1) cannot be evaluated at the
     *        unknown E_n+1^k+1 of course )
     *
     *      Iteration k+1:
     *
     *      div( [eps' + deltaP^k/deltaE^k]*deltaE^k+1 )
     *            = f_lin + f_nonlin(E_n+1^k) - div( eps0*E_n+1^k ) - div( P_n+1^k)
     *            = f_lin + f_nonlin(E_n+1^k) - div( eps0*E_n+1^0 ) - div( P_n+1^0)
     *              - sum_jj_from_1_to_k ( [eps0 + deltaP^jj/deltaE^jj]*deltaE^jj )
     *            = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 ) - div( P_n+1^0)
     *              - sum_jj_from_1_to_k ( [eps0 + deltaP^jj/deltaE^jj]*deltaE^jj )
     *              + sum_jj_from_1_to_k ( f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1) )
     *
     *       with res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 ) - div( P_n+1^0)
     *       and deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *       div( [eps' + deltaP^k/deltaE^k]*deltaE^k+1 )
     *            = res0 + sum_jj_from_1_to_k ( deltaRes^k )
     *
     *      Residuals
     *           res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 ) - div( P_n+1^0)
     *           res^k+1 = res0 + sum_jj_from_2_to_k+1 ( deltaRes^jj )
     *
     *           deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *    3.2 P/deltaE approach
     *      -> during first iteration (k=1) compute deltaMatrix using
     *            eps' + P_n+1^0/deltaE^0
     *
     *      -> during all other iterations, use deltaP/deltaE approach
     *            eps' + deltaP^k/deltaE^k (eps' basically is eps0, but during
     *                                    the first iteration a larger value
     *                                    might be reasonable for better convergence)
     *
     *      Actual system that can be solved
     *      Iteration 1:

     *      div( [eps' + P_n+1^0/deltaE^0]*deltaE^1 )
     *            = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
     *            = res0
     *
     *      Iteration 2:
     *
     *      div( [eps' + deltaP^1/deltaE^1]*deltaE^2 )
     *            = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *            = res1
     *
     *      Iteration k+1>2
     *
     *      div( [eps' + deltaP^k/deltaE^k]*deltaE^k+1 )
     *            = res1 + sum_jj_from_2_to_k ( deltaRes^k )
     *
     *      with res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *      and deltaRes^k = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                      - [eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *      Residuals
     *           res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
     *           res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *           res^k+1 = res1 + sum_jj_from_2_to_k+1 ( deltaRes^jj )
     *
     *           deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *    3.3 P/E approach
     *      -> during first iteration (k=1) compute deltaMatrix using
     *            eps' + P_n+1^0/E_n+1^0
     *
     *      -> during all other iterations, use deltaP/deltaE approach
     *            eps' + deltaP^k/deltaE^k (eps' basically is eps0, but during
     *                                    the first iteration a larger value
     *                                    might be reasonable for better convergence)
     *
     *      Actual system that can be solved
     *      Iteration 1:
     *
     *      div( [eps' + P_n+1^0/E_n+1^0]*E_n+1^1 )
     *            = f_lin + f_nonlin(E_n+1^0)
     *            = res0
     *
     *      Iteration 2:
     *
     *      div( [eps' + deltaP^1/deltaE^1]*deltaE^2 )
     *            = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *            = res1
     *
     *      Iteration k+1>2
     *
     *      div( [eps' + deltaP^k/deltaE^k]*deltaE^k+1 )
     *            = res1 + sum_jj_from_2_to_k ( deltaRes^k )
     *
     *      with res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *      and deltaRes^k = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                      - [eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *      Residuals
     *           res0 = f_lin + f_nonlin(E_n+1^0)
     *           res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
     *           res^k+1 = res1 + sum_jj_from_2_to_k+1 ( deltaRes^jj )
     *
     *           deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
     *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
     *
     *
     */
     /*
      * NOTE:
      *   1. this function is bascially meant for the evaluation of the residual
      *   2. as the evaluation of the residual requires the setup of the equation
      *       system, this function acts as helper to build up the equation
      *       system during prior to the solution step, too
      *   3. during the first iteration, the residual is calculated by a different
      *       equation system than the one to be solved; in that case a second
      *       assemble-step is performed; otherwise, the system can be reused
      *       Example:
      *         version1:
      *         res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 ) - div( P_n+1^0)
      *         -> to subtract div( eps0*E_n+1^0 ) we have to assemble the system
      *             with a material tensor of eps0 only
      *         lhs = div( [eps' + deltaP^0/deltaE^0]*deltaE^1 )
      *         -> system to be solved requires material tensor of eps'+deltaP/deltaE
      *         -> reassemble system after computation of residual
      */

    std::cout << "EvalVersion: " << evalVersion << std::endl;

    /*
     * 1. CALCULATE RESIDUAL
     */
    if(callingCnt == 0){
      /*
       * compute res0
       *
       * 3.0
       *  res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
       * 3.1
       *  res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 ) - div( P_n+1^0 )
       * 3.2
       *  res0 = f_lin + f_nonlin(E_n+1^0) - div( eps0*E_n+1^0 )
       * 3.3
       *  res0 = f_lin + f_nonlin(E_n+1^0)
       *
       * -> all versions have f_lin and f_nonlin in common
       * -> 3.3 does not update with old solution, but it might need system
       *     with eps0 to update with stage solution -> assemble system, too
       * -> 3.1 is the only version that puts div( P_n+1^0 ) on rhs
       */

      algsys_->InitRHS(RhsLinVal_); //f_lin

      if(evalVersion == 1){
        // call to hysteresis function will return P
        // -> assembleNonLinRHS will also add P
        PDE_.SetFlagInCoefFncHyst("returnZeroValues",false);
      } else {
        // call to hysteresis function will return 0
        // -> assembleNonLinRHS will not add P
        PDE_.SetFlagInCoefFncHyst("returnZeroValues",true);
      }
      assemble_->AssembleNonLinRHS();

      //setup system with eps0
      PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
      assemble_->AssembleMatrices(false);

      if(evalVersion != 3){
        /*
         * subtract solution
         */
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(STIFFNESS,solVec_,true);

        if(trans){
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
        }

        solVec_.ScalarMult(-1.0);
      }

    } else if (callingCnt == 1){
      /*
       * compute res1
       *
       * 3.0
       *  res^k+1 = res0 + sum_jj_from_2_to_k+1 ( deltaRes^jj )
       *
       *  deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
       *                           - usedEta^jj*eps0*deltaE^jj
       * 3.1
       *  res^k+1 = res0 + sum_jj_from_1_to_k+1 ( deltaRes^jj )
       *
       *  deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
       *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
       * 3.2
       *  res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
       * 3.3
       *  res1 = f_lin + f_nonlin(E_n+1^1) - div( eps0*E_n+1^1 ) - div( P_n+1^1 )
       *
       *  -> version 3.2 and 3.3 build up res completely anew whereas version
       *       3.1 and 3.0 already can start with updating
       *  -> the updating has the advantage, that the system which is used for
       *       the update ( [eps0 + deltaP^jj/deltaE^jj] for 3.1 and eps0 for 3.0 )
       *       is the same as the one which has to be solved during the next iteration
       *       -> avoid reassemble
       *  -> note: when using delta for updating, do not forget usedEta as scaling
       *        (this is included in the other versions, too as they use
       *          E_n+1^1 = E_n+1^0 + usedEta * deltaE^1
       */
      if( (evalVersion == 2)||(evalVersion == 3) ){
        /*
         * completely setup rhs
         */
        algsys_->InitRHS(RhsLinVal_); //f_lin

        // call to hysteresis function will return P
        // -> assembleNonLinRHS will also add P
        PDE_.SetFlagInCoefFncHyst("returnZeroValues",false);
        assemble_->AssembleNonLinRHS();

        //setup system with eps0
        PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
        assemble_->AssembleMatrices(false);

        /*
         * subtract solution
         */
        solVec_.ScalarMult(-1.0);
        algsys_->UpdateRHS(STIFFNESS,solVec_,true);

        if(trans){
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
        }

        solVec_.ScalarMult(-1.0);

      } else {
        /*
         * update residual with
         * 3.0
         *        deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
         *                        - usedEta^jj*eps0*deltaE^jj
         * 3.1
         *        deltaRes^jj = f_nonlin(E_n+1^jj) - f_nonlin(E_n+1^jj-1)
         *                        - usedEta^jj*[eps0 + deltaP^jj/deltaE^jj]*deltaE^jj
         */
        /*
         * subtract f_nonlin(E_n+1^jj-1) from res0
         */
        resVec_.Add( 1.0, resVec_, -1.0, nonLinRHS_);
        algsys_->InitRHS(resVec_);

        /*
         * add f_nonlin(E_n+1^jj); disable P
         */
        PDE_.SetFlagInCoefFncHyst("returnZeroValues",true);
        assemble_->AssembleNonLinRHS();

        if(evalVersion == 0){
          //setup system with eps0
          PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
          assemble_->AssembleMatrices(false);
        } else {
          //setup system with eps0 + deltaP^jj/deltaE^jj
          PDE_.SetFlagInCoefFncHyst("addFreeFieldTensorToDeltaMat",true);
          PDE_.SetFlagInCoefFncHyst("useNextToLastTS",false);
          PDE_.SetFlagInCoefFncHyst("useDeltaX",true);
          PDE_.SetFlagInCoefFncHyst("useDeltaY",true);
          assemble_->AssembleMatrices(false);
        }

        /*
         * subtract usedEta * solution update
         */
        if(usedEta != 0){
          solIncrement.ScalarMult(-1.0*usedEta);
          algsys_->UpdateRHS(STIFFNESS,solIncrement,true);

          if(trans){
            algsys_->UpdateRHS(STIFFNESS_UPDATE,solIncrement,true);
          }

          solIncrement.ScalarMult(-1.0/usedEta);
        }
      }

    } else {
      /*
       * during all later calls, residual is computed using
       *  res^k = res1 + sum_jj_from_2_to_k ( deltaRes^jj )
       *          = res^k-1 + deltaRes^k
       *
       * version0
       *  deltaRes^k = f_nonlin(E_n+1^k) - f_nonlin(E_n+1^k-1)
       *                        - usedEta^k*eps0*deltaE^k
       * else
       *  deltaRes^k = f_nonlin(E_n+1^k) - f_nonlin(E_n+1^k-1)
       *                        - usedEta^k*[eps0 + deltaP^k/deltaE^k]*deltaE^k
       */
      /*
       * subtract f_nonlin(E_n+1^jj-1) from res0
       */
      std::cout << "Calling cnt: " << callingCnt << std::endl;
      std::cout << "old residual: " << resVec_.ToString() << std::endl;
      resVec_.Add( 1.0, resVec_, -1.0, nonLinRHS_);
      std::cout << "old residual - nonLinRHS_: " << resVec_.ToString() << std::endl;
      algsys_->InitRHS(resVec_);

      /*
       * add f_nonlin(E_n+1^jj); disable P
       */
      PDE_.SetFlagInCoefFncHyst("returnZeroValues",true);
      assemble_->AssembleNonLinRHS();

      if(evalVersion == 0){
        //setup system with eps0
        PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
        assemble_->AssembleMatrices(false);
      } else {
        //setup system with eps0 + deltaP^jj/deltaE^jj
        PDE_.SetFlagInCoefFncHyst("addFreeFieldTensorToDeltaMat",true);
        PDE_.SetFlagInCoefFncHyst("useNextToLastTS",false);
        PDE_.SetFlagInCoefFncHyst("useDeltaX",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaY",true);
        assemble_->AssembleMatrices(false);
      }

      /*
       * subtract usedEta * solution update
       */
      if(usedEta != 0){
        solIncrement.ScalarMult(-1.0*usedEta);
        algsys_->UpdateRHS(STIFFNESS,solIncrement,true);

        if(trans){
          algsys_->UpdateRHS(STIFFNESS_UPDATE,solIncrement,true);
        }

        solIncrement.ScalarMult(-1.0/usedEta);
      }
    }

    /*
     * common steps and all methods and callingCnt (if needed at all)
     */
    /*
     * update with stage rhs
     */
    if(trans){
      /*
       * TODO: a) find out what stageRHS_ is used for
       *       b) is this update needed for hysteresis case?
       *       c) what matrices have to be used for this kind of update?
       *          do we need the linear matrices (with eps0) or the ones we
       *          finally solve (with eps0+delta)?
       */
      std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
      std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
      std::map<FEMatrixType,Integer>::iterator matIt;
      UInt pos = 0;
      //now update RHS according to time stepping
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){

        if(matIt->second < 0)
          continue;
        for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(pos));
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
      }
    }

    /*
     * save residual
     */
    std::cout << "Save residual" << std::endl;
    algsys_->GetRHSVal( resVec_ );
    std::cout << "residual: " << resVec_.ToString() << std::endl;

    /*
     * 2. ASSEMBLE SYSTEM FOR NEXT ITERATION
     */
    if(callingCnt == 0){
      /*
       * 1. the first call to this function in each time step is done BEFORE solving
       *     the system
       * 2. to (hopefully) improve convergence, the first system gets assembled
       *      with eps' = eps_initial instead of eps0
       */
      if(evalVersion == 0){
        /*
         * div( [eps' * deltaE^k+1 )
         * -> only eps' = eps_initial needed
         */
        PDE_.SetFlagInCoefFncHyst("returnInitialTensor",true);
        //PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
      } else if (evalVersion == 1){
        /*
         * div( [eps' + deltaP^k/deltaE^k]*deltaE^k+1 )
         * -> compute deltaMatrix and add eps_initial
         * -> during the first call to this function, there exists no
         *    old iteration value, yet
         *    -> compute deltaP and deltaE using the value from the
         *        last and next to last timestep
         *        deltaE^0 = E_n^final - E_n-1^final
         */
        PDE_.SetFlagInCoefFncHyst("addInitialTensorToDeltaMat",true);
        PDE_.SetFlagInCoefFncHyst("useNextToLastTS",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaX",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaY",true);
      } else if (evalVersion == 2){
        /*
         * div( [eps' + P_n+1^0/deltaE^0]*deltaE^1 )
         * -> compute deltaMatrix and add eps_initial
         * -> during the first call to this function, there exists no
         *    old iteration value, yet
         *    -> compute deltaP and deltaE using the value from the
         *        last and next to last timestep
         *        deltaE^0 = E_n^final - E_n-1^final
         */
        PDE_.SetFlagInCoefFncHyst("addInitialTensorToDeltaMat",true);
        PDE_.SetFlagInCoefFncHyst("useNextToLastTS",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaX",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaY",false);
      } else if (evalVersion == 3){
        /*
         * div( [eps' + P_n+1^0/E_n+1^0]*E_n+1^1 )
         * -> compute deltaMatrix and add eps_initial
         */
        PDE_.SetFlagInCoefFncHyst("addInitialTensorToDeltaMat",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaX",false);
        PDE_.SetFlagInCoefFncHyst("useDeltaY",false);
      } else {
        std::cout << "EvalVersion: " << evalVersion << std::endl;
        EXCEPTION("Unknown evalVersion");
      }
      assemble_->AssembleMatrices(false);
    } else if (callingCnt == 1){
      if((evalVersion == 2)||(evalVersion == 3)){
        //setup system with eps0 + deltaP^jj/deltaE^jj
        PDE_.SetFlagInCoefFncHyst("addFreeFieldTensorToDeltaMat",true);
        PDE_.SetFlagInCoefFncHyst("useNextToLastTS",false);
        PDE_.SetFlagInCoefFncHyst("useDeltaX",true);
        PDE_.SetFlagInCoefFncHyst("useDeltaY",true);
        assemble_->AssembleMatrices(false);
      }
      /*
       * for evalVersion0 and evalVersion1, the system was already
       * build up during the calculation of the residual
       */
      //assemble_->AssembleMatrices(false);
    }
    /*
     * for all later iterations, we need no reassmeble
     * as the residual gets computed using deltaP/deltaE, too
     */

    /*
     * 3. WRAP UP
     */

    /*
     * set value for nonlinear rhs without hysteresis which is needed
     * during the next iteration
     */
    nonLinRHS_ = tmpRHS;

    std::cout << "wrap up residual: " << resVec_.ToString() << std::endl;

  }

  void StdSolveStep::ConfigureSystemForHysteresisResidual(SBM_Vector& oldSolution,SBM_Vector& solIncrement,Double usedEta, UInt stage, bool trans){
    /*
     *  Normal way to calculate residual:
     *  1. get new solution by updating old one with solution update
     *      u_new = u_old + deltaU
     *  2. setup nonlinear system
     *      K(u_new)
     *  3. check if K(u_new)*u_new = f
     *
     *  Problem with delta hysteresis stepping (for example of electrostatics):
     *    Our system looks like
     *      div (eps_0 + P/deltaE)*deltaE = f- div(eps_0 * E_old)
     *    -> we have different matrices on rhs and lhs
     *    If we compute rhs as above, the residual check would be
     *      res = f - div(eps_0 * E_new) = f - div(eps_0*E_old) - div(eps_0*deltaE)
     *    Obviously, we are missing the term div(P) on the rhs
     *
     *    If we compute rhs using the lhs matrix, we would get
     *      res = f - div (eps_0 + P/deltaE)*(E_old+deltaE)
     *          = f - div (eps_0 + P/deltaE)*E_old + div (eps_0 + P/deltaE)*deltaE
     *          = f - div (eps_0*E_old) + div(eps_0+deltaE) + div(P) + div(P/deltaE*E_old)
     *    Now the term + div(P/deltaE*E_old) is too much
     *
     *    Idea (and reason for this extra function):
     *      1. build up rhs matrix using only eps0
     *      2. calculate (with old solution only! Otherwise we have the term eps0*deltaE twice!)
     *          res = f - div(eps_0 * E_old)
     *      3. build up lsh matrix using E_new
     *      (so far, these steps are all done in ConfigureSystemForHysteresis)
     *      4. (NEW) update res with deltaE (has to be passed to this function)
     *          res -= div(eps_0 + P/deltaE)*deltaE
     *
     *      -> res = f - div(eps_0*E_old) - div(eps_0*deltaE) - div(P_new)
     *
     *
     */

     /*
      * get current new solution
      */
     if(usedEta != 0){
       solVec_.Add( 1.0, oldSolution, usedEta, solIncrement);
     }

     SBM_Vector actRHS_TMP(BaseMatrix::DOUBLE);

     algsys_->InitRHS(RhsLinVal_);

     assemble_->AssembleNonLinRHS();

     std::cout << "RHS after AssembleNonLinRHS: " << std::endl;
     algsys_->GetRHSVal( actRHS_TMP );
     std::cout << actRHS_TMP.ToString() << std::endl;

     //std::cout << "Deactivate Delta Computation" << std::endl;
     PDE_.ActivateDeactivateDeltaMat(false);

     // here we get the system matrix using only eps0
     assemble_->AssembleMatrices(false);

     oldSolution.ScalarMult(-1.0);
     algsys_->UpdateRHS(STIFFNESS,oldSolution,true);

     if(trans){
       algsys_->UpdateRHS(STIFFNESS_UPDATE,oldSolution,true);
     }

     oldSolution.ScalarMult(-1.0);

     if(trans){

       std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
       std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
       std::map<FEMatrixType,Integer>::iterator matIt;
       UInt pos = 0;
       //now update RHS according to time stepping
       //TODO: find out for which matrices this updating has to be done, i.e.
       //      do we need the temporary matrices (with eps0) or the ones we finally solve (with eps0+delta)
       for(matIt = matrices.begin();matIt != matrices.end();matIt++){

         if(matIt->second < 0)
           continue;
         for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
           fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(pos));
         }
         algsys_->UpdateRHS(matIt->first,stageRHS_,true);
       }

     }

     PDE_.ActivateDeactivateDeltaMat(true);

     /*
      * new we get the matrices with a deltaMaterial tensor
      */
     assemble_->AssembleMatrices(false);

     std::cout << "RHS after Update RHS with sol 2: " << std::endl;
     algsys_->GetRHSVal( actRHS_TMP );
     std::cout << actRHS_TMP.ToString() << std::endl;
     std::cout << "L2Norm: " << actRHS_TMP.NormL2() << std::endl;

     /*
      * NEW step 4
      */
     solIncrement.ScalarMult(-1.0);
     algsys_->UpdateRHS(STIFFNESS,solIncrement,true);

     if(trans){
       algsys_->UpdateRHS(STIFFNESS_UPDATE,solIncrement,true);
     }

     std::cout << "RHS after Update RHS with sol increment: " << std::endl;
     algsys_->GetRHSVal( actRHS_TMP );
     std::cout << actRHS_TMP.ToString() << std::endl;
     std::cout << "L2Norm: " << actRHS_TMP.NormL2() << std::endl;
     std::cout << "SolInc: " << solIncrement.ToString() << std::endl;

     solIncrement.ScalarMult(-1.0);

   }


  void StdSolveStep::ConfigureSystemForHysteresis(UInt stage, bool trans, bool firstTime){
    /*!
     * Important note: The following system will assemble the rhs and lhs
     * of the equation system as needed for StepTransNonLinHysteresis().
     *
     * Both lhs and rhs will depend on the current solution as well as on the
     * previous solution (from the previous iteration).
     *
     * The current solution is thereby taken from the class variable
     *  solVec_
     * that is linked to the solution storage in the PDEs.
     *
     * The previous solution is stored by CoefFunctionHyst that is stored by
     * the PDEs. To set the previous values, call the function
     *  PDE_.SetPreviousHystVals()
     *
     */

    /*
     * System to be configured (here at the example of electrostatics)
     *    div([eps0 + deltaP^k/deltaE^k]*deltaE^k+1) = f_n+1 - div(eps0*E_n+1^k) - div(P_n+1^k)
     */
    /*
     *    a,1) build up rhs
     *      algsys_->InitRHS(RhsLinVal_)             -> f_n+1
     *      assemble_->AssembleNonLinRHS()           -> -div(P_n+1^k), has to be set as a rhs-term in elecPDE!
     *
     *      PDE_->ActivateDeactivateDeltaMat(false)  -> tell CoefFncHyst to return only eps0 as material tensor
     *
     *      assemble_->AssembleMatrices(false)       -> build up system matrix with eps0 as material relation
     *      solVec_.ScalarMult(-1.0);
     *      algsys_->UpdateRHS(STIFFNESS,solVec_,true);
     *      algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
     *      solVec_.ScalarMult(-1.0);
     *                                               -> -div(eps0*E_n+1^k)
     */

//    std::cout << "RHS at begin of ConfigureSystemForHysteresis: " << std::endl;
    SBM_Vector actRHS_TMP(BaseMatrix::DOUBLE);
//    algsys_->GetRHSVal( actRHS_TMP );
//    std::cout << actRHS_TMP.ToString() << std::endl;

    //std::cout << "Setup Lin RHS" << std::endl;

    algsys_->InitRHS(RhsLinVal_);

    std::cout << "rhsVec_: " << rhsVec_.ToString() << std::endl;
    std::cout << "RhsLinVal_: " << RhsLinVal_.ToString() << std::endl;



//    std::cout << "RHS after InitRHS: " << std::endl;
//    algsys_->GetRHSVal( actRHS_TMP );
//    std::cout << actRHS_TMP.ToString() << std::endl;

    //std::cout << "Setup NonLin RHS" << std::endl;

    assemble_->AssembleNonLinRHS();

    std::cout << "RHS after AssembleNonLinRHS: " << std::endl;

    algsys_->GetRHSVal( actRHS_TMP );
    std::cout << actRHS_TMP.ToString() << std::endl;

    //std::cout << "Deactivate Delta Computation" << std::endl;
    PDE_.ActivateDeactivateDeltaMat(false);

    /*
     * Note: we do not assemble Newton parts at all, as we do not have information about the
     * derivative of the HysteresisOperator
     */
    assemble_->AssembleMatrices(false);

//    std::cout << "solVec_ before updating rhs" << std::endl;
//    std::cout << solVec_.ToString() << std::endl;

    std::cout << "solVec_: " << solVec_.ToString() << std::endl;
    solVec_.ScalarMult(-1.0);
    //std::cout << "Update RHS" << std::endl;
    algsys_->UpdateRHS(STIFFNESS,solVec_,true);

    //std::cout << "RHS after Update RHS with sol 1: " << std::endl;
    //algsys_->GetRHSVal( actRHS_TMP );
    //std::cout << actRHS_TMP.ToString() << std::endl;

    if(trans){
      algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
    }

    std::cout << "RHS after Update RHS with sol 2: " << std::endl;
    algsys_->GetRHSVal( actRHS_TMP );
    std::cout << actRHS_TMP.ToString() << std::endl;

    solVec_.ScalarMult(-1.0);

    if(trans){

      std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
      std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
      std::map<FEMatrixType,Integer>::iterator matIt;
      UInt pos = 0;
      //now update RHS according to time stepping
      //TODO: find out for which matrices this updating has to be done, i.e.
      //      do we need the temporary matrices (with eps0) or the ones we finally solve (with eps0+delta)
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){

        if(matIt->second < 0)
          continue;
        for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
          fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(pos));
        }
        algsys_->UpdateRHS(matIt->first,stageRHS_,true);
      }

//      std::cout << "RHS after Update RHS with stageRHS: " << std::endl;
//      algsys_->GetRHSVal( actRHS_TMP );
//      std::cout << actRHS_TMP.ToString() << std::endl;

    }

    /*    a,2) build up lhs
     *      PDE_->ActivateDeactivateDeltaMat(true)   -> tell CoefFncHyst to return only eps0+deltaP^k/deltaE^k
     *                                                    as material tensor
     *      assemble_->AssembleMatrices(false)       -> rebuild up system matrix with eps0+deltaP^k/deltaE^k
     *                                                    as material relation
     */
    //std::cout << "Activate Delta Computation" << std::endl;
    if(firstTime){
      PDE_.UseNextToLastTSForDeltaMat(true);
    } else {
      PDE_.UseNextToLastTSForDeltaMat(false);
    }

    //if(!firstTime){
    /*
     * use always deltaMatrix (also for first iteration), but for first iteration
     * we use the nextToLast value instead of the last iteration value
     */
      PDE_.ActivateDeactivateDeltaMat(true);
    //}

    //std::cout << "Re-Assemble Matrices" << std::endl;
    assemble_->AssembleMatrices(false);

  }

  Double StdSolveStep::LineSearchHyst(SBM_Vector& solIncrement, Double& etaLineSearch, UInt evalVersion,
                                      UInt callingCnt, bool trans, bool performLineSearch)  {

    /*
     * backup nonlinear rhs, residual and solution vector
     */
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = solVec_;

    SBM_Vector resOld(BaseMatrix::DOUBLE);
    resOld = resVec_;

    SBM_Vector nonLinOld(BaseMatrix::DOUBLE);
    nonLinOld = nonLinRHS_;

    const UInt nrEtas = 9;
    Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01,0.00316,0.001,0.000316,0.0001};

    bool assumeQuadratic = true;

    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    UInt nrEtasTest = 9;
    if(performLineSearch == false){
      /*
       * only take first eta = 1
       * -> we do not have to lock hysteresis nor do we have to backup
       *    the solution or the residual as we directly take the first value
       */
      etaOpt = 1.0;
    } else {
      /*
       * lock hysteresis memory as we otherwise would store the
       * effects of all tested eta
       */
      PDE_.SetFlagInCoefFncHyst("lockHysteresisMemory",true);

      /*
       * test different eta
       */
      for( UInt i=0; i<nrEtasTest; i++) {

        std::cout << "CalcResidualAndConfigSystemForHysteresis during linesearch" << std::endl;
        CalcResidualAndConfigSystemForHysteresis(solVec_,solIncrement,eta[i], 0, callingCnt, evalVersion, trans);

        /*
         * calculate error norm as criterion
         */
        Double residualL2Norm = resVec_.NormL2();

        /*
         * restore old solution vector, old residual and nonlinear rhs
         * -> needed to perform a correct next call to CalcResidual...
         */
        solVec_ = solOld;

        resVec_ = resOld;

        nonLinRHS_ = nonLinOld;

        if (residualL2Norm < residualL2NormOpt) {
          residualL2NormOpt = residualL2Norm;
          etaOpt = eta[i];
        } else {
          if(assumeQuadratic){
            /*
             * here we assume that the residual follows a quadratic curve
             * i.e. we only have one minimum which is the global minimum;
             * in that case, if the residual increases between two iterations,
             * we can assume that we passed the minimum already and stop
             */
            break;
          }
        }
      }

      /*
       * unlock hysteresis memory
       */
      PDE_.SetFlagInCoefFncHyst("lockHysteresisMemory",false);
    }

    /*
     * now we have found etaOpt, so we can call CalcResidualAndConfigSystemForHysteresis
     * one last time to setup the correct system for next computation
     *
     */
    etaLineSearch = etaOpt;

    std::cout << "CalcResidualAndConfigSystemForHysteresis after linesearch" << std::endl;
    std::cout << "residual at end of linesearch pre call: " << resVec_.ToString() << std::endl;
    CalcResidualAndConfigSystemForHysteresis(solVec_,solIncrement,etaLineSearch, 0, callingCnt, evalVersion, trans);
    std::cout << "residual at end of linesearch after call: " << resVec_.ToString() << std::endl;


    return resVec_.NormL2();
  }

  void StdSolveStep::StepTransNonLinHysteresis() {
    /*
     * Idea: for deltaMaterial Hysteresis, use base code of StepTransNonLin() as it also uses
     *       incremental stepping
     *       However, it is adapted as described below
     */
    /*! Basic idea shown for electrostatics (D = flux density, E = electric field, P = polarization
     *      eps0 = permittivity of vacuum, n = time step counter, k = iteration counter, f = rhs load vector):
     *
     *    D_n+1 = D_n + deltaD
     *          = eps0*E_n + P_n + esp0*deltaE + deltaP
     *          = (eps0 + P_n/E_n)*E_n + (eps0 + deltaP/deltaE)*deltaE
     *          = eps(E_n)*E_n + eps(deltaE)*deltaE
     *
     *   > Note that unlike the normal nonlinear case, we have two different
     *     "eps" here due to the different splitting approach;
     *     for comparison
     *        a) nonlin,nonhyst (see StepTransNonLin()):
     *                D_n+1 = eps(E_n+1)*E_n+1
     *                      = eps(E_n+1)*(deltaE + E_n)
     *                      ~ eps(E_n)*(deltaE + E_n)
     *        b) nonlin,hyst:
     *                D_n+1 = eps(E_n)*E_n + eps(deltaE)*deltaE
     *
     *
     *  Solution strategy:
     *    0. Initialize D_0^0, E_0^0, P_0^0, E_0^-1,P_0^-1
     *
     *    1. At beginning of time step n+1 retrieve initial values
     *        E_n+1^0 = E_n^kfinal     -> stored as actSol (in fact potential is stored, but E derives directly from it)
     *        E_n+1^-1 = E_n^kfinal-1  -> stored in CoefFncHyst as previousHystValue
     *        deltaE^0 = E_n+1^0 - E_n+1^-1
     *        P_n+1^0 = HystOperator(E_n+1^0)
     *        P_n+1^-1 = HystOperator(E_n+1^-1)
     *        deltaP^0 = P_n+1^0 - P_n+1^-1
     *
     *    2. Compute E_n+1^k+1
     *    a) build up system using actSol E_n+1^k and previousSol E_n+1^k-1
     *
     *        div(D_n+1^k+1) = f_n+1
     *     >  div(deltaD^k+1 + D_n+1^k) = f_n+1
     *     >  div(deltaD^k+1) = f_n+1 - div(D_n+1^k)
     *     >  div([eps0 + deltaP^k/deltaE^k]*deltaE^k+1) = f_n+1 - div(eps0*E_n+1^k + P_n+1^k)
     *     >  div([eps0 + deltaP^k/deltaE^k]*deltaE^k+1) = f_n+1 - div(eps0*E_n+1^k) - div(P_n+1^k)
     *
     *    a,1) build up rhs
     *      algsys_->InitRHS(RhsLinVal_)             -> f_n+1
     *      assemble_->AssembleNonLinRHS()           -> div(P_n+1^k), has to be set as a rhs-term in elecPDE!
     *
     *      PDE_->ActivateDeactivateDeltaMat(false)  -> tell CoefFncHyst to return only eps0 as material tensor
     *
     *      assemble_->AssembleMatrices(false)       -> build up system matrix with eps0 as material relation
     *      solVec_.ScalarMult(-1.0);
     *      algsys_->UpdateRHS(STIFFNESS,solVec_,true); // we also or only need the updated version
     *      algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
     *      solVec_.ScalarMult(-1.0);
     *           -> div(eps0*E_n+1^k)
     *
     *    a,2) build up lhs
     *      PDE_->ActivateDeactivateDeltaMat(true)   -> tell CoefFncHyst to return only eps0+deltaP^k/deltaE^k
     *                                                    as material tensor
     *      assemble_->AssembleMatrices(false)       -> rebuild up system matrix with eps0+deltaP^k/deltaE^k
     *                                                    as material relation
     *
     *    b) store E_n+1^k as old solution by calling SetPreviousHystVals (overwriting E_n+1^k-1)
     *       -> after the solution of the system, it will be too late to call this function!
     *       -> Reason: SetPreviousHystVals will use the CURRENT solution as input to the hysteresis
     *                  operator and store it together with the output as previousHystValues
     *                  -> has to be called as long as E_n+1^k is still the CURRENT solution!
     *
     *    c) solve system
     *
     *    c,1) if k = 1 (first iteration)
     *
     *        call ClearIDBCFromSolutionVal to remove the inhom. Dirichlet values from old solution
     *          -> without calling this function, we would add up the IDBC values from all time steps
     *          -> see comment at start of StepTransNonLin(), implementation III
     *
     *        set IDBC=true then solve system -> solution will have current IDBC values
     *
     *        add new solution (increment) to old solution
     *
     *    c,2) if k > 1
     *
     *        set IDBC=false then solve system -> solution will have 0 at IDBC nodes
     *
     *        add new solution (increment) to old solution
     *
     *    3. Use linesearch to fine optimal stepping value eta (except for first iteration where
     *        we need a stepping of 1 to get correct IDBC values)
     *
     *        E_n+1^k+1 = E_n+1^k + eta*deltaE^k+1
     *
     *    4. Check for convergence
     *
     *        div([eps0 + deltaP^k+1/deltaE^k+1]*deltaE^k+1) = f_n+1 - div(eps0*E_n+1^k+1) - div(P_n+1^k+1)
     *
     *    5. If convergence = true
     *        PDE_->ActivateDeactivateDeltaMat(false)
     *          -> tell CoefFncHyst to return only eps0 as material tensor
     *          -> this is done for the calculation of D_n+1^kfinal which is done
     *              at the end of the timestep; if we would have the delta computation activated,
     *              D would be calculated using a wrong eps;
     *              Note that we have to add the value of P to D in elecPDE!
     *
     *
     *  Note: The splitted computation approach, i.e. having deltaP/deltaE on lhs and eps0 E + P on rhs,
     *        instead of using deltaP/deltaE or P/E on rhs is the following:
     *
     *    a) using deltaP/deltaE simply is not correct from the mathematical side
     *    b) using P/E has the issue that it cannot handle remanence; even if E is not 0
     *        but only close, the value of P/E will be enormous and most probably lead to bad behavior
     *
     */

    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransNonLinHysteresis";
    bool performOneMoreStep;

    SBM_Vector solInc(BaseMatrix::DOUBLE);

    //get actual solution
    // not used -> write directly to solVec_
//    SBM_Vector  actSol(BaseMatrix::DOUBLE);
//    actSol = solVec_;

    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 ){
      /*
       * maybe it is quite easy to enable multiple stages, maybe it does not work;
       * feel free to find out
       */
       EXCEPTION("StepTransNonLinHysteresis: just one stage time-stepping allowed");
    }

    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;

    UInt pos = 0;

    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }

    //do initialization
    rhsVec_.Init();

    LOG_DBG(stdsolvestep) << "StepTransNonLinHysteresis: Stage: " << 0 ;

    /*
     * stage solution not needed, as we only allow one stage, stageSol will
     * be solVec_
     */
//    //we obtain a reference to the stage vectors of the scheme
//    SBM_Vector stageSol;
//    stageSol.Resize(feFunctions_.size());
//    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(0),pos);
//      fncIt->second->GetTimeScheme()->InitStage(0,actTime_,PDE_.GetDomain());
//    }
//    stageSol.SetOwnership(false);
//
//    /*
//     * During the later computations, we will add all incremental solutions to
//     * stageSol; unlike the case of StepTransNonLin, we add the results to the
//     * current solution instead of only gathering the increments.
//     * Note that we initialize stageSol with actSol in StepTransNonLin, too, but
//     * there, we will Init stageSol with 0.0 during the first iteration
//     */
//    stageSol = actSol;

    /*
     * at this point, stageSol, actSol and solVec_ all have the same value!
     * Note that solVec_ is pointing to the results which are stored in the PDEs
     * i.e. changes done to this vector will update the actual solution which
     * is handeled by the PDEs itself!
     */

    // setup right hand side
    Double loadFactor = 1.0;
    Double RhsLinL2Norm = SetLinRHS(loadFactor);

    // set iteration counter
    UInt iterationCounter=0;
    UInt callingCnt = 0;

    UInt evalVersion;
    if(PDE_.IsHysteresis_Fixpoint() == true){
      evalVersion = 0;
    } else {
      evalVersion = 1;
    }

    bool deltaIDBC = false;
    //Double oldIncError = 0.0;

    /*
     * disable output of switching and rotation state as BMP images (even if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",false);

    do {
      iterationCounter++;

      LOG_DBG(stdsolvestep) << "StepTransNonLinHysteresis: Iteration: " << iterationCounter ;
      std::cout << "StepTransNonLinHysteresis: Iteration: " << iterationCounter << std::endl;

//      std::cout << "RHS at begin of iteration: " << std::endl;
//      SBM_Vector actRHS_TMP(BaseMatrix::DOUBLE);
//      algsys_->GetRHSVal( actRHS_TMP );
//      std::cout << actRHS_TMP.ToString() << std::endl;
//

      if( iterationCounter == 1){
        /*
         * During the first iteration, we need an initial setup of the algebraic system.
         * Otherwise, we would have nothing to solve.
         * Note that we rebuild the system later on when we compute the residual.
         * This is the reason, why we do not have to build it at the beginning of all
         * further iterations (we simply resuse the system from the residual computation).
         */
        PDE_.SetFlagInCoefFncHyst("lockHysteresisMemory",false);

        /*
         * configure system via CalcResidualAndConfigSystemForHysteresis event though
         * we do not need the residual right now
         */
        CalcResidualAndConfigSystemForHysteresis(solVec_,solInc,0.0, 0, callingCnt, evalVersion, true);
        callingCnt++;

        if(evalVersion == 3){
          /*
           * here we perform a complete full step during the first iteration
           * -> see CalcResidualAndConfigSystemForHysteresis for details
           * -> we have to set full idbc values for solving the system
           *
           */
          deltaIDBC = false;
        } else {
          /*
           * as we subtract the old solution from the rhs, we also have to
           * subtract the old IDBC values
           * -> set deltaIDBC to true
           */
          deltaIDBC = true;
        }

      }
      /*
       *    b) store E_n+1^k as old solution by calling SetPreviousHystVals (overwriting E_n+1^k-1)
       *       -> after the solution of the system, it will be too late to call this function!
       *       -> Reason: SetPreviousHystVals will use the CURRENT solution as input to the hysteresis
       *                  operator and store it together with the output as previousHystValues
       *                  -> has to be called as long as E_n+1^k is still the CURRENT solution!
       */

      /*
       * to set previous hysteresis values, we have to enable their calculation first
       */
      PDE_.SetFlagInCoefFncHyst("returnZeroValues",false);
      if(iterationCounter == 1){
        /*
         * During the first iteration, we store the current input (=solution of last timestep)
         * for the next timestep as value nextToLastTS (after all it will be the next to last
         * vector then)
         */
        PDE_.SetPreviousHystVals(true);
      }
      /*
       * store current solution before solving the system
       * -> previousItValue
       */
      PDE_.SetPreviousHystVals(false);

      // set system matrix to zero initially, as ConstructEffectiveMatrix only
      // sums up the contributions
      matrix_factor_.clear();
      algsys_->InitMatrix(SYSTEM);
      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
        FeFctIdType fctId = fncIt->second->GetFctId();
        fncIt->second->GetTimeScheme()->AddMatFactors(0,matrices,matrix_factor_[fctId]);
        algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
      }
      // setup the matrices to compute correct error norms

      PDE_.SetBCs();
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond();
      algsys_->SetupSolver();

      // just set inh. Dirichlet BCs for the first iteration
      bool setIDBC = false;

      if ( iterationCounter == 1 && couplingIter_ == 0 ){
        setIDBC = true;
      }

      /*
       * solve system; if IDBC=true, the values of the IDBC nodes are set in solInc
       *
       * ATTENTION - POSSIBLE HUGE ERROR - SOURCE:
       *  we must not set the original IDBC values but the delta to the previous
       *  time step!
       *  we set up the system such that we only solve for the delta (i.e. we
       *  subtract the old solution via the rhs); this will gives us an update
       *  for all inner nodes BUT ONLY if we solve with a delta of the IDBC
       *  values; if we simply set the IDBC values of the new time step, that
       *  will be seen as the delta values at the IDBC nodes;
       *  as long as we do not have other loadings (i.e. f = 0) the resulting
       *  value for sol_inc will simply be the full solution
       *  (-> as we subtracted nothing from the rhs for f = 0 and assuming that
       *  the old solution also satisfied this, i.e. divD = 0, we simply solve
       *  for divD_inc = 0; the resulting solution will of course satisfy that
       *  condition but for the original IDBC value, so that D_inc is in fact
       *  D_new; that is ok as we simply can take the new solution then; BUT
       *  if we have a rhs load, than the results won't match!
       *  we would solve for
       *  divD_new = delta_f
       *  two solutions:
       *  1. do not subtract old solution during first iteration, i.e.
       *    solve divD_new = f in first step, then update
       *  2. subtract old solution during first iteration, but then we also
       *    have to subtract the old idbc values!l
       *
       *  -> 1. used for evalVerson == 3
       *  -> 2. used otherwise
       *
       */
      algsys_->Solve(setIDBC,deltaIDBC);
      algsys_->GetSolutionVal(solInc, setIDBC, deltaIDBC );

      Double residualL2Norm = 0.0;
      Double etaLineSearch  = 1.0;

      /*
       * In StepTransNonLin, we would set stageSol to 0 now, but here, we don't!
       * -> as we have eliminated the old IDBC values from solVec_ we can simply add
       * the solution increment to it
       */
//        if ( iterationCounter == 1 && couplingIter_ == 0 )
//          stageSol.Init();

      /*
       * Note: LineSearchHyst will call CalcResidualAndConfigSystemForHysteresis
       * -> solVec_ will be set to actual solution
       * -> resVec_ and nonLinRHS_ will be set to current values
       * -> system will be configured for next iteration
       */

      if ( lineSearch_ == "none" || iterationCounter == 1) {
        /*
         * set flag performLineSearch to false
         * -> LineSearchHyst will only test eta=1 and keep that value
         */
        residualL2Norm = LineSearchHyst(solInc, etaLineSearch, evalVersion, callingCnt, true,false);
      } else {
        /*
         * set flag performLineSearch to true
         * -> LineSearchHyst will test various eta and keep the one which mini-
         *      mizes the residual
         */
        residualL2Norm = LineSearchHyst(solInc, etaLineSearch, evalVersion, callingCnt, true,true);
      }
      /*
       * callingCnt is the number of actual residuals to be kept
       * -> temporal residuals during lineSearch shall not increase the counter
       */
      callingCnt++;

      /*
       * now we have the updated solution
       */
      std::cout << "New solVec_ " << solVec_.ToString() << std::endl;

      // calculation of relative residual error =======================================
      Double residualErr;
      if ( RhsLinL2Norm > 1.0 ){
        residualErr = residualL2Norm / RhsLinL2Norm;
      } else {
        residualErr = residualL2Norm;
      }

      // calculate incremental error ========================================
      Double incrementalErr;
      //TODO: we should have eta*solInc.NormL2() as we do not perform full steps in general
      //Double solIncrL2Norm = solInc.NormL2();
      Double solIncrL2Norm = etaLineSearch*solInc.NormL2();
      Double actSolL2Norm  = solVec_.NormL2();

//      std::cout << "UsedEta: " << etaLineSearch << std::endl;
//      std::cout << "Norm of incremental: " << etaLineSearch*solIncrL2Norm << std::endl;
//      std::cout << "Norm of new solution: " << actSolL2Norm << std::endl;

      if ( actSolL2Norm > 1.0 ){
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      } else {
        incrementalErr = solIncrL2Norm;
      }

      //Double incStopCritABS_ = 1e-2;
      //Double incrementalErrABS = solIncrL2Norm;
      std::cout << "incrementalErr_abs: " << solIncrL2Norm << std::endl;
      std::cout << "incrementalErr_rel: " << incrementalErr << std::endl;
      std::cout << "residualErr: " << residualErr << std::endl;

      // output of norms and data
      if ( nonLinLogging_ == true ) {
        // get current step
        UInt actStep = PDE_.GetSolveStep()->GetActStep();

        if (PDE_.IsIterCoupled()) {
          WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
        } else {
          //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
          WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
        }
        // write norm to file
        logFile_ <<  iterationCounter << "\t"
            << residualErr << "\t"
            << incrementalErr << "\t"
            << etaLineSearch << std::endl;
      }

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

    /*
     * At the end of the time step, we have to tell the hysteresis operator
     * to return only eps0 (nu0 for magnetic) as we compute
     *  D = eps0*E + P
     * during post processing
     */
    PDE_.SetFlagInCoefFncHyst("returnFreeFieldTensor",true);
    /*
     * to calculate P during output, we have to enable the output of P of course
     */
    PDE_.SetFlagInCoefFncHyst("returnZeroValues",false);
    /*
     * allow to output switching and rotation state as BMP images (if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",true);

    /*
     * Store IDBC values from the current time step in form of its rhs representation
     * (i.e. the effect that it will have on the rhs)
     */
    algsys_->SetOldDirichletValues();

    std::cout << "EndOfStepTransNonLinHysteresis" << std::endl;

  }
  
// old backup
//  void StdSolveStep::StepTransNonLinHysteresis() {
//    /*
//     * Idea: for deltaMaterial Hysteresis, use base code of StepTransNonLin() as it also uses
//     *       incremental stepping
//     */
//
//    LOG_TRACE(stdsolvestep) << "StdSolveStep::StepTransNonLinHysteresis";
//    bool performOneMoreStep;
//    bool isNewton = false;
//
//    SBM_Vector solInc(BaseMatrix::DOUBLE);
//
//    //get actual solution
//    SBM_Vector  actSol(BaseMatrix::DOUBLE);
//    actSol = solVec_;
//
//    //obtain the number of stages
//    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
//
//    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
//    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
//    std::map<FEMatrixType,Integer>::iterator matIt;
//
//    UInt pos = 0;
//
//    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
//    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
//      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
//    }
//
//    for(UInt i=0;i<numStages;i++){
//      //do initialization
//      rhsVec_.Init();
//      LOG_DBG(stdsolvestep) << "StepTransNonLinHysteresis: Stage: " << i ;
//
//      //we obtain a reference to the stage vectors of the scheme
//      SBM_Vector stageSol;
//      stageSol.Resize(feFunctions_.size());
//      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),pos);
//        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
//      }
//      stageSol.SetOwnership(false);
//
//
//      //initialize solution vector for each stage
//      if ( i > 0 )
//        actSol = stageSol;
//      else{
//        //special case of incremental non-linearity, we set the stage vector to the solution vector
//        stageSol = actSol;
//      }
//
//      /*
//       * Note: solVec_ stores pointer to FeFunctions_ which itself store the results
//       *       -> although assemble and algsys never access this variable, its value influences
//       *          the algebraic system as assemble (and algsys?) access the FeFunctions to retrieve
//       *          the results needed for the computation of material non-linearities for example
//       */
//      solVec_  = actSol;
//
//      // setup right hand side
//      Double loadFactor = 1.0;
//      Double RhsLinL2Norm = SetLinRHS(loadFactor);
//
//      // set iteration counter
//      UInt iterationCounter=0;
//      Double oldIncError = 0.0;
//
//      // FinalizeAfterTimestep sets old hysteresis values to 0 so that deltaP and deltaE is basically
//      // P_current/E_current
//      //PDE_.FinalizeAfterTimeStep();
//
//      do {
//        iterationCounter++;
//
//        std::cout << "############# ITERATION " << iterationCounter << "############" << std::endl;
//
//        /*
//         * in StepTransNonLin, we perform the setup step at multiple points of this loop
//         * 1. the first time is at the beginning of the loop but only for two cases:
//         *      - first overall iteration:
//         *          The assembly during first iteration is clear, as we do not have any system to solve yet
//         *
//         *      - no lineSearch activated:
//         *          The question is why not also for the case without lineSearch?
//         *          The reason is, that we later on update the whole system with the new solution in
//         *          order to calculate the residual, in that sense, it would be unnecessary to update
//         *          it with the same values a second time.
//         *          In case of lineSearch, we set up the system and calculate the residual in the lineSearch
//         *          function, again after the solution was computed. However, we try different step width
//         *          in the lineSearch function and first after testing all possible steps, we calculate the
//         *          real new solution. Only if the chosen stepwidth is the last in the list of possibilities,
//         *          the algebraic system will be updated with the correct value at the end of the do loop.
//         *          We have now two possibilities:
//         *            a) directly update the system with the new found solution
//         *            b) update system with new found solution at the beginning of the next iteration
//         *          Clearly, option b) was chosen for StepTransNonLin. For this function, however, option
//         *          a) was taken for the following reason:
//         *            If a solution was found and, thus, no new iteration is needed, the solution and the
//         *            last assembled equation system do not fit. This is not problem as we would update it
//         *            either way during the next timestep, but in order to check if the solution has changed
//         *            since last time, it would be nice if the system was updated with the correct solution.
//         *            In that case, the assemble function would be called twice with the same solution vector,
//         *            once at the end of the last iteration and a second time during the first iteration of
//         *            the nex timestep. The hysteresis operator could thus see, that the same input was used
//         *            and reuse the old values, whereas it would need to recompute it if the last evaluation
//         *            was performed with a temporal solution during the lineSearch algorithm.
//         *            In fact this first reason is already considered by the fact, that the update during
//         *            the lineSearch only use a temporal copy of the hysteresis operator to not destroy
//         *            the memory with temporal inputs.
//         *            A second reason is simply, that both the version with lineSearch and the one without
//         *            will have an updated system at the end.
//         */
//
//       // option b) if ( lineSearch_ != "none" || iterationCounter == 1) {
//        if ( iterationCounter == 1) {
//          /*
//           * Setup system for the initial iteration (otherwise we could not solve nothing)
//           */
//
//          // set linear part of RHS
//          algsys_->InitRHS(RhsLinVal_);
//
//          // set nonlinear part of RHS (if any)
//          /*
//           * this is not included in StepTransNonLin()
//           * -> does it not work or is it simply not needed there?
//           */
//          assemble_->AssembleNonLinRHS();
//
//          // setup the matrices
//          isNewton = false;
//          assemble_->AssembleMatrices(isNewton);
//
//          //now update RHS according to time stepping
//          for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//
//            if(matIt->second < 0)
//              continue;
//            for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//              fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
//            }
//            algsys_->UpdateRHS(matIt->first,stageRHS_,true);
//          }
//          //substract from RHS the term K*sol
//          solVec_.ScalarMult(-1.0);
//          algsys_->UpdateRHS(STIFFNESS,solVec_,true); // we also or only need the updated version
//          algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
//          solVec_.ScalarMult(-1.0);
//
//          /*
//           * Now we should have the system (neglecting M and C-Matrices):
//           * K(old_solution) * solution_increment = f_lin + f_nonlin(old_solution) - K(old_solution) * old_solution
//           * which is equal to
//           * K(old_solution) * (old_solution + solution_increment)  = f_lin + f_nonlin(old_solution)
//           * which is equal to
//           * K(old_solution) * new_solution = f_lin + f_nonlin(old_solution)
//           */
//        }
//
//        //now assemble the Newton bilinear forms
//        /*
//         * this is not needed for the hysteresis case, as we do not have the derivative
//         * of the Hysteresis operator;
//         * it might be of interest, if some regions have hysteresis and other regions have
//         * 'simple' nonlinearities, which allow a treatment with Newtons methods.
//         * The big question is, if such a combination will work at all.
//         */
////        isNewton = true;
////        assemble_->AssembleMatrices(isNewton);
//
//        matrix_factor_.clear();
//
//        // set system matrix to zero initially, as ConstructEffectiveMatrix only
//        // sums up the contributions
//        algsys_->InitMatrix(SYSTEM);
//        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
//          FeFctIdType fctId = fncIt->second->GetFctId();
//          fncIt->second->GetTimeScheme()
//            ->AddMatFactors(i,matrices,matrix_factor_[fctId]);
//          algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
//        }
//        // setup the matrices to compute correct error norms
//
//        PDE_.SetBCs();
//        algsys_->BuildInDirichlet();
//        algsys_->SetupPrecond();
//        algsys_->SetupSolver();
//
//        // just set inh. Dirichlet BCs for the first iteration
//        /*
//         * I do not fully understand this by now.
//         * Why do we set inh. DBC only during the first iteration?
//         */
//        bool setIDBC = false;
//        if ( iterationCounter == 1 && couplingIter_ == 0 )
//          setIDBC = true;
//
//        algsys_->Solve(setIDBC);
//        algsys_->GetSolutionVal(solInc, setIDBC );
//
//        //std::cout << solInc.ToString(1,',') << std::endl;
//
//        Double residualL2Norm = 0.0;
//        Double etaLineSearch  = 1.0;
//
//        //necessary due to inh. Dirichlet BCs!!
//        /*
//         * Why?
//         */
//        if ( iterationCounter == 1 && couplingIter_ == 0 )
//          stageSol.Init();
//
//        if ( lineSearch_ == "none" || iterationCounter == 1) {
//          //to incooperate the inhomog. Dirichlet BCs we need a full
//          //step for the first iteration
//          /*
//           * Why?
//           */
//
//          stageSol.Add(1.0, solInc);
//          solVec_  = stageSol;
//        } else {
//          /*
//           * lock hysteresis operator to ensure that the tested lineSearch steps do not
//           * disturb the memory of the operator
//           */
//          PDE_.LockHysteresis();
//
//          residualL2Norm = LineSearch(solInc, stageSol, etaLineSearch,true);
//          solVec_  = stageSol;
//          /*
//           * unlock afterwards
//           */
//
//          PDE_.UnlockHysteresis();
//        }
//
//        /*
//         * according to the previously stated option a) we update the algebraic system
//         * with the correct next solution at this point (also for the lineSearch case
//         * (in StepTransNonLin, the following lines follow only in the
//         * if ( lineSearch_ == "none" || iterationCounter == 1) {
//         * case
//         */
//
//        // setup the matrices with new solution
//        /*
//         * As stated previously, by setting solVec_ to the new solution (by updating the old solution
//         * with the solution increment), the results are also set in the corresponding feFunctions,
//         * such that assemble and algsys can access it.
//         */
//        isNewton = false;
//        assemble_->AssembleMatrices(isNewton);
//
//        // set linear part of RHS
//        algsys_->InitRHS(RhsLinVal_);
//        /*
//         * again, this is not included in StepTransNonLin()
//         * -> does it not work or is it simply not needed there?
//         */
//        assemble_->AssembleNonLinRHS();
//
//        //now update RHS according to time stepping
//        for(matIt = matrices.begin();matIt != matrices.end();matIt++){
//          //std::cout << "Matrix: " << matIt->first << std::endl;
//          if(matIt->second < 0)
//            continue;
//          for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//            fncIt->second->GetTimeScheme()->ComputeStageRHS(i,matIt->second,stageRHS_.GetPointer(pos));
//          }
//          algsys_->UpdateRHS(matIt->first,stageRHS_,true);
//        }
//
//        //substract from RHS the term K*sol
//        solVec_.ScalarMult(-1.0);
//        algsys_->UpdateRHS(STIFFNESS,solVec_,true);
//        algsys_->UpdateRHS(STIFFNESS_UPDATE,solVec_,true);
//        solVec_.ScalarMult(-1.0);
//
//        if ( lineSearch_ == "none" || iterationCounter == 1) {
//          /*
//           * this part is only for the case that we did not perform a linesearch
//           * in that case, we have to evaluate the residual which otherwise is
//           * done during linesearch
//           */
//          //get RHS vector
//          SBM_Vector actRHS(BaseMatrix::DOUBLE);
//          algsys_->GetRHSVal( actRHS );
//
//          // calculation of residual error =======================================
//          residualL2Norm = actRHS.NormL2();
//        }
////        else {
////           //solVec_  = stageSol;
////          residualL2Norm = LineSearch(solInc, stageSol, etaLineSearch,true);
////          solVec_  = stageSol;
////        }
//
//        // calculation of residual error =======================================
//        Double residualErr;
//        if ( RhsLinL2Norm > 1.0 )
//          residualErr = residualL2Norm / RhsLinL2Norm;
//        else
//          residualErr = residualL2Norm;
//
//        // calculate incremental error ========================================
//        Double incrementalErr;
//        Double solIncrL2Norm = solInc.NormL2();
//        Double actSolL2Norm  = stageSol.NormL2();
//
//        if ( actSolL2Norm)
//          incrementalErr = solIncrL2Norm / actSolL2Norm;
//        else {
//          incrementalErr = solIncrL2Norm;
//          //WARN("Zero solution vector!! ");
//        }
//
//        std::cout << "incrementalErr: " << incrementalErr << std::endl;
//        std::cout << "residualErr: " << residualErr << std::endl;
//
//
////        if(incrementalErr > oldIncError){
////          std::cout << "Error increased since last iteration; do a reset of deltaMat" << std::endl;
////          // FinalizeAfterTimestep sets old hysteresis values to 0 so that deltaP and deltaE is basically
////          // P_current/E_current
////          PDE_.FinalizeAfterTimeStep();
////        } else {
////          std::cout << "Error decreased since last iteration; continue as normal" << std::endl;
////          // FinalizeAfterIteration sets old hysteresis values so that deltaP and deltaE can be computed
////          PDE_.FinalizeAfterIteration();
////        }
////
////        oldIncError = incrementalErr;
////
//
//        // output of norms and data
//        if ( nonLinLogging_ == true ) {
//          // get current step
//          UInt actStep = PDE_.GetSolveStep()->GetActStep();
//
//          if (PDE_.IsIterCoupled()) {
//            WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
//          } else {
//                WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLineSearch);
//          }
//          // write norm to file
//          logFile_ <<  iterationCounter << "\t"
//              << residualErr << "\t"
//              << incrementalErr << "\t"
//              << etaLineSearch << std::endl;
//        }
//
//        // boolean variable, holds condition if another iteration step is necessary
//        performOneMoreStep =
//          (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
//
//       // if (performOneMoreStep){
//          PDE_.FinalizeAfterIteration();
//       // }
//
//        if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
//          EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
//                  << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
//                  << "' at iteration '" << iterationCounter
//                  << "'.\n ==> incremental error: " << incrementalErr
//                  << "\n ==> residual error: " << residualErr);
//        }
//      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);
//
//    } //stages
//
//    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
//    limitFeFctIt = feFunctions_.find(solutionLimit_);
//    if (limitFeFctIt != feFunctions_.end() ) {
//      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
//        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
//          //const SingleVector * subv = solVec_.GetPointer(pos);
//          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
//          for (UInt j=0; j < dsubVec.GetSize(); j++) {
//            if (dsubVec[j] >= maxValidValue_) {
//              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
//                      "' is larger than the allowed maximum limit set in the XML: "
//                      << maxValidValue_);
//            }
//            if (dsubVec[j] <= minValidValue_) {
//              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
//                      "' is smaller than the allowed minimum limit set in the XML: "
//                      << minValidValue_);
//            }
//          }
//        }
//      }
//    }
//
//    //update stage
//    for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
//      fncIt->second->GetTimeScheme()->FinishStep();
//    }
//
//    //PDE_.FinalizeAfterTimeStep();
//    //PDE_.FinalizeAfterIteration();
//  }

  void StdSolveStep::StepTransNonLinHysteresisTotal() {

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
       if(PDE_.IsHysteresis_Fixpoint() == true){

         // nach dem löschen muss er eigentlich sowohl die linearen als auch die nichtlinearen teile nochmals hinzufügen
         // die unterscheidung zwischen linear und nichtlinear ist also obsolet

         RhsLinVal_.Init();
         algsys_->InitRHS(RhsLinVal_);

         std::cout << "SetLinRHS" << std::endl;
         incrementalErr = SetLinRHS(loadFactor,true);
         std::cout << "Done with setting RHS" << std::endl;
         /*
          * Fixpoint iteration used to work, however, P was just postprocessed upon E
          * as we RhsLinVal_.Init() followed by SetLinRHS( ,true) we just had nonlinear terms on the rhs
          * although P depends on E it was not marked as solDependent and such we solved div(eps0*E) = 0
          * this worked of course
          * Doing a real updating on the rhs will however not work
          */
         //incrementalErr = SetLinRHS(loadFactor,false);
       } else {
         //incrementalErr = SetLinRHS(loadFactor,false);
       }

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

       //just dummy things
       Double etaLineSearch = 1.0;
       Double residualErr = 1.0;
       // why not just try the linesearch algorithm?
       // actual linesearch is meant for the case of delta formulation
       //   i.e. we compute deltaSol using the equation system
       //   K(oldSol) * deltaSol = f - K(oldSol)*oldSol
       //   then we find out eta such that residual
       //   f - K(oldSol + eta*deltaSol)*(oldSol + eta*deltaSol) is minimal
       //   finally newSol = oldSol + eta*deltaSol
       // in our case, we compute newSol directly via
       //   K(oldSol)*newSol = f
       //   we can however try the following here
       //     compute diffSol = newSol - oldSol
       //     use linesearch to find a better eta than 1.0
       //     set newSol = oldSol + eta*diffSol
       // Note: oldSol will be overwritten!
       stageSol = oldSol;

       // before performing linesearch, make hysteresis operator read only,
       // i.e. do not store the effects of temporal input vectors (as they are
       // used by the linsearch function to test which eta is the best)
       //PDE_.LockHysteresis();

       residualErr = LineSearch(diffSol,stageSol,etaLineSearch,true);
       std::cout << "Used etaLineSearch: " << etaLineSearch << std::endl;
       std::cout << "Old solution: " << oldSol.ToString() << std::endl;
       std::cout << "New solution: " << stageSol.ToString() << std::endl;

       //PDE_.UnlockHysteresis();
       // Why is residualErr = incrementalErr?
       //Double residualErr = incrementalErr;

       // calculate new incremental error ========================================
       diffSol = stageSol; // this is the solution we just calculated via linesearch
       diffSol.Add( -1.0, oldSol);
       solIncrL2Norm = diffSol.NormL2();
       solNewL2Norm = newSol.NormL2();

       if (solNewL2Norm > 1)
         incrementalErr = solIncrL2Norm / solNewL2Norm;
       else
         incrementalErr = solIncrL2Norm;

       std::cout << "Incremental error: " << std::endl;
       std::cout << "rel error: " << incrementalErr << std::endl;

       std::cout << "Residual error: " << std::endl;
       std::cout << "rel error: " << residualErr << std::endl;

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
       //stageSol = newSol;
       //Double relaxfac = 0.1;
       //stageSol.Add( -relaxfac, newSol);
       //stageSol.Add( relaxfac, oldSol);
       //stageSol = newSol;
       solVec_  = stageSol;

       //store new solution for next iteration
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


    //Set special RHS Values
    //std::cout << "Do Apply Loads" << std::endl;
    PDE_.SetRhsValues();

    //this has to be done each frequency!
    assemble_->AssembleLinRHS();

    assemble_->AssembleMatrices( );
    PDE_.SetBCs();

    // store rhs vector back to PDE
    algsys_->GetRHSVal( rhsVec_ );

    // Where should we get the matrix factors from in a harmonic case?
    // In my opinion this method
    //if( assemble_->IsMatrixUpdated() ) {
    std::map<FEMatrixType,Double> empty;
    algsys_->ConstructEffectiveMatrix(NO_FCT_ID,  empty );

    // Check if the AMG-framework is used (if so, we have
    // to gather some geometry information at this point)
    // needs only be built once, doesn't change over frequency
    if( ((algsys_->UseAMG()) && (auxSet_ == false)) ||
        ((algsys_->UseAMG()) && (algsys_->GetAMGType() != AMGType::EDGE)) ){
      // only works for elimination
      if( solStrat_->UseDirichletPenalty() ) EXCEPTION("AMG only works for Dirichlet elimination!");
      PDE_.SetGeomInfo();
      algsys_->BuildAMGAuxMatrix();
      auxSet_ = true;
    }

    // Incorporate Boundary conditions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
    }

    algsys_->Solve();
    algsys_->GetSolutionVal(solVec_);

    if ( adjointSource_ ) {
    	//check if adjoint PDE has been solved in case of source localization
    	//if yes, we have to multiply the solution with a standard mass matrix
    	std::cout << "DO multiply with MASS-matrix" << std::endl;
    	//solVec_.Export("sol1.dat",BaseMatrix::MATRIX_MARKET);
    	algsys_->InitRHS();
    	algsys_->UpdateRHS(AUXILIARY,solVec_,true);
    	algsys_->GetRHSVal( solVec_ );
    	//solVec_.Export("sol2.dat",BaseMatrix::MATRIX_MARKET);
    	//std::cout << "SOL after: \n " << solVec_ << std::endl;
    	adjointSource_ = false;
    }
  }


  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds,
                                           UInt numFreq, double shift, bool sort) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver(numFreq, shift, false, sort, false);

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies(frequencies, errBounds);

    return frequencies.GetSize();
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex>& frequencies, Vector<Double>& errBounds,
                                           UInt numFreq, Double shift, bool sort, bool bloch) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver  - we cannot be quadratic and bloch concurrently!
    algsys_->SetupEigenSolver(numFreq, shift, !bloch, sort, bloch);

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

    //std::cout << "SetLinRHS with bool nonlin = " << nonlin << std::endl;
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
      //std::cout << "Testing eta = " << eta[i] << std::endl;
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

    //std::cout << "Optimal eta = " << etaOpt << std::endl;
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
