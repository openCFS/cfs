#include <fstream>
#include <iostream>
#include <string>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "SolveStepEB.hh"

#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Utils/Timer.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"

namespace CoupledField
{

  DEFINE_LOG(solvestepeb, "solvestepeb")

  SolveStepEB::SolveStepEB(StdPDE &apde) : StdSolveStep(apde)
  {
    matModelCoef_ = apde.GetModelCoef();
  }

  SolveStepEB::~SolveStepEB()
  {
  }

  // ======================================================
  // Solve Step Transient SECTION
  // ======================================================

  void SolveStepEB::SolveStepTrans()
  {
    if (nonLin_)
    {
      StepTransNonLin();
    }
    else
    {
      StepTransLin();
    }
  }

  void SolveStepEB::StepTransNonLin()
  {
 
    mParser_->SetExpr(MathParser::GLOB_HANDLER,"iterationCounter");
    mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", 0);

    LOG_TRACE(solvestepeb) << "SolveStepEB::StepTransNonLin";
    bool performOneMoreStep;
    bool isNewton = false;
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    SBM_Vector stageSol(BaseMatrix::DOUBLE);
    SBM_Vector actRHS(BaseMatrix::DOUBLE);
    SBM_Vector actSol(BaseMatrix::DOUBLE);
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;
    
    for(UInt i=0;i<numStages;i++){
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        FeFctIdType fncId = fncIt->second->GetFctId();

        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),fncId);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false);

/*       //do initialization 
      rhsVec_.Init();
      LOG_DBG(solvestepeb) << "StepTransNonLin: Stage: " << i ;
      
      // setup right hand side
      algsys_->InitRHS();
      assemble_->AssembleLinRHS(); */

      // set iteration counter
      UInt iterationCounter=0;

/*       LOG_DBG3(solvestepeb) << "\n===========================================================================\n"
                             << "============================ TIMESTEP: "<< PDE_.GetSolveStep()->GetActStep()<<" ==========================\n"
                             << "============================================================================"; */


      iterationCounter++;
      mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);
                
      // setup the matrices, which automatically evaluates the EBHysteresis model
      // in the EBHysteresis class, we have a switch to detect the iterationstep=1
      // and use the values from the previous timestep
      if(PDE_.GetSolveStep()->GetActStep() == 1){
        isNewton = false;
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
      }

/*       PDE_.SetBCs();
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond();
      algsys_->SetupSolver();

      // just set inh. Dirichlet BCs for the first iteration
      bool setIDBC = true;        
      algsys_->Solve(setIDBC);

      // if setIDBC is true, solInc will contain the inhom. Dirichlet values
      // Since the entries of solVec_ are pointers to the SingleVector
      // of the FE function, it automatically inserts the values there
      algsys_->GetSolutionVal(solInc, setIDBC );
      
      stageSol.Init();
      stageSol.Add(1.0, solInc); */

      solVec_  = stageSol;
      
      LOG_DBG3(solvestepeb) << "================ Vectors after stageSol.Add(1.0, solInc): ======================="
            << "\n ===== iterationCounter:"<< iterationCounter
            << "\n solInc:"<< solInc.ToString()
            << "\n actSol:"<< actSol.ToString()
            << "\n solVec_:"<< solVec_.ToString()
            << "\n stageSol:"<< stageSol.ToString()
            << "\n rhsVec_:"<< rhsVec_.ToString()
            << "\n RhsLinVal_:"<< RhsLinVal_.ToString()
            << "\n stageRHS_:"<< stageRHS_.ToString();




      // ===================================================================================
      // =================== START NONLINEAR ITERATION =====================================
      // ===================================================================================
      do {
        iterationCounter++;
        mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);
                  
        algsys_->InitRHS();
        // if the RHS depends on the nonlinearity, we have to re-assemble it
        assemble_->AssembleLinRHS();
        assemble_->AssembleNonLinRHS();

        // setup the matrices, which automatically evaluates the EBHysteresis model
        // since this is iterationstep > 1, we actually perform the material tensor evaluation
        isNewton = false;
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

        algsys_->GetRHSVal( actRHS );
        LOG_DBG3(solvestepeb) << "rhsVec_:"<< rhsVec_.ToString()
                              << "\n RhsLinVal_:"<< RhsLinVal_.ToString()
                              << "\n stageRHS_:"<< stageRHS_.ToString()
                              << "\n actRHS:"<< actRHS.ToString();

        PDE_.SetBCs();
        algsys_->BuildInDirichlet();
        algsys_->SetupPrecond();
        algsys_->SetupSolver();
        bool setIDBC = true;
        algsys_->Solve(setIDBC);


        // if setIDBC is true, solInc will contain the inhom. Dirichlet values
        // Since the entries of solVec_ are pointers to the SingleVector
        // of the FE function, it automatically inserts the values there
        algsys_->GetSolutionVal(solInc, setIDBC );

        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none"){
          stageSol.Add(etaLineSearch, solInc);
        }else if ( lineSearch_ == "minEnergy"){
          LineSearchHeavy(solInc, stageSol, etaLineSearch);
        }else if ( lineSearch_ == "Armijo"){
          //LineSearchArmijo(solInc, stageSol, etaLineSearch, iterationCounter);
          LineSearchArmijo(solInc, stageSol, etaLineSearch);
        }else if ( lineSearch_ == "ArmijoRegularization"){
          LineSearchArmijoRegularization(solInc, stageSol, etaLineSearch, iterationCounter);
        }

        solVec_  = stageSol;

        LOG_DBG3(solvestepeb) << "================ Vectors after stageSol.Add(1.0, solInc): ======================="
            << "\n ===== iterationCounter:"<< iterationCounter
            << "\n solInc:"<< solInc.ToString()
            << "\n actSol:"<< actSol.ToString()
            << "\n solVec_:"<< solVec_.ToString()
            << "\n stageSol:"<< stageSol.ToString()
            << "\n rhsVec_:"<< rhsVec_.ToString()
            << "\n RhsLinVal_:"<< RhsLinVal_.ToString()
            << "\n stageRHS_:"<< stageRHS_.ToString();

        
        // ======= residual computation ============
        // start with assembling the nonlinear rhs, which is div(B)
        // set nonlinear rhs
        algsys_->InitRHS();
        assemble_->AssembleLinRHS();
        // if the RHS depends on the nonlinearity, we have to re-assemble it
        assemble_->AssembleNonLinRHS();

        algsys_->GetRHSVal( actRHS );
        LOG_DBG3(solvestepeb) << "rhsVec_:"<< rhsVec_.ToString()
                              << "\n RhsLinVal_:"<< RhsLinVal_.ToString()
                              << "\n stageRHS_:"<< stageRHS_.ToString()
                              << "\n actRHS:"<< actRHS.ToString();
        
        
        // calculation of residual error =======================================
        Double residualL2Norm = 0.0;
        residualL2Norm = actRHS.NormL2();
        Double residualErr = residualL2Norm;
        
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

        OutputNonLinIterInfo(pdename_, PDE_.GetSolveStep()->GetActStep(),iterationCounter, residualErr, incrementalErr, etaLineSearch, PDE_.IsIterCoupled() ? couplingIter_ : -1);


        if ( nonLinLogging_ == true ) { 
          // write norm to file
          logFile_ <<  iterationCounter << "\t"
                  << residualErr << "\t"
                  << incrementalErr << "\t"
                  << etaLineSearch << std::endl;
        }
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
        
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
        FeFctIdType fncId = fncIt->second->GetFctId();
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(fncId)));
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



  void SolveStepEB::LineSearchArmijoRegularization(SBM_Vector& solIncrement, SBM_Vector& actSol,
  Double& gamma, UInt iterationCounter)  {
    
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    gamma = 1.0;
    Double sigma1 = 1.0e-4;
    Double sigma2 = 1.0e-4;
    Double rho = 0.01;
    Double eta = 1.1;
    Double phi = 2.0;
    Integer i = (-1.0) * (iterationCounter-1);
    Double eta_k = std::pow(eta, i);
    UInt numLSIter = 0;
    Double res_x_trial = 0.0;

    Double solIncNorm = solIncrement.NormL2();
    algsys_->InitRHS();
    assemble_->AssembleLinRHS();
    assemble_->AssembleNonLinRHS();
    SBM_Vector actRHS(BaseMatrix::DOUBLE);
    algsys_->GetRHSVal( actRHS );
    Double startingResidual = actRHS.NormL2();

    do{
      actSol.Add( 1.0, solOld, gamma, solIncrement);
      //store new solution
      solVec_ = actSol;

      // ======= residual computation ============
      algsys_->InitRHS();
      assemble_->AssembleLinRHS();
      assemble_->AssembleNonLinRHS();
      algsys_->GetRHSVal( actRHS );
      // calculation of residual error =======================================
      res_x_trial = actRHS.NormL2();

      if( res_x_trial <= (rho * startingResidual - sigma2 * std::pow(startingResidual,2)) ){
        break;
      }else{
        Double t = gamma * solIncNorm;
        if( res_x_trial <= startingResidual - sigma1 * std::pow(t,2) + eta_k* startingResidual){
        break;
        }else if( gamma < 1.0e-16){
          break;
        }else{
          gamma = 1.0/std::pow(phi,numLSIter);
        }
      }
      numLSIter++;
    }while(numLSIter < 100); 
  }



  void SolveStepEB::LineSearchArmijo(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch)   {

    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 11;
    const Double eta[nrEtas] = {1,0.6667,0.4444,0.2963,0.1975,0.1317,0.0878,0.0585,0.039,0.026,0.01733333333};
    //const Double eta[nrEtas] = {1,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001};
    const Double eta0 = 0;

    Double delta = 1e-8;
    
    // initialize etaOpt
    Double etaOpt = 0.0;

    // obtain reference residual (start)
    actSol.Add( 1.0, solOld, eta0, solIncrement);
    solVec_ = actSol;
    algsys_->InitRHS();
    assemble_->AssembleLinRHS();
    assemble_->AssembleNonLinRHS();
    SBM_Vector actRHS(BaseMatrix::DOUBLE);
    algsys_->GetRHSVal( actRHS );
    Double residualL2Norm_eta0 = actRHS.NormL2();
    // obtain reference residual (end)
    
    for( UInt i=0; i<nrEtas; i++) {

      // obtain trial residual (start)
      actSol.Add( 1.0, solOld, eta[i], solIncrement);
      solVec_ = actSol;
      algsys_->InitRHS();
      assemble_->AssembleLinRHS();
      assemble_->AssembleNonLinRHS();
      SBM_Vector actRHS(BaseMatrix::DOUBLE);
      algsys_->GetRHSVal( actRHS );
      Double residualL2Norm = actRHS.NormL2();
      // obtain trial residual (end)

      // check if the trial residual has a sufficient decrease
      if (residualL2Norm < (1.0 - delta*eta[i])*residualL2Norm_eta0) {
        //printf("residualL2Norm: %f; residualL2Norm_eta0: %f\n",residualL2Norm,residualL2Norm_eta0);
        etaOpt = eta[i];
        break;
      }
      // after all possible line search parameter are tried take the last one and break the for loop
      if (i == nrEtas-1){
        etaOpt = eta[i];
        break;
      }
    }
    etaLineSearch = etaOpt;
    // Set new solution
    actSol.Add( 1.0, solOld, etaOpt, solIncrement );
    
  }


  void SolveStepEB::LineSearchHeavy(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch)  {
    
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 10;
    //const Double eta[nrEtas] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    const Double eta[nrEtas] = {0.026,0.039,0.0585,0.0878,0.1317,0.1975,0.2963,0.4444,0.6667,1};
    printf("eta[0]: %f\n",eta[0]);
    
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;


    for( UInt i=0; i<nrEtas; i++) {
      actSol.Add( 1.0, solOld, eta[i], solIncrement);
      //store new solution
      solVec_ = actSol;

      // ======= residual computation ============
      // start with assembling the nonlinear rhs, which is div(B)
      // set nonlinear rhs
      algsys_->InitRHS();
      assemble_->AssembleLinRHS();
      // if the RHS depends on the nonlinearity, we have to re-assemble it
      assemble_->AssembleNonLinRHS();

      // calculation of error norms
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
  }


} // end of namespace
