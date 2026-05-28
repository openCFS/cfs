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
    matModelCoefm_ = apde.GetModelCoefm();
    pseudo_time_stepping_ = 1;
  }
  SolveStepEB::SolveStepEB(StdPDE &apde, UInt is_pseudo_time_stepping) : StdSolveStep(apde)
  {
    matModelCoefm_ = apde.GetModelCoefm();
    pseudo_time_stepping_ = is_pseudo_time_stepping;
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

  void SolveStepEB::StepTransNonLin() {
 
    mParser_->SetExpr(MathParser::GLOB_HANDLER,"iterationCounter");
    mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", 0);

    bool performOneMoreStep;
    bool isNewton = false;
    
    SBM_Vector solInc(BaseMatrix::DOUBLE);
    SBM_Vector stageSol(BaseMatrix::DOUBLE);
    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE);
    SBM_Vector actRHS(BaseMatrix::DOUBLE);
    SBM_Vector Linform_nm1(BaseMatrix::DOUBLE);
    SBM_Vector Linform_nm1_temp(BaseMatrix::DOUBLE);
    
    //obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    std::map<FEMatrixType,Integer>::iterator matIt;
    
    UInt pos = 0;

    LOG_DBG2(solvestepeb) << "STARTING SetTransNonLin() =============================================";
    LOG_DBG2(solvestepeb) << "numStages:" << numStages;
    LOG_DBG2(solvestepeb) << "actTime_:" << actTime_;

    for(UInt i=0;i<numStages;i++){
      stageSol.Resize(feFunctions_.size());
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        FeFctIdType fncId = fncIt->second->GetFctId();

        stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(i),fncId);
        fncIt->second->GetTimeScheme()->InitStage(i,actTime_,PDE_.GetDomain());
      }
      stageSol.SetOwnership(false); 

      LOG_DBG2(solvestepeb) << "stageSol before while:" << stageSol.ToString();

      // set iteration counter
      UInt iterationCounter=0;

      Double residualErr = 0.0;
      Double residualErr0 = 0.0;
      Double incrementalErr = 0.0;
      Double solIncrL2Norm = 0.0;
      Double stageSolL2Norm  = 0.0;
      // ===================================================================================
      // =================== START NONLINEAR ITERATION =====================================
      // ===================================================================================
      std::cout << "iteration" << "     " << "residual 2-norm" << "     " <<"step 2-norm" << "           " <<"    eta    " << std::endl;
      std::cout << "---------" << "     " << "---------------" << "     " <<"-----------" << "           " <<"-----------" << std::endl;
      // ###############################################################
      // PSEUDO TIME STEPPING (START)
      // ###############################################################
      if (pseudo_time_stepping_ == 1) {

        Double residualErrPrev = 1e15;
        const UInt windowSize = stagnationWindowSize_;
        std::vector<Double> improvementHistory;
        const Double stagnationTol = stagnationTolerance_;

        bool useStagnationDetection = useStagnationDetection_;
        while (true){
          iterationCounter++; mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);
          stageSol_temp = stageSol;

          LOG_DBG2(solvestepeb) << "=============== Start iteration " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();

          // set up RHS
          algsys_->InitRHS();
          assemble_->AssembleLinRHS();
          assemble_->AssembleNonLinRHS();
          algsys_->GetRHSVal( actRHS );
          if(iterationCounter == 1){
            residualErr0 = actRHS.NormL2();
            std::cout <<"residual error @ iteration 0: " << residualErr0 << std::endl;
            if ( residualErr0 < 1.0 )
              residualErr0 = 1.0;
          }



          LOG_DBG2(solvestepeb) << "\n\t\t =============== after setup RHS " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();


          // set up matrix
          assemble_->AssembleMatrices(isNewton);
          matrix_factor_.clear();
          algsys_->InitMatrix(SYSTEM);
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
            algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
          }

          LOG_DBG2(solvestepeb) << "\n\t\t =============== after setup SYS matrix " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();


          // solve system
          PDE_.SetBCs();
          algsys_->BuildInDirichlet();
          algsys_->SetupPrecond();
          algsys_->SetupSolver();
          bool setIDBC = true;
          algsys_->Solve(setIDBC);
          algsys_->GetSolutionVal(solInc, setIDBC );

          LOG_DBG2(solvestepeb) << "\n\t\t =============== after SOLVE " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();

          // apply line search
          Double etaLineSearch = 1.0;
          if ( lineSearch_ == "none"){
            stageSol.Add(etaLineSearch, solInc);
          }else if ( lineSearch_ == "minEnergy"){
            etaLineSearch = ExactLineSearch(solInc, stageSol);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }else if ( lineSearch_ == "minEnergyHeavy"){
            LineSearchHeavy(solInc, stageSol, etaLineSearch);
          }else if ( lineSearch_ == "Armijo"){
            etaLineSearch = LineSearchArmijo(solInc, stageSol);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }else if ( lineSearch_ == "Inexact"){
            etaLineSearch = InexactLineSearch(solInc, stageSol);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }

          LOG_DBG2(solvestepeb) << "\n\t\t =============== after LINESEARCH " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t etaLineSearch " << etaLineSearch;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();
        
          std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
          for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->AllowUpdates(true);  
          } 
          
          // residual
          algsys_->InitRHS();
          assemble_->AssembleLinRHS();
          assemble_->AssembleNonLinRHS();
          algsys_->GetRHSVal( actRHS );
          residualErr = actRHS.NormL2();
          residualErr = std::abs(residualErr)/residualErr0;

          Double maxImprovement = 1.0; // default: "not stagnating"

          if (useStagnationDetection) {
              Double improvement = (residualErrPrev - residualErr)
                                    / (residualErrPrev + 1e-15);

                improvementHistory.push_back(improvement);

                if (improvementHistory.size() > windowSize) {
                    improvementHistory.erase(improvementHistory.begin());
                }

                maxImprovement = *std::max_element(
                    improvementHistory.begin(),
                    improvementHistory.end()
                );
          }

          residualErrPrev = residualErr;
          LOG_DBG2(solvestepeb) << "\n\t\t =============== after RESIDUAL " << iterationCounter;
          LOG_DBG2(solvestepeb) << "\n\t\t solInc:" << solInc.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol:" << stageSol.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t stageSol_temp:" << stageSol_temp.ToString();
          LOG_DBG2(solvestepeb) << "\n\t\t actRHS:" << actRHS.ToString();

          
          // calculate incremental error ========================================
          solIncrL2Norm = solInc.NormL2();
          stageSolL2Norm  = stageSol.NormL2();
          if ( stageSolL2Norm )
            incrementalErr = solIncrL2Norm/stageSolL2Norm;
          else {
            incrementalErr = solIncrL2Norm;
          }
          std::cout <<"    " << iterationCounter << "           " << residualErr << "       " << incrementalErr << "       " << etaLineSearch <<"\n" << std::scientific;
          OutputNonLinIterInfo(pdename_, PDE_.GetSolveStep()->GetActStep(),iterationCounter, residualErr, incrementalErr, etaLineSearch, PDE_.IsIterCoupled() ? couplingIter_ : -1);

          // boolean variable, holds condition if another iteration step is necessary
          //performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

          bool residualOk    = (residualErr    <= residualStopCrit_);
          bool incrementalOk = (incrementalErr <= incStopCrit_);
          
          bool stagnating = (maxImprovement < stagnationTol)
                          && (improvementHistory.size() == windowSize)
                          && (incrementalErr <= incStopCrit_);

          performOneMoreStep = !( (residualOk && incrementalOk) || stagnating );

          // normal convergence
          if ( performOneMoreStep == 0 && !stagnating ) {
            std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
            for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->UpdateHistoryValues(); 
              it->second->AllowUpdates(false);
            }
            break;
          }

          // stagnation: warning + continue with next time step
          if ( stagnating ) {
            std::cout << "\n***********************************************************************\n"
                      << " WARNING: Stagnation detected in PDE '" << pdename_
                      << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                      << "' at iteration '" << iterationCounter
                      << "'.\n ==> incremental error: " << incrementalErr
                      << "\n ==> residual error: "    << residualErr
                      << "\n ==> max improvement over last " << windowSize << " iterations: " 
                      << maxImprovement
                      << "\n continuing with next time step..."
                      << "\n***********************************************************************\n";
            std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
            for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->UpdateHistoryValues(); 
              it->second->AllowUpdates(false);
            }
            break;
          }
          
          // hard limit - real error
          if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
            EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                    << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                    << "' at iteration '" << iterationCounter 
                    << "'.\n ==> incremental error: " << incrementalErr
                    << "\n ==> residual error: " << residualErr);
          }
        }
      }
      // ###############################################################
      // PSEUDO TIME STEPPING (END)
      // ###############################################################
      // ###############################################################
      // REAL TIME STEPPING (START)
      // ###############################################################
      if (pseudo_time_stepping_ == 0) {
        Double residualErrPrev = 1e15;
        const UInt windowSize = stagnationWindowSize_;
        std::vector<Double> improvementHistory;
        const Double stagnationTol = stagnationTolerance_;
        bool useStagnationDetection = useStagnationDetection_;

        while (true){
          iterationCounter++; mParser_->SetValue(MathParser::GLOB_HANDLER, "iterationCounter", iterationCounter);
          stageSol_temp = stageSol;

          // set up RHS in two steps:
          // 1) if Newton iteration counter is 1, then assemble the linear forms having quantities from the last time step (declared LINEAR)
          // 2) assemble all linear forms that only contain values from the last NEWTON iteration (declared NONLINEAR)
          if (iterationCounter == 1) {
            algsys_->InitRHS();
            assemble_->AssembleLinRHS();
            algsys_->GetRHSVal(Linform_nm1);
            assemble_->AssembleNonLinRHS();
          } else {
            algsys_->InitRHS(Linform_nm1);
            assemble_->AssembleNonLinRHS();
          }
          if(iterationCounter == 1){
            residualErr0 = actRHS.NormL2();
            std::cout <<"residual error @ iteration 0: " << residualErr0 << std::endl;
            if ( residualErr0 < 1.0 )
              residualErr0 = 1.0;
          }

          // set up matrix
          assemble_->AssembleMatrices(isNewton);
          matrix_factor_.clear();
          algsys_->InitMatrix(SYSTEM);
          for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
            FeFctIdType fctId = fncIt->second->GetFctId();
            fncIt->second->GetTimeScheme()->AddMatFactors(i,matrices,matrix_factor_[fctId]);
            algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
          }

          // solve system
          PDE_.SetBCs();
          algsys_->BuildInDirichlet();
          algsys_->SetupPrecond();
          algsys_->SetupSolver();
          bool setIDBC = true;
          algsys_->Solve(setIDBC);
          algsys_->GetSolutionVal(solInc, setIDBC );
          string nr_unknowns_string = solInc.ToString();
          UInt nr_unknowns = std::count(nr_unknowns_string.begin(), nr_unknowns_string.end(), ',');
          // The number of vector entries is the number of commas + 1
          nr_unknowns = nr_unknowns + 1;
          //std::cout << "nr_unknowns: " << nr_unknowns << "\n";

          // apply line search
          Double etaLineSearch = 1.0;
          if ( lineSearch_ == "none"){
            stageSol.Add(etaLineSearch, solInc);
          }else if ( lineSearch_ == "minEnergy"){
            Linform_nm1_temp = Linform_nm1;
            etaLineSearch = ExactLineSearch(solInc, stageSol, Linform_nm1_temp);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }else if ( lineSearch_ == "minEnergyHeavy"){
            Linform_nm1_temp = Linform_nm1;
            LineSearchHeavy(solInc, stageSol, etaLineSearch, Linform_nm1_temp);
          }else if ( lineSearch_ == "Armijo"){
            etaLineSearch = LineSearchArmijo(solInc, stageSol);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }else if ( lineSearch_ == "Inexact"){
            etaLineSearch = InexactLineSearch(solInc, stageSol);
            stageSol_temp.Add(etaLineSearch, solInc);
            stageSol = stageSol_temp;
          }

          std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
          for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->AllowUpdates(true);  
          } 

          // residual
          algsys_->InitRHS(Linform_nm1);
          //assemble_->AssembleLinRHS();
          assemble_->AssembleNonLinRHS();
          algsys_->GetRHSVal( actRHS );
          residualErr = actRHS.NormL2();
          //residualErr = std::abs(residualErr)/(residualErr0*nr_unknowns);
          residualErr = std::abs(residualErr)/(residualErr0);

          Double maxImprovement = 1.0; // default: "not stagnating"

          if (useStagnationDetection) {
              Double improvement = (residualErrPrev - residualErr)
                                    / (residualErrPrev + 1e-15);

                improvementHistory.push_back(improvement);

                if (improvementHistory.size() > windowSize) {
                    improvementHistory.erase(improvementHistory.begin());
                }

                maxImprovement = *std::max_element(
                    improvementHistory.begin(),
                    improvementHistory.end()
                );
          }

          residualErrPrev = residualErr;

          // calculate incremental error ========================================
          solIncrL2Norm = solInc.NormL2();
          stageSolL2Norm  = stageSol.NormL2();
          if ( stageSolL2Norm )
            incrementalErr = solIncrL2Norm/stageSolL2Norm;
          else {
            incrementalErr = solIncrL2Norm;
          }
          std::cout <<"    " << iterationCounter << "           " << residualErr << "       " << incrementalErr << "       " << etaLineSearch <<"\n" << std::scientific;
          OutputNonLinIterInfo(pdename_, PDE_.GetSolveStep()->GetActStep(),iterationCounter, residualErr, incrementalErr, etaLineSearch, PDE_.IsIterCoupled() ? couplingIter_ : -1);

          // boolean variable, holds condition if another iteration step is necessary
          bool residualOk    = (residualErr    <= residualStopCrit_);
          bool incrementalOk = (incrementalErr <= incStopCrit_);
          
          bool stagnating = (maxImprovement < stagnationTol)
                          && (improvementHistory.size() == windowSize)
                          && (incrementalErr <= incStopCrit_);

          performOneMoreStep = !( (residualOk && incrementalOk) || stagnating );

          // normal convergence
          if ( performOneMoreStep == 0 && !stagnating ) {
            // now that we have reached our convergence threshold for this timestep, let's save the states in our model
            std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
            for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->AllowUpdates(true);  
              assemble_->AssembleLinRHS();
              it->second->UpdateHistoryValues(); 
              it->second->AllowUpdates(false);
            } 
            break;
          }

          // stagnation: warning + continue with next time step
          if ( stagnating ) {
            std::cout << "\n***********************************************************************\n"
                      << " WARNING: Stagnation detected in PDE '" << pdename_
                      << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                      << "' at iteration '" << iterationCounter
                      << "'.\n ==> incremental error: " << incrementalErr
                      << "\n ==> residual error: "    << residualErr
                      << "\n ==> max improvement over last " << windowSize << " iterations: " 
                      << maxImprovement
                      << "\n continuing with next time step..."
                      << "\n***********************************************************************\n";
            std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >::iterator it;
            for (it = matModelCoefm_.begin(); it != matModelCoefm_.end(); it++) {
              it->second->UpdateHistoryValues(); 
              it->second->AllowUpdates(false);
            }
            break;
          }
          
          if (performOneMoreStep && iterationCounter == nonLinMaxIter_ && abortOnMaxIter_) {
            EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_ 
                    << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
                    << "' at iteration '" << iterationCounter 
                    << "'.\n ==> incremental error: " << incrementalErr
                    << "\n ==> residual error: " << residualErr);
          }
        }
      }
      // ###############################################################
      // REAL TIME STEPPING (END)
      // ###############################################################
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
      fncIt->second->GetTimeScheme()->FinishStep();
    }
  }

  double SolveStepEB::ExactLineSearch(SBM_Vector& solIncrement, SBM_Vector& stageSol){

    Double a = 0.0 + 1e-3; //bottom interval of eta
    Double b = 1.0; //top interval of eta
    Double gamma = 0.0;

    Double Fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a);
    Double Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b);
    
    //for the case that there is no sign change in the line search interval (0,1], we fall back to an inexact (i.e., heuristic) line search method.
    //sign change is crucial bc otherwise there might be no root, a double root or an even number of roots in (0,1]
    //behvior of Brent's method is simply not defined if there is no sign change (see also https://en.wikipedia.org/wiki/Brent's_method)
    if (Fa * Fb >= 0.0) {
        std::cout << "WARNING: No sign change of the energy functional in the interval [" << a << ", " << b
                  << "] - falling back to InexactLineSearch\n";
        gamma = InexactLineSearch(solIncrement, stageSol); //note that Inexact is computationally more expensive
    } else {
        gamma = BrentMethod(solIncrement, stageSol, a, b);
    }
    return gamma;
  }

  double SolveStepEB::ExactLineSearch(SBM_Vector& solIncrement, SBM_Vector& stageSol, SBM_Vector& Linform_nm1){

    Double a = 0.0 + 1e-3; //bottom interval of eta
    Double b = 1.0; //top interval of eta
    Double gamma = 0.0;

    Double Fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a, Linform_nm1);
    Double Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b, Linform_nm1);
    
    if (Fa * Fb >= 0.0) {
        std::cout << "WARNING: No sign change of the energy functional in the interval [" << a << ", " << b
                  << "] - falling back to InexactLineSearch\n";
        gamma = InexactLineSearch(solIncrement, stageSol, Linform_nm1);
    } else {
        gamma = BrentMethod(solIncrement, stageSol, a, b, Linform_nm1);
    }
    return gamma;
  }

  double SolveStepEB::GetLineSearchDerivativeFunctionValue(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double eta){
    
    SBM_Vector residual_vector(BaseMatrix::DOUBLE);
    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE); stageSol_temp = stageSol;
    double EnergyDerivative = 0.0;

    stageSol.Add(eta,solIncrement);
    algsys_->InitRHS();
    assemble_->AssembleLinRHS();
    assemble_->AssembleNonLinRHS();
    algsys_->GetRHSVal( residual_vector );
    residual_vector.Inner(solIncrement,EnergyDerivative);
    stageSol = stageSol_temp;
    return EnergyDerivative;
  }

  double SolveStepEB::GetLineSearchDerivativeFunctionValue(SBM_Vector& solIncrement, SBM_Vector& stageSol,Double eta, SBM_Vector& Linform_nm1){
    
    SBM_Vector residual_vector(BaseMatrix::DOUBLE);
    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE); stageSol_temp = stageSol;
    double EnergyDerivative = 0.0;


    stageSol.Add(eta,solIncrement);
    algsys_->InitRHS(Linform_nm1);
    assemble_->AssembleNonLinRHS();
    algsys_->GetRHSVal( residual_vector );
    residual_vector.Inner(solIncrement,EnergyDerivative);
    stageSol = stageSol_temp;
    return EnergyDerivative;
  }

  double SolveStepEB::LineSearchArmijo(SBM_Vector& solIncrement, SBM_Vector& stageSol)   {

    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE); stageSol_temp = stageSol;
    SBM_Vector refRHS(BaseMatrix::DOUBLE);
    SBM_Vector stepRHS(BaseMatrix::DOUBLE);
    Double residual_ref = 0.0;
    Double residual_step = 0.0;

    const UInt nrEtas = 4;
    //const Double eta[nrEtas] = {1,0.6667,0.4444,0.2963,0.1975,0.1317,0.0878,0.0585,0.039,0.026,0.01733333333};
    const Double eta[nrEtas] = {1,0.1,0.01,0.001};
    const Double eta0 = 0;
    Double etaOpt = 0.0;

    Double delta = 1e-4;
    
    // obtain reference residual (start)
    stageSol.Add(eta0,solIncrement);
    algsys_->InitRHS();
    assemble_->AssembleLinRHS();
    assemble_->AssembleNonLinRHS();
    algsys_->GetRHSVal( refRHS );
    residual_ref = refRHS.NormL2();
    // obtain reference residual (end)
    
    for( UInt i=0; i<nrEtas; i++) {

      // obtain trial residual (start)
      stageSol.Add(eta[i],solIncrement);
      algsys_->InitRHS();
      assemble_->AssembleLinRHS();
      assemble_->AssembleNonLinRHS();
      algsys_->GetRHSVal( stepRHS );
      residual_step = stepRHS.NormL2();
      // obtain trial residual (end)

      // check if the trial residual has a sufficient decrease
      if (residual_step < (1.0 - delta*eta[i])*residual_ref) {
        etaOpt = eta[i];
        break;
      }
      if (i == nrEtas-1){ // after all possible line search parameter are tried take the last one and break the for loop
        etaOpt = eta[i];
        break;
      }
    }
    stageSol = stageSol_temp;
    return etaOpt;
  }

  double SolveStepEB::InexactLineSearch(SBM_Vector& solIncrement, SBM_Vector& stageSol)   {

    Double nr_gammas = 5;
    double gamma = 1.0;
    std::vector<Double> gamma_trial = {1,0.5,0.25,0.125,0.1}; 
    std::vector<Double> Energy;

    for( UInt idx=0; idx<nr_gammas; idx++) {
      Energy.push_back(GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, gamma_trial[idx]));
    }

    int closest_index = 0;
    Double closest_to_zero = Energy[0]; // Initialize with the first element
    Double min_distance = std::abs(Energy[0]);

    for (int idx = 1; idx < Energy.size(); ++idx) {
      Double distance = std::abs(Energy[idx]);
      if (distance < min_distance) {
        min_distance = distance;
        closest_index = idx;
      }
    }
    gamma = gamma_trial[closest_index];

    return gamma;
  }

  double SolveStepEB::InexactLineSearch(SBM_Vector& solIncrement, SBM_Vector& stageSol, SBM_Vector& Linform_nm1)  {
    
    double gamma = 1.0;
    std::vector<Double> gamma_trial = {1, 0.5, 0.25, 0.125, 0.1};
    std::vector<Double> Energy;

    for (UInt idx = 0; idx < gamma_trial.size(); idx++) {
        Energy.push_back(GetLineSearchDerivativeFunctionValue(
            solIncrement, stageSol, gamma_trial[idx], Linform_nm1));
    }

    int closest_index = 0;
    Double min_distance = std::abs(Energy[0]);
    for (int idx = 1; idx < (int)Energy.size(); ++idx) {
        Double distance = std::abs(Energy[idx]);
        if (distance < min_distance) {
            min_distance = distance;
            closest_index = idx;
        }
    }
    gamma = gamma_trial[closest_index];

    return gamma;
}

double SolveStepEB::BrentMethod(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double a, Double b)
{

    //Brent's method finds the root of the 1D function E'(eta) = r(theta_k + eta*eta_phi)^T * delta_eta, see also
    //DOI: 10.1109/TMAG.2024.3387354 and https://en.wikipedia.org/wiki/Brent's_method for Brent's method
    //or https://apps.dtic.mil/sti/tr/pdf/AD0726170.pdf

    const double tolInterval  = lineSearchTolerance_;   //lineSearchTolerance_ = 1e-3 default, for |b-a|
    const double tolFunction  = 1e-10;                  //for |f(b)| and |f(s)|    
    const int    maxIter = static_cast<int>(lineSearchMaxIter_);
    const double eps     = std::numeric_limits<double>::epsilon(); //previously min (10e-308) was used instead of epsilon (2.2x10e-16)

    //evaluation of the function at the interval endpoints
    double fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a);
    double fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b);

    //std::cout << "Fa=" << fa<< " Fb=" << fb<< " -> " << ((fa * fb > 0.0) ? "no root in [a,b]" : "ok")<< std::endl;

    //Initial swap: b should always be the best current approximation of root (i.e., |f(b)| <= |f(a)|)
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c  = a; //c... "oldest" of the three points
    double fc = fa;
    double d  = c; // d is the previous value of c, used to check if the step is a bisection step or an interpolation step
    bool mflag = true; //mflag: true if last step was bisection, false if interpolation (=IQI)
    double s  = 0.0; //new candidate point
    double fs = 0.0;

    for (int iter = 0; iter < maxIter; ++iter) {

        //convergence check on b: function value small enough
        if (std::abs(fb) <= tolFunction) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (function tolerance on b)" << std::endl;
            return b;
        }

        //convergence check on interval: interval small enough
        if (std::abs(b - a) <= 2.0 * eps * std::abs(b) + tolInterval) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (interval tolerance)" << std::endl;
            return b;
        }

        //interpolation step (IQI): only if function values at a,b,c are distinct enough to avoid
        //division small numbers or zero
        if (std::abs(fa - fc) > eps && std::abs(fb - fc) > eps) {
            s =   a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa); //secant method
        }

        double lower = std::min((3.0 * a + b) / 4.0, b);
        double upper = std::max((3.0 * a + b) / 4.0, b);

        bool cond1 = (s < lower || s > upper);
        bool cond2 = ( mflag && std::abs(s - b) >= std::abs(b - c) / 2.0);
        bool cond3 = (!mflag && std::abs(s - b) >= std::abs(c - d) / 2.0);
        bool cond4 = ( mflag && std::abs(b - c) < tolInterval);
        bool cond5 = (!mflag && std::abs(c - d) < tolInterval);

        //bisection if any condition holds, guarantees sufficient intervall reduction
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = 0.5 * (a + b); // bisection method
            mflag = true;
        } else {
            mflag = false;
        }

        fs = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, s); //function at new candidate point s

        if (std::abs(fs) <= tolFunction) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (function tolerance on s)" << std::endl;
            return s;
        }

        d  = c; //update d before overwriting c
        c  = b;
        fc = fb;

        //Update bracket: replace the endpoint that has the same sign as s,
        //preserving the sign change in [a,b].
        if (fa * fs < 0.0) {
            b  = s;
            fb = fs;
        } else {
            a  = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    //std::cout << "Brent reached maximum iterations (" << maxIter<< ") without convergence. Returning best of b/s."<< std::endl;
    //per wikipedia "output b or s": return the point with smaller |f|.
    return (std::abs(fb) <= std::abs(fs)) ? b : s;

}

double SolveStepEB::BrentMethod(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double a, Double b, SBM_Vector& Linform_nm1)  {

    //Brent's method finds the root of the 1D function E'(eta) = r(theta_k + eta*eta_phi)^T * delta_eta, see also
    //DOI: 10.1109/TMAG.2024.3387354 and https://en.wikipedia.org/wiki/Brent's_method for Brent's method
    //or https://apps.dtic.mil/sti/tr/pdf/AD0726170.pdf

    const double tolInterval  = lineSearchTolerance_;   //lineSearchTolerance_ = 1e-3 default, for |b-a|
    const double tolFunction  = 1e-10;                  //for |f(b)| and |f(s)|    
    const int    maxIter = static_cast<int>(lineSearchMaxIter_);
    const double eps     = std::numeric_limits<double>::epsilon(); //previously min (10e-308) was used instead of epsilon (2.2x10e-16)

    //evaluation of the function at the interval endpoints
    double fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a, Linform_nm1);
    double fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b, Linform_nm1);

    //std::cout << "Fa=" << fa<< " Fb=" << fb<< " -> " << ((fa * fb > 0.0) ? "no root in [a,b]" : "ok")<< std::endl;

    //Initial swap: b should always be the best current approximation of root (i.e., |f(b)| <= |f(a)|)
    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c  = a; //c... "oldest" of the three points
    double fc = fa;
    double d  = c; // d is the previous value of c, used to check if the step is a bisection step or an interpolation step
    bool mflag = true; //mflag: true if last step was bisection, false if interpolation (=IQI)
    double s  = 0.0; //new candidate point
    double fs = 0.0;

    for (int iter = 0; iter < maxIter; ++iter) {

        //convergence check on b: function value small enough
        if (std::abs(fb) <= tolFunction) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (function tolerance on b)" << std::endl;
            return b;
        }

        //convergence check on interval: interval small enough
        if (std::abs(b - a) <= 2.0 * eps * std::abs(b) + tolInterval) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (interval tolerance)" << std::endl;
            return b;
        }

        //interpolation step (IQI): only if function values at a,b,c are distinct enough to avoid
        //division small numbers or zero
        if (std::abs(fa - fc) > eps && std::abs(fb - fc) > eps) {
            s =   a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb));
        } else {
            s = b - fb * (b - a) / (fb - fa); //secant method
        }

        double lower = std::min((3.0 * a + b) / 4.0, b);
        double upper = std::max((3.0 * a + b) / 4.0, b);

        bool cond1 = (s < lower || s > upper);
        bool cond2 = ( mflag && std::abs(s - b) >= std::abs(b - c) / 2.0);
        bool cond3 = (!mflag && std::abs(s - b) >= std::abs(c - d) / 2.0);
        bool cond4 = ( mflag && std::abs(b - c) < tolInterval);
        bool cond5 = (!mflag && std::abs(c - d) < tolInterval);

        //bisection if any condition holds, guarantees sufficient intervall reduction
        if (cond1 || cond2 || cond3 || cond4 || cond5) {
            s = 0.5 * (a + b); // bisection method
            mflag = true;
        } else {
            mflag = false;
        }

        fs = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, s, Linform_nm1); //function at new candidate point s

        if (std::abs(fs) <= tolFunction) {
            //std::cout << "Brent converged after " << (iter + 1)<< " iterations (function tolerance on s)" << std::endl;
            return s;
        }

        d  = c; //update d before overwriting c
        c  = b;
        fc = fb;

        //Update bracket: replace the endpoint that has the same sign as s,
        //preserving the sign change in [a,b].
        if (fa * fs < 0.0) {
            b  = s;
            fb = fs;
        } else {
            a  = s;
            fa = fs;
        }

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    //std::cout << "Brent reached maximum iterations (" << maxIter<< ") without convergence. Returning best of b/s."<< std::endl;
    //per wikipedia "output b or s": return the point with smaller |f|.
    return (std::abs(fb) <= std::abs(fs)) ? b : s;

}

  void SolveStepEB::LineSearchHeavy(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch) {

    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 10;
    const Double eta[nrEtas] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

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

  void SolveStepEB::LineSearchHeavy(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch, SBM_Vector& Linform_nm1) {

    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = actSol;
    const UInt nrEtas = 10;
    const Double eta[nrEtas] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      actSol.Add( 1.0, solOld, eta[i], solIncrement);

      //store new solution
      solVec_ = actSol;

      // ======= residual computation ============
      algsys_->InitRHS(Linform_nm1);
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

}
