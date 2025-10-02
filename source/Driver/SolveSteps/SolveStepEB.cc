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
            //std::cout <<"residualErr0: " << residualErr0 << std::endl;
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
          performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
          if ( performOneMoreStep == 0){
            // now that we have reached our convergence threshold for this timestep, let's save the states in our model
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
      // PSEUDO TIME STEPPING (END)
      // ###############################################################
      // ###############################################################
      // REAL TIME STEPPING (START)
      // ###############################################################
      if (pseudo_time_stepping_ == 0) {
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
          performOneMoreStep = (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
          if ( performOneMoreStep == 0){
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

    Double bottom_interval = 0.0 + 1e-3;
    Double top_interval = 1; 
    Double gamma = 0;
    
    gamma = BrentMethod(solIncrement, stageSol, bottom_interval, top_interval);
    return gamma;
  }
  double SolveStepEB::ExactLineSearch(SBM_Vector& solIncrement, SBM_Vector& stageSol, SBM_Vector& Linform_nm1){

    Double bottom_interval = 0.0 + 1e-3;
    Double top_interval = 1; 
    Double gamma = 0;
    
    gamma = BrentMethod(solIncrement, stageSol, bottom_interval, top_interval, Linform_nm1);
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



  double SolveStepEB::BrentMethod(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double a, Double b){

    double Fa, Fb, Fc;
    double c;
    double tolerance = 1e-3;
    double max_iter = 1e3;
    double iter_counter = 0;
    double prev_step = 0;
    double tol_act = 0;
    double eps = std::numeric_limits<Double>::min();
    double new_step = 0;
    double cb = 0;
    double t1 = 0;
    double t2 = 0;
    double p = 0;
    double q = 0;

    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE); stageSol_temp = stageSol;

    Fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a);
    Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b);
    c = a; 
    Fc = Fa;

    while(iter_counter < max_iter){
        iter_counter  = iter_counter + 1;

        prev_step = b - a;

        if (std::abs(Fc) < std::abs(Fb)){ // swap for b to be best approximation
            a = b ;
            b = c ;
            c = a ;
            Fa = Fb ;
            Fb = Fc ;
            Fc = Fa ;
        }

        tol_act =  2*eps*std::abs(b) + tolerance/2 ;
        new_step = (c-b)/2 ;

        if ((std::abs(new_step) <= tol_act) || (std::abs(Fb) < eps)){
            return b;
        }
        if ((std::abs(prev_step) >= tol_act) && ((std::abs(Fa-Fb))<eps)){
            cb = c - b;
            if (std::abs(a-c) < eps){ // linear interpolation, only two points available
                t1 = Fb/Fa;
                p = cb*t1;
                q = 1.0 - t1;
            }else { // three points, do quadratic inverse  interpolation
                a = Fa/Fc;
                t1 = Fb/Fc;
                t2 = Fb/Fa;
                p = t2*( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                q = (q-1.0)*(t1-1.0)*(t2-1.0);
            }

            if (p > 0){
                q = -q;
            }else {
                p = -p;
            }

            if ((p < ( 0.75*cb*q-std::abs(tol_act*q)/2.0 )) && (p < abs(prev_step*q/2.0))){
                new_step = p/q;
            }
        }

        // step must be at least as large as tolerance
        if (std::abs(new_step) < tol_act){
            if (new_step > 0){
                new_step = tol_act ;
            } else{
                new_step = -tol_act ;
            }
        }

        a = b;
        Fa = Fb;
        b = b + new_step;
        Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b);
        if (( Fb > 0.0 && Fc > 0.0 ) || ( Fb < 0.0 && Fc < 0.0 )){
            c = a;
            Fc = Fa;        
        }
    }
    return b;
  }
  double SolveStepEB::BrentMethod(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double a, Double b, SBM_Vector& Linform_nm1){

    double Fa, Fb, Fc;
    double c;
    double tolerance = 1e-3;
    double max_iter = 1e3;
    double iter_counter = 0;
    double prev_step = 0;
    double tol_act = 0;
    double eps = std::numeric_limits<Double>::min();
    double new_step = 0;
    double cb = 0;
    double t1 = 0;
    double t2 = 0;
    double p = 0;
    double q = 0;

    SBM_Vector stageSol_temp(BaseMatrix::DOUBLE); stageSol_temp = stageSol;

    Fa = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, a, Linform_nm1);
    Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b, Linform_nm1);
    c = a; 
    Fc = Fa;

    while(iter_counter < max_iter){
        iter_counter  = iter_counter + 1;

        prev_step = b - a;

        if (std::abs(Fc) < std::abs(Fb)){ // swap for b to be best approximation
            a = b ;
            b = c ;
            c = a ;
            Fa = Fb ;
            Fb = Fc ;
            Fc = Fa ;
        }

        tol_act =  2*eps*std::abs(b) + tolerance/2 ;
        new_step = (c-b)/2 ;

        if ((std::abs(new_step) <= tol_act) || (std::abs(Fb) < eps)){
            return b;
        }
        if ((std::abs(prev_step) >= tol_act) && ((std::abs(Fa-Fb))<eps)){
            cb = c - b;
            if (std::abs(a-c) < eps){ // linear interpolation, only two points available
                t1 = Fb/Fa;
                p = cb*t1;
                q = 1.0 - t1;
            }else { // three points, do quadratic inverse  interpolation
                a = Fa/Fc;
                t1 = Fb/Fc;
                t2 = Fb/Fa;
                p = t2*( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                q = (q-1.0)*(t1-1.0)*(t2-1.0);
            }

            if (p > 0){
                q = -q;
            }else {
                p = -p;
            }

            if ((p < ( 0.75*cb*q-std::abs(tol_act*q)/2.0 )) && (p < abs(prev_step*q/2.0))){
                new_step = p/q;
            }
        }

        // step must be at least as large as tolerance
        if (std::abs(new_step) < tol_act){
            if (new_step > 0){
                new_step = tol_act ;
            } else{
                new_step = -tol_act ;
            }
        }

        a = b;
        Fa = Fb;
        b = b + new_step;
        Fb = GetLineSearchDerivativeFunctionValue(solIncrement, stageSol, b, Linform_nm1);
        if (( Fb > 0.0 && Fc > 0.0 ) || ( Fb < 0.0 && Fc < 0.0 )){
            c = a;
            Fc = Fa;        
        }
    }
    return b;
  }
}
