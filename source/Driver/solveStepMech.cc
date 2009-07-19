#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMech.hh"
#include "assemble.hh"
#include "PDE/StdPDE.hh"
#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  SolveStepMech::SolveStepMech(StdPDE& apde) : StdSolveStep(apde) {
    //predRelax_= PDE_.GetPredictorRelax();
    //calcPredOutputCoupling_ = false;

    //for aitken relaxation in nonlinear mechanics
    //aitkenMu_=0.1;
    oldInc_.Resize(numEqns_);
    oldInc_.Init();
  }
  
  
  SolveStepMech::~SolveStepMech() {
  }

//  void SolveStepMech::PredictorStep(){
//
//    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
//
//    if(iterCoupledCounter==0){
//      Vector<Double> solHelp; 
//      sol_->GetAlgSysVector(solHelp);
//
//      Vector<Double> solderiv1Help=TS_alg_->GetDeriv1(); 
//
//      Vector<Double> solpredRelax, solderiv1predRelaxed;
//
//      if(predRelax_>1e-3){
//        TS_alg_->RelaxedPredictor(solHelp,solpredRelax,solderiv1predRelaxed,predRelax_);
//
//        sol_->SetAlgSysVector(solpredRelax);
//        TS_alg_->SetDeriv1(solderiv1predRelaxed);
//
//        if(calcPredOutputCoupling_)
//          PDE_.CalcPredictorOutputCoupling();
//      
//        sol_->SetAlgSysVector(solHelp);
//        TS_alg_->SetDeriv1(solderiv1Help);
//      }
//      else{
//        TS_alg_->Predictor(solHelp);
//        if(calcPredOutputCoupling_)
//          PDE_.CalcPredictorOutputCoupling();
//      }
//    }
//  }


  void SolveStepMech::StepTransLin(InfoNode* analysis_base) {

    //account for RHS
    assemble_->AssembleLinRHS();

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
    bool isIterCoupled    = PDE_.IsIterCoupled();

    InfoNode* analysis_id = !isIterCoupled ? analysis_base :
          BaseDriver::CreateAnalysisIdChild(analysis_base, "iterCoupled", iterCoupledCounter);    

    // perform predictor step: if we have an iterative coupled 
    // PDE-system, we should perform the predictor state just 
    // in the first iteration
    if ( isIterCoupled == false || iterCoupledCounter == 0 ) {        
      Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      TS_alg_->Predictor(solHelp);
    }

    assemble_->AssembleMatrices();

    if (assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }

    // update Right Hand Side
    TS_alg_->UpdateRHS();

    // set boundary conditions
    PDE_.SetBCs();
    algsys_->BuildInDirichlet();

    if (assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
    }
    algsys_->Solve(analysis_id);

    Vector<Double> solHelp;
    algsys_->GetSolutionVal(solHelp);
    PDE_.SaveSolution(solHelp.GetPointer(), solHelp.GetSize());
   
    if( isIterCoupled == false ) {
      Vector<Double> & solHelp = 
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      TS_alg_->Corrector(solHelp);
    }
    if ( isIterCoupled ) {
      iterCoupledCounter++;
    }
  }

  // Newton-Raphson Iteration
//  void SolveStepMech::StepTransNonLin() {
//
//    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
//    bool isIterCoupled    = PDE_.IsIterCoupled();
//
//    Double *solPtr;
//
//    bool performOneMoreStep;
//
//    Vector<Double> solInc( numEqns_ );
//    Vector<Double> actSol( numEqns_ );
//    Vector<Double> oldSol( numEqns_ );
//    actSol.Init();
//
//    UInt numNodes;
//    numNodes = sol_->GetNumNodes();
//  
//    // Cast BaseStoreSol into StoreSol<Double>,
//    // since this function is only called
//    // in the transient case
//
//    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
//
//    //get actual solution  
//    actSol = solhelp->GetAlgSysVector();
//    oldSol=actSol;
//
//    // currently just for testing!!
//    // loop over load factor
//    Double loadFactor=0;
//
//    UInt nrLoadSteps = SetDeltaLinRHS();
//    Double loadInc=1.0/Double(nrLoadSteps);
//    Info->PrintF(pdename_, "\n");
//    Info->PrintF(pdename_, "Nr. of Load Steps: %d \t increment : %g\n",
//                 nrLoadSteps, loadInc);
//
//
//    if( TS_alg_ != NULL && nrLoadSteps > 1) {
//      TS_alg_->setSubSteps(nrLoadSteps);
//      matrix_factor_ = TS_alg_->GetEffSysMatFactors();
//    }
//
//    for ( UInt iload=0; iload<nrLoadSteps; iload++ ) {
//      loadFactor += loadInc;
//      //loadFactor = 1.0;
//      Info->PrintF(pdename_, "\n");
//      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);
//      
//      // perform predictor step
//      if ( TS_alg_== NULL ) {
//        Error( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin",
//               __FILE__, __LINE__ );
//      }
//      else {
//          if(nrLoadSteps>1)
//            TS_alg_->Predictor(actSol);
//      }
//
//      //! account for Dirichlet BCs
//      PDE_.SetBCs();
//
//      // setup right hand side 
//      Double RhsLinL2Norm = SetLinRHS(loadFactor); 
//
//      algsys_->InitRHS(RhsLinVal_.GetPointer());
//
//      // inner forces due to nonlin formulation
//      assemble_->AssembleNonLinRHS( actTime_ );  
//      //Update RHS (mass matrix on right hand side)
//      assemble_->AssembleMatrices();
//      TS_alg_->UpdateRHS(actSol);
//
//      // set iteration counter
//      UInt iterationCounter=0;
//      Double iterTime=0;
//      firstEi_=true; 
//      do {
//        iterationCounter++;
//        iterTime+=0.0001;
//        
//        // RHS is already set up!!
//
//        if(iterationCounter !=1)
//          assemble_->AssembleMatrices();
//        algsys_->ConstructEffectiveMatrix(matrix_factor_);
//        algsys_->BuildInDirichlet();
//  
//        algsys_->SetupPrecond();
//        algsys_->SetupSolver();
//        algsys_->Solve();
//
//        // new solution is only an increment of the full solution =============
//        algsys_->GetSolutionVal( solPtr );
//        StoreAlgsysToVec(solInc, solPtr );
//
//        Double residualL2Norm;
//        Double etaLineSearch = 1.0;
//        Double Ei=0.0;
//        if ( lineSearch_ == "none" ) {
//          oldInc_ = solInc;
//          actSol += solInc;
//        }
//        else {
//          residualL2Norm = LineSearch(solInc,actSol,etaLineSearch,true);
//        }
//      
//        //store A_(n+1) in the solution-object sol_
//        sol_->SetAlgSysVector(actSol);
//
//	//PDE_.WriteResultsInFile(iterationCounter, actTime_+iterTime, 0, 0);
//
//        if ( lineSearch_ == "none" ) {
//          // calculation of error norms
//          // recalculate RHS with new values to get new residual (f^(k+1))=======
//          algsys_->InitRHS(RhsLinVal_.GetPointer());
//
//          // inner forces due to nonlin formulation
//          assemble_->AssembleNonLinRHS( actTime_ );
//
//          //Update RHS (mass matrix on right hand side)
//          TS_alg_->UpdateRHS(actSol);
//
//          Vector<Double> actRHS;
//          Double *rhsPtr;
//          algsys_->GetRHSVal( rhsPtr );
//          StoreAlgsysToVec( actRHS, rhsPtr );
//          
//          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
//          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
//
//          solInc.Inner(actRHS,Ei);
//          Ei=abs(Ei);
//          //Ei=std::sqrt(Ei);
//          //std::cerr << "Ei=" << Ei << std::endl;
//
//        }
//        if(firstEi_){
//        	firstEi_=false;
//        	Ei0_=Ei;
//        }
//        Ei/=Ei0_;
//
//
//        Double residualErr;
//        if ( RhsLinL2Norm > 1.0 )
//          residualErr    = residualL2Norm /  RhsLinL2Norm;
//        else{
//          residualErr    = residualL2Norm;
//          //std::cerr << "RhsLinL2Norm < 1" << std::endl;
//        }
//        // calculate incremental error
//        Double solIncrL2Norm = solInc.NormL2();
//        Double actSolL2Norm = actSol.NormL2();
//        Double incrementalErr;
//      
//        if ( actSolL2Norm > 1.0)
//          incrementalErr = solIncrL2Norm / actSolL2Norm;
//        else
//          incrementalErr = solIncrL2Norm;
//
//        // --------------------------------------------------------------------
//        // output of norms and data
//        // --------------------------------------------------------------------
//        if ( nonLinLogging_ == true ) {
//          Info->WriteNonLinIter2(pdename_, iterationCounter, residualErr,
//                                incrementalErr, Ei, etaLineSearch);
//        }
//
//        // boolean variable, holds condition if another iteration step
//        // is necessary
//        performOneMoreStep = 
//          (incrementalErr > incStopCrit_)||(Ei > residualStopCrit_);
//        //if( Ei < 1.0e-6)
//        //  performOneMoreStep = false;
//      
//      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
//      
//      //SetOldRHS(actRHS);
//
//      //} // load step loop
//      
//      //if( isIterCoupled == false ){
//      //perform corrector step  
//      if(nrLoadSteps>1)
//        TS_alg_->Corrector(actSol);
//      //}
//    }
//    if( TS_alg_ != NULL && nrLoadSteps>1 ) {
//      TS_alg_->resetDeltaT();
//      matrix_factor_ = TS_alg_->GetEffSysMatFactors();
//    }
//    
//    if ( isIterCoupled ) {
//      iterCoupledCounter++;
//    }
//    
//  }

  
  
  
  
  
  
  
  
  
  
  
  
  //Fixpoint Iteration
//   void SolveStepMech::StepTransNonLin() {

//     ENTER_FCN( "SolveStepMech::StepTransNonLin" );

//     assemble_->AssembleLinRHS( actTime_ );

//     UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
//     bool isIterCoupled    = PDE_.IsIterCoupled();

//     Double *solPtr;

//     bool performOneMoreStep;

//     Vector<Double> solInc( numEqns_ );
//     Vector<Double> actSol( numEqns_ );
//     Vector<Double> newSol( numEqns_ );
//     actSol.Init();
//     newSol.Init();

//     UInt numNodes;
//     numNodes = sol_->GetNumNodes();
  
//     // Cast BaseStoreSol into StoreSol<Double>,
//     // since this function is only called
//     // in the transient case

//     NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

//     //get actual solution  
//     actSol = solhelp->GetAlgSysVector();
//     //oldSol=actSol;

//     // perform predictor step
//     if ( TS_alg_== NULL ) {
//       Error( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin",
//              __FILE__, __LINE__ );
//     }
//     else {
//       if ( isIterCoupled == false) {        
//       //if ( isIterCoupled == false || iterCoupledCounter == 0 ) {
//         //compute predictors
//         TS_alg_->Predictor(actSol);
//       }
//     }

//     //! account for Dirichlet BCs
//     PDE_.SetBCs( actTime_ );

//     // inner forces due to nonlin formulation
//     //assemble_->AssembleNonLinRHS( actTime_ );  
//     //Update RHS (mass matrix on right hand side)
//     assemble_->AssembleMatrices();
//     TS_alg_->UpdateRHS();

//     // set iteration counter
//     UInt iterationCounter=0;
    
    
//     do {
//       iterationCounter++;

//       // RHS is already set up!!

//       //assemble_->AssembleMatrices();
//       if ( iterationCounter != 1){
//         assemble_->AssembleMatrices();
//       }
//       algsys_->ConstructEffectiveMatrix(matrix_factor_);
//       algsys_->BuildInDirichlet();
  
//       algsys_->SetupPrecond();
//       algsys_->SetupSolver();
//       algsys_->Solve();

//       algsys_->GetSolutionVal( solPtr );
//       StoreAlgsysToVec(newSol, solPtr );
//       sol_->SetAlgSysVector(newSol);

//       solInc=newSol-actSol;

//       actSol=newSol;

//       // calculate incremental error
//       Double solIncrL2Norm = solInc.NormL2();
//       Double actSolL2Norm = actSol.NormL2();
//       Double incrementalErr;
      
//       //         if ( actSolL2Norm > 1.0)
//       //           incrementalErr = solIncrL2Norm / actSolL2Norm;
//       //         else
//       incrementalErr = solIncrL2Norm;

//       // --------------------------------------------------------------------
//       // output of norms and data
//       // --------------------------------------------------------------------
//       if ( nonLinLogging_ == true ) {
//         Info->WriteNonLinIter(pdename_, iterationCounter, actSolL2Norm,
//                               incrementalErr, 1.0);
//       }

//       // boolean variable, holds condition if another iteration step
//       // is necessary
//       performOneMoreStep = 
//         (incrementalErr > incStopCrit_);
      
//     } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  


//     if( isIterCoupled == false ){
//       //perform corrector step  
//       TS_alg_->Corrector(actSol);
//     }
//     if ( isIterCoupled ) {
//       iterCoupledCounter++;
//     }
    
//   }

} // end of namespace

