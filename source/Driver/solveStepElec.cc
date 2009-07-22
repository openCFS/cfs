// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepElec.hh"
#include "Utils/preisach.hh"
#include "assemble.hh"
#include "Forms/linearForm.hh"
#include "Forms/gradfieldop.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  SolveStepElec::SolveStepElec(StdPDE& apde) : StdSolveStep(apde) {
    doInit_ = true;
  }
  
  
  SolveStepElec::~SolveStepElec() {
  }
  

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepElec::PreStepStatic() {


    // Set right-hand side to zero
    //
    // Note: Though this is PreStepStatic, this is most important
    //       for transient analysis ;-)
    algsys_->InitRHS();

  }



  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void SolveStepElec::SolveStepStatic(InfoNode* analysis_id) {
    if ( isHyst_ ) 
      StepStaticNonLinEpsDiff(analysis_id);
    else 
      StepStaticLin(analysis_id);
  }


  void SolveStepElec::StepStaticNonLinEpsDiff(InfoNode* analysis_base) {
    REFACTOR;

//    bool performOneMoreStep;
//    UInt iterationCounter=0;
//  
//    Vector<Double> newSol( numEqns_ ); 
//    Vector<Double> oldSol( numEqns_ );
//    Vector<Double> solPrev( numEqns_ );
//    Vector<Double> incrSol( numEqns_ );
//    Vector<Double> coeff;
//
//    oldSol.Init(0);
//    newSol.Init(0);
//  
//    //save solution of previous time step  
//    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
//    solhelp->GetAlgSysVector(solPrev);
//
//    //clear RHS
//    algsys_->InitRHS();
//
//    //set BCs
//    PDE_.SetBCs();
//
//    // stores this as linear part of RHS
//    algsys_->GetRHSVal( RhsLinVal_ );
//
//    do {
//      iterationCounter++;
//      // for every time step write out number of iteration loops to standard out
//      if (iterationCounter == 1)
//        std::cout << std::endl << "Time step:   "  << actStep_
//                  << "  ,Iterations: " << iterationCounter << std::endl;
//      else 
//        std::cout << "Iter:  " << iterationCounter << std::endl;
//
//      InfoNode* analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);
//      
//      // set solution of previous iteration
//      if (iterationCounter == 1 ) {
//        oldSol = solPrev;
//      }
//      else {
//        oldSol = newSol;
//      }
//    
//      // Set linear part of RHS
//      //    algsys_->InitRHS(RhsLinVal_.GetPointer());
//
//      //clear RHS
//      algsys_->InitRHS();
//
//      //perform new assembly
//      assemble_->AssembleMatrices();
//      algsys_->ConstructEffectiveMatrix(matrix_factor_);
//
//      // set changing part of RHS
//      algsys_->UpdateRHS(SYSTEM,solPrev);
//
//      // build in the Dirichlet vales in system mmatrix and rhs
//      algsys_->BuildInDirichlet();
//
//      //get RHS
//      Vector<Double> RHS;
//      algsys_->GetRHSVal( RHS );
//      Double residualNorm = PDE_.GetRhsL2Norm( RHS );
//
//      algsys_->SetupSolver(analysis_id);
//      algsys_->SetupPrecond();
//    
//      algsys_->Solve(analysis_id);
//      algsys_->GetSolutionVal( incrSol );
//
//      //Double alpha = 1;
//      //    newSol = oldSol + incrSol * alpha;
//      newSol = incrSol;
//      
//      //put new solution to sol_
//      sol_->SetAlgSysVector(newSol);  
//    
//      // compute L2-Norm of error between last incremental solution and
//      // actual incremental solution
//      Double solIncrL2Norm=0;
//      for (UInt i=0; i<newSol.GetSize(); i++)
//        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
//    
//      solIncrL2Norm = sqrt(solIncrL2Norm);
//      Double actSolL2Norm = newSol.NormL2();
//    
//      Double incrementalErr;
//      if (actSolL2Norm > 1)
//        incrementalErr = solIncrL2Norm / actSolL2Norm;
//      else
//        incrementalErr = solIncrL2Norm;
//    
//      // output of norms and data
//      nonLinLogging_ = true;
//      if ( nonLinLogging_ == true ) {
//        Info->WriteNonLinIter(pdename_, iterationCounter, residualNorm,
//                              incrementalErr);
//      }
//
//      //    std::cout << "ResNorm=" << residualNorm << "  incrNorm=" 
//      //        << incrementalErr << std::endl;
//    
//      // boolean variable, holds condition if another iteration step is necessary
//      performOneMoreStep = 
//	(incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
//      
//    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
//
//    if ( iterationCounter >= nonLinMaxIter_ ) {
//      EXCEPTION( "Number of nonlinear iterations exceeds limit "
//               << "nonLinearMaxIter_ = "
//               << nonLinMaxIter_ );
//    }
//
//    //set the current values to the previous for the next time step!
//    SetPreviousVals4Hyst();
  }


  void SolveStepElec::SetPreviousVals4Hyst() {
  
    REFACTOR;
    
//    //we assume, that the actual solution is stored in sol_!
//    NodeStoreSol<Double> * solhelp = 
//      dynamic_cast<NodeStoreSol<Double>*>(sol_);
//
//    GradientFieldOp<Double> * FieldOp = 
//      new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnMap_,
//                                  *solhelp, ELEC_POTENTIAL, 
//                                  results_[0], isaxi_);
//
//    Vector<Double> LCoord, Efield;
//    Double Ecomp;
//    UInt pdeElem;
//
//    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
//
//      RegionIdType actRegion = subdoms_[actSD];
//      Hysteresis* hyst = materialData_[actRegion]->getHysteresis();
//      if ( hyst!= NULL ) {
//        //get direction of polarization
//        std::string str;
//        materialData_[actRegion]->GetScalar(str, P_DIRECTION);
//        Directions dirP;
//        String2Enum(str,dirP);
//
//        ElemList actSDList(ptgrid_ );
//        actSDList.SetRegion( actRegion );
//
//        EntityIterator it = actSDList.GetIterator();
//        UInt iel = 0;
//        for ( it.Begin(); !it.IsEnd(); it++, iel++) {
//
//          //compute the electric field intensity
//          it.GetElem()->ptElem->GetCoordMidPoint(LCoord);
//          FieldOp->CalcElemGradField( Efield, it, LCoord, 1);
//
//          //get correct component of electric field for scalar Preisach model
//          Ecomp   = Efield[dirP]; 	  
//          pdeElem = it.GetElem()->elemNum;
//
//          //set the values
//          materialData_[actRegion]->SetPreviousHystVal( pdeElem, Ecomp );
//        }  
//      }
//    }
  }



} // end of namespace

