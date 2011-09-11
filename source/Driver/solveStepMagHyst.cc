// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMagHyst.hh"
#include "Utils/preisach.hh"
#include "assemble.hh"
#include "Forms/linearForm.hh"
#include "Forms/gradfieldop.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

#include "OLAS/algsys/algebraicSys.hh"

namespace CoupledField {

  SolveStepMagHyst::SolveStepMagHyst(StdPDE& apde) : StdSolveStep(apde) {
  }
  
  
  SolveStepMagHyst::~SolveStepMagHyst() {
  }
  
  void SolveStepMagHyst::StepTransNonLinHysteresis(PtrParamNode analysis_base) {

    UInt timeStep = (UInt) mParser_->Eval(mHandle_);
    std::cout << "In StepTransNonLinHysteresis, step:" << timeStep << std::endl;

    bool performOneMoreStep;
    Vector<Double> newSol( numEqns_ );
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> prevSol( numEqns_ );
    Vector<Double> prevSolDeriv, helpVec;
    Vector<Double> helpSol( numEqns_ );
  
    if ( prevLoadRHS_.GetSize() == 0 ) {
      prevLoadRHS_.Resize( numEqns_ );
      prevLoadRHS_.Init();
    }

    // get solution of last time step
    algsys_->GetSolutionVal( prevSol );      

    //get previous 1st order derivative of solution
    prevSolDeriv = TS_alg_->GetDeriv(FIRST_DERIV);

    // perform predictor step
    if ( TS_alg_== NULL ) {
      EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMagHyst::StepTransNonLin" );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(prevSol);
    }
    
    //! account for Dirichlet BCs
    PDE_.SetBCs();

    // handle loads
    // 1)  assemble to RHS loads for time step (n+1) and save it to RhsLinVal_ 
    assemble_->AssembleLinRHS(); 
    algsys_->GetRHSVal( RhsLinVal_ );

    //2) substract the loads of time step n
    helpVec     = RhsLinVal_;
    RhsLinVal_ -= prevLoadRHS_;
    algsys_->InitRHS( RhsLinVal_ );

    //3) store loads of time step (n+1) for the next time step
    prevLoadRHS_ = helpVec;

    //4) add mass matrix multiplied by derivativ of previous solution
    algsys_->UpdateRHS(MASS, prevSolDeriv );

    //5) update RHS according to time stepping
    TS_alg_->UpdateRHS();

    //6) this will be our constant RHS for this time step
    algsys_->GetRHSVal( RhsLinVal_ );

    //L2 of linear RHS
    Double RhsLinL2Norm = RhsLinVal_.NormL2();
    //std::cout << "RHSconst:\n" << RhsLinVal_ << std::endl;


    //initialize old solution
    oldSol = prevSol;

    // set iteration counter
    UInt iterationCounter=0;

    do {
      iterationCounter++;

      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
        std::cout << std::endl << "Time step:   "  << actStep_
                  << "  ,Iterations: " << iterationCounter << std::endl;
      else 
        std::cout << "Iter:  " << iterationCounter << std::endl;

      PtrParamNode analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);
      
      // set solution of previous iteration
      if (iterationCounter > 1 )
        oldSol = newSol;


      // RHS is already set up!!
      if ( iterationCounter > 1 ) {
        // setup linear part of right hand side 
        algsys_->InitRHS(RhsLinVal_);
      }

      if ( timeStep == 1 || iterationCounter > 1 ) {
        //assemble matrices
        assemble_->AssembleMatrices();
      }

      // add stiffness multiplied by solution of previous time step;
      algsys_->UpdateRHS(STIFFNESS, prevSol );

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->BuildInDirichlet();
  
      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
      algsys_->Solve(analysis_id);

      algsys_->GetSolutionVal( newSol );

      std::cout << "New Solution: iter  " << iterationCounter << "\n" << newSol << std::endl;
      //    std::cout << "Old Solution:\n" << oldSol << std::endl;

      //put new solution to sol_
      sol_->SetAlgSysVector(newSol);  
    

      // compute L2-Norm of error between last incremental solution and
      // actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newSol.GetSize(); i++)
        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
    
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newSol.NormL2();
    
      Double incrementalErr;
      if (actSolL2Norm > 1)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;


      // compute residual norm

      // 1) Set linear part to RHS
      algsys_->InitRHS(RhsLinVal_ );

      //2) assemble matrices
      assemble_->AssembleMatrices();

      //3) add stiffness multiplied by solution of previous time step;
      algsys_->UpdateRHS(STIFFNESS, prevSol );

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(newSol);

      //4) stiffness matrix
      helpSol = -newSol;
      algsys_->UpdateRHS(STIFFNESS, helpSol );


      Vector<Double> actRHS;
      Double residualL2Norm;
      algsys_->GetRHSVal( actRHS );
      
      // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
      residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
     
       Double residualErr=0;
      if ( RhsLinL2Norm > 1.0 )
        residualErr    = residualL2Norm /  RhsLinL2Norm;
      else
        residualErr    = residualL2Norm;

 
      nonLinLogging_ = true;
      if ( nonLinLogging_ == true )
        WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr);

      //    std::cout << "ResNorm=" << residualNorm << "  incrNorm=" 
      //        << incrementalErr << std::endl;
    
      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_ ); //|| residualErr > residualStopCrit_);      
      
    } while( performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    if ( iterationCounter >= nonLinMaxIter_ ) {
      EXCEPTION( "Number of nonlinear iterations exceeds limit "
               << "nonLinearMaxIter_ = "
               << nonLinMaxIter_ );
    }


    std::cout << "\n Do compute Prev Val\n" << std::endl;
    //set the current values to the previous for the next time step!
    SetPreviousVals4Hyst();

    //perform corrector step  
    TS_alg_->Corrector(newSol);
  }


  void SolveStepMagHyst::SetPreviousVals4Hyst() {
  
    EXCEPTION("Magnetics with hysteresis currently not supported");

//     Vector<Double> LCoord, field;
//     UInt pdeElem;

//     for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {

//       RegionIdType actRegion = subdoms_[actSD];
//       Hysteresis* hyst = materialData_[actRegion]->getHysteresis();
 
//       if ( hyst!= NULL ) {

//         CurlCurlNode2DInt *curlOp = new CurlCurlNode2DInt( materialData_[actRegion], isaxi_);

//         ElemList actSDList(ptgrid_ );
//         actSDList.SetRegion( actRegion );

//         EntityIterator it = actSDList.GetIterator();
//         UInt iel = 0;
//         for ( it.Begin(); !it.IsEnd(); it++, iel++) {
//           //compute the magnetic field intensity
//           it.GetElem()->ptElem->GetCoordMidPoint(LCoord);

//           Matrix<Double> CornerCoords;
//           ptgrid_->GetElemNodesCoord( CornerCoords, it.GetElem()->connect, 
//                                       true );
//           // get element solution
//           Vector<Double> elSol;
//           sol_->GetElemSolution(elSol, it);

//           Vector<double> field;
//           if( LCoord.GetSize() ==3 ) {
//             field.Resize(3);
//             EXCEPTION("3D magnetic with hysteresis not implemented!");
//           }
//           else {
//             Matrix<Double> bMat;
//             //set element info
//             curlOp->ExtractElemInfo( it );
//             curlOp->SetIntPoint( LCoord );
//             curlOp->calcBMat(bMat, 0, CornerCoords);
//             curlOp->UnsetIntPoint();     
//             field = bMat * elSol;

//             // Account for curl 
//             Double temp = field[0];
//             if ( isaxi_ ) {
//               field[0] = -field[1];
//               field[1] = temp;
//             } else {
//               field[0] = field[1];
//               field[1] = -temp;
//             }
//           }

//           pdeElem = it.GetElem()->elemNum;

//           // std::cout << "New Bfield: \n " << field << std::endl << std::endl;
//           //set the values
//           materialData_[actRegion]->SetPreviousInvVecHystVal( pdeElem, field );
//         }
//         delete curlOp;
//       }
//     }
  }


  void SolveStepMagHyst::StepTransNonLinHysteresisDiff(PtrParamNode analysis_base) {

    std::cout << "In :StepTransNonLinHysteresisDiff" << std::endl;

    bool performOneMoreStep;

    Vector<Double> actSol( numEqns_ );
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> helpSol( numEqns_ );
  

    //get actual solution  
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    actSol = solhelp->GetAlgSysVector();

    // perform predictor step
    if ( TS_alg_== NULL ) {
      EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMagHyst::StepTransNonLin" );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs();


    // handle loads

    // 1)  assemble to RHS loads for time step (n+1) and save it to RhsLinVal_ 
    assemble_->AssembleLinRHS(); 
  
    //2) update RHS according to time stepping
    TS_alg_->UpdateRHS();

    //3) this will be our constant RHS for this time step
    algsys_->GetRHSVal( RhsLinVal_ );

    //L2 of linear RHS
    Double RhsLinL2Norm = RhsLinVal_.NormL2();
    //std::cout << "RHSconst:\n" << RhsLinVal_ << std::endl;

    // set iteration counter
    UInt iterationCounter=0;

    //initialize old solution
    oldSol = actSol;

    do {
      iterationCounter++;

      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
        std::cout << std::endl << "Time step:   "  << actStep_
                  << "  ,Iterations: " << iterationCounter << std::endl;
      else 
        std::cout << "Iter:  " << iterationCounter << std::endl;

      PtrParamNode analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);
            
      // set solution of previous iteration
      if (iterationCounter > 1 )
        oldSol = actSol;

      // RHS is already set up!!
      if ( iterationCounter > 1 ) {
        // setup linear part of right hand side 
        algsys_->InitRHS(RhsLinVal_);
      }

      //asemble matrices
      assemble_->AssembleMatrices();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->BuildInDirichlet();
  
      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
      algsys_->Solve(analysis_id);

      algsys_->GetSolutionVal( actSol );
      //      std::cout << "New Solution:\n" << newSol << std::endl;
      //    std::cout << "Old Solution:\n" << oldSol << std::endl;

      //put new solution to sol_
      sol_->SetAlgSysVector(actSol);  
    

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


      // compute residual norm

      // 1) Set linear part to RHS
      algsys_->InitRHS(RhsLinVal_ );

      //2) assemble matrices
      assemble_->AssembleMatrices();

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(actSol);

      //4) stiffness matrix
      helpSol = -actSol;
      algsys_->UpdateRHS(STIFFNESS, helpSol );


      Vector<Double> actRHS;
      Double residualL2Norm;
      algsys_->GetRHSVal( actRHS );
      
      // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
      residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
     
      Double residualErr;
      if ( RhsLinL2Norm > 1.0 )
        residualErr    = residualL2Norm /  RhsLinL2Norm;
      else
        residualErr    = residualL2Norm;

 
      nonLinLogging_ = true;
      if ( nonLinLogging_ == true )
        WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr);

      //    std::cout << "ResNorm=" << residualNorm << "  incrNorm=" 
      //        << incrementalErr << std::endl;
    
      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_ || residualErr > residualStopCrit_);      
      
    } while( performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    if ( iterationCounter >= nonLinMaxIter_ ) {
      EXCEPTION( "Number of nonlinear iterations exceeds limit "
               << "nonLinearMaxIter_ = "
               << nonLinMaxIter_ );
    }

    //set the current values to the previous for the next time step!
    SetPreviousVals4Hyst();

    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }

} // end of namespace

