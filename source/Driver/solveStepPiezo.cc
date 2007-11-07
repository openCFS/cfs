// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepPiezo.hh"
#include "Utils/preisach.hh"
#include "Utils/jiles.hh"
#include "assemble.hh"
#include "Forms/gradfieldop.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "Utils/result.hh"

namespace CoupledField {

  SolveStepPiezo::SolveStepPiezo(StdPDE& apde) : StdSolveStep(apde)
  {

    doInit_ = true;
  }


  SolveStepPiezo::~SolveStepPiezo() {
  }
 

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepPiezo:: PreStepTrans()
  {


    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();


    // save solution of time step n as previous solution
    //save solution of previous time step 
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

    //store solution for (n)
    PDE_.SavePrevSolution( solHelp.GetPointer(), solHelp.GetSize() );

  }



  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void SolveStepPiezo::SolveStepTrans() {


    if (isHyst_) {
      StepTransNonLinEpsDiff();
    }
    else if (nonLin_){
      StepTransNonLin();
    }
    else if (nonLinMaterial_){
      EXCEPTION("StepTransNonLinMaterial() not supported in SolveStepPiezo");
    }

    else {
      StepTransLin();
    }

  }

  void SolveStepPiezo::StepTransMaterialNonLin(){
    
    //  std::cout<<"We do a step of transient nonlinear calculation"<<std::endl;

  }

  void SolveStepPiezo::StepTransNonLinEpsDiff() {


    //    std::cout << "\n In :StepTransNonLinEpsDiff  \n " << std::endl;
 
    bool performOneMoreStep;
    UInt iterationCounter=0;
  
    Vector<Double> newSol( numEqns_ ); 
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> solPrev( numEqns_ );
 
    // get second order derivative of previous time step;
    Vector<Double> solDeriv2Prev = PDE_.getS2();

    Double* actRHS;
    Double* solPtr;

    Vector<Double> coeff;

    oldSol.Init(0);
    newSol.Init(0);
  
    //save solution of previous time step  
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    solPrev = solHelp;
    //    std::cout << "Previous Solution:\n" << solPrev << std::endl;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Error( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin",
             __FILE__, __LINE__ );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(solPrev);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( actTime_ );    //clear RHS

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
//       if (iterationCounter == 1)
//         std::cout << std::endl << "Time step:   "  << actStep_
//                   << "  ,Iterations: " << iterationCounter << std::endl;
//       else 
//         std::cout << "Iter:  " << iterationCounter << std::endl;
    

#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear Elecs: Perform internal loop no. "
             << iterationCounter << std::endl;      
#endif
        

      // set solution of previous iteration
      if (iterationCounter == 1 ) {
        oldSol = solPrev;
      }
      else {
        oldSol = newSol;
      }
    
      //clear RHS
      algsys_->InitRHS();

      // setup RHS to incorporate loads and linear-Forms
      assemble_->AssembleLinRHS(actTime_); 

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS();

      //perform new assembly
      assemble_->AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      // K*(u_n,Ve_n)
      algsys_->UpdateRHS(STIFFNESS, solPrev.GetPointer());

      // M*(d2u_n,Ve_n)
      algsys_->UpdateRHS(MASS, solDeriv2Prev.GetPointer());

      // build in the Dirichlet vales in system mmatrix and rhs
      algsys_->BuildInDirichlet();

      //get RHS
      Vector<Double> RHS;
      algsys_->GetRHSVal( actRHS );
      StoreAlgsysToVec(RHS, actRHS );       
      Double residualNorm = PDE_.GetRhsL2Norm( RHS );
      residualNorm = 0.0;

      algsys_->SetupSolver();
      algsys_->SetupPrecond();
    
      algsys_->Solve();
      algsys_->GetSolutionVal( solPtr );
      StoreAlgsysToVec(newSol, solPtr );

      //store solution for (n+1)
      PDE_.SaveSolution( newSol.GetPointer(), newSol.GetSize() );

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
    
      // output of norms and data
      nonLinLogging_ = true;
      if ( nonLinLogging_ == true ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualNorm,
                              incrementalErr);
      }

      //    std::cout << "ResNorm=" << residualNorm << "  incrNorm=" 
      //        << incrementalErr << std::endl;
    
      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    if ( iterationCounter >= nonLinMaxIter_ ) {
      (*error) << "Number of nonlinear iterations exceeds limit "
               << "nonLinearMaxIter_ = "
               << nonLinMaxIter_;
      Error( __FILE__, __LINE__ );
    }

    //set the current values to the previous for the next time step!
    SetPreviousVals4Hyst();

    //perform corrector step  
    TS_alg_->Corrector(newSol);
  }


  void SolveStepPiezo::SetPreviousVals4Hyst() {
  
    // get solution object of PDE elec
    SinglePDE * pdeElec = domain->GetSinglePDE( "electrostatic" );

    BaseNodeStoreSol* solPDE = pdeElec->getPDESolution();
    NodeStoreSol<Double> * sol = 
          dynamic_cast<NodeStoreSol<Double>*>(solPDE);
   
    ResultInfoList results =  pdeElec->GetResultInfos();

    GradientFieldOp<Double> * FieldOp = 
      new GradientFieldOp<Double>(ptgrid_, pdeElec, pdeElec->GetEqnMap(),
                                  *sol, ELEC_POTENTIAL, 
                                  results[0], isaxi_);

    std::map<RegionIdType, BaseMaterial*> elecMat = 
      pdeElec->getPDEMaterialData();

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, Dval, dE, dD, eps;
    UInt pdeElem;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = elecMat.begin(); it != elecMat.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      Hysteresis* hyst = elecMat[actRegion]->getHysteresis();
      if ( hyst!= NULL ) {
	//get direction of polarization
	std::string str;
	elecMat[actRegion]->GetScalar(str, P_DIRECTION);
	Directions dirP;
	String2Enum(str,dirP);

	ElemList actSDList(ptgrid_ );
	actSDList.SetRegion( actRegion );
	
	EntityIterator it = actSDList.GetIterator();
	UInt iel = 0;
	for ( it.Begin(); !it.IsEnd(); it++, iel++) {
	  
	  //compute the electric field intensity
	  it.GetElem()->ptElem->GetCoordMidPoint(LCoord);
	  FieldOp->CalcElemGradField( Efield, it, LCoord, 1);
	  
	  //get correct component of electric field for scalar Preisach model
	  Ecomp   = Efield[dirP]; 	  
	  pdeElem = it.GetElem()->elemNum;
	  
	  //set the values
	  elecMat[actRegion]->SetPreviousHystVal( pdeElem, Ecomp );
	}  
      }
    }
  }


  void SolveStepPiezo::ComputeDiffEpsilon() {
  

    //we assume, that the actual solution is stored in sol_!
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    
    GradientFieldOp<Double> * FieldOp = 
      new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnMap_,
                                  *solhelp, ELEC_POTENTIAL, 
                                  results_[0],isaxi_);

    Vector<Double> LCoord, Efield;
    Double Ecomp, Pval, Dval, dE, dD, eps;
    UInt comp;
    UInt pdeElem=1;

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      ElemList actSDList(domain->GetGrid() );
      actSDList.SetRegion( subdoms_[actSD] );
      EntityIterator it = actSDList.GetIterator();
      
      //get direction of polarization
      materialData_[actSD]->GetScalar((Integer&)comp,P_DIRECTION,INTEGER);
      comp -= 1;
      UInt iel = 0;
      for ( it.Begin(); !it.IsEnd(); it++, iel++ ) {
        
        //compute the electric field intensity
        it.GetElem()->ptElem->GetCoordMidPoint(LCoord);
        FieldOp->CalcElemGradField( Efield, it, LCoord, 1);

        //get correct component of electric field for scalar Preisach model
        Ecomp = Efield[comp]; //.NormL2(); //[comp]; 

        //compute polarization
        //      pdeElem = 0;
        Pval = hyst_[actSD]->computeValue(Ecomp,pdeElem);

        //      Pval = Ecomp*2e-7;

        //compute dielectric displacement
        Dval = Pval; //8.85419E-12 * Ecomp + Pval;

        //compute differential epsilon
        dE = Ecomp - Eprevious_[pdeElem-1];
        dD = Dval - Dprevious_[pdeElem-1];
        if ( (abs(dD) < 1e-12) || (abs(dE) < 1e-10) ) {
          materialData_[actSD]->GetScalar(eps,ELEC_PERMITTIVITY,REAL);
          if (eps < 8.854e-12) {
            eps = 8.854e-12;
          }
          epsDiff_[actSD][iel] = eps;
        }
        else {
          epsDiff_[actSD][iel] = dD / dE;
        }

        pdeElem++;
      }  
    }

    //  std::cout << "EpsDiff: " << epsDiff_ << std::endl;
  }



} // end of namespace

