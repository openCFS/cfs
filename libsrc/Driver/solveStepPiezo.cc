#include <fstream>
#include <iostream>
#include <string>

#include "solveStepPiezo.hh"
#include "Forms/forms_header.hh"
#include "Utils/preisach.hh"

namespace CoupledField {

SolveStepPiezo::SolveStepPiezo(BasePDE& apde) : BaseSolveStep(apde) {
  ENTER_FCN( "SolveStepPiezo::SolveStepPiezo" );

  doInit_ = TRUE;
}

SolveStepPiezo::~SolveStepPiezo() {
  ENTER_FCN( "SolveStepPiezo::~SolveStepPiezo" );
}
 
// ======================================================
// Solve Step Transient SECTION  
// ======================================================

void SolveStepPiezo::PreStepTrans( const Integer kstep, const Double asteptime,
				   const Integer level, const Boolean reset ) {

  ENTER_FCN( "SolveStepPiezo::PreStepTrans" );

  lasttimecalc_= asteptime;
  
  // due to coupling-pdes, the RHS has to be initialized BEFORE 
  // the coupling forces are assembled to the RHS
  algsys_->InitRHS();

  if (doInit_) {
    //set the Preisach values; should be obtained from the xml_file
    Esat_     = 2.0e6;
    Psat_     = 0.04;
    Ec_       = 0.9e6;
    Dir_      = 2;
    isVirgin_ = TRUE;

    Integer numElems = PDE_.getPDE_numElems();
    hyst_ = new Preisach(numElems, Esat_, Psat_, Ec_, isVirgin_);

    doInit_ = FALSE;
  }
  else {
    DoUpdateHyst();
  }
}
  
void SolveStepPiezo::StepTransNonLin(const Integer kstep, const Double asteptime,
				     const Integer level, const Boolean reset) {

  ENTER_FCN( "SolveStepPiezo::StepTransNonLin" );

  laststepcalc_ = kstep;
  lasttimecalc_ = asteptime;

  Integer job;
  Boolean performOneMoreStep;
  Integer iterationCounter=0;
  
  Vector<Double> newSol(numPDENodes_);
  Vector<Double> oldSol(numPDENodes_);
  
  // just update dirichlet values
  job = 3;
  
  // if first time step, setup system matrix
  if (laststepcalc_ == 1) {
    assemble_->AssembleMatrices(level);
    algsys_->ConstructEffectiveMatrix(matrix_factor_);
    
    //set job to 1: build in dirichlet BCs and compute preconditioner
    job = 1;
  }  
  
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

  //compute predictors
  TS_alg_->Predictor(solhelp->GetAlgSysVector());

  // set BCs, if effective mass matrix formulation, values of BCs depend on 
  //  predictors, so predictors have to be computed beforehand
  SetBCs(level, lasttimecalc_);

  // set old solution  
  newSol = solhelp->GetAlgSysVector();

  // Update RHS (mass matrix and damping matrix on right hand side)
  TS_alg_->UpdateRHS();
  // stores this as linear part of RHS
  StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

  do {
    iterationCounter++;
    // for every time step write out number of iteration loops to standard out
    if (iterationCounter == 1)
      std::cout << std::endl << "Time step:   "  << kstep 
		<< "  ,Iterations: " << iterationCounter;
    else 
      std::cout	<< "  " << iterationCounter;
    
#ifdef DEBUG
    *debug << std::endl
	   << "====================================================== "
	   << std::endl
	   <<   "Nonlinear Piezos: Perform internal loop no. "
	   << iterationCounter << std::endl;      
#endif
        
    // set solution of previous iteration
    oldSol = newSol;
    
    if ( iterationCounter>1 ) 
      job = 3;
    
    // Set linear part of RHS
    algsys_->InitRHS(RhsLinVal_.GetPointer());
    
    // put nonlinear part to RHS
    AddPolarizationToRHS();
    
    algsys_->BuildInDirichlet();
    if ( job == 1 ) {
      algsys_->SetupSolver(job);
      algsys_->SetupPrecond(job);
    }
    
    algsys_->Solve();
    
    // store new solution in newSol
    StoreAlgsysToVec(newSol, algsys_->GetSolutionVal() );
    
    // perform corrector step, if effective mass formulation is used,
    //   we need the Corrector step before we store newsol to sol_,
    //   because newsol is second time derivative at first!
    TS_alg_->Corrector(newSol);
    
#ifdef DEBUG
    *debug << std::endl
	   << "New Solution:" << std::endl << newSol << std::endl;
#endif

    //put new solution to sol_
    sol_->SetAlgSysVector(newSol);  
    
    // compute L2-Norm of error between last incremental solution and
    //   actual incremental solution
    Double solIncrL2Norm=0;
    for (Integer i=0; i<newSol.GetSize(); i++)
      solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
    
    solIncrL2Norm = sqrt(solIncrL2Norm);
    Double actSolL2Norm = newSol.NormL2();
    
    Double incrementalErr;
    if (actSolL2Norm > 1)
      incrementalErr = solIncrL2Norm / actSolL2Norm;
    else
      incrementalErr = solIncrL2Norm;
    
    // output of norms and data
    if ( nonLinLogging_ == TRUE ) {
      Info->WriteNonLinIter(pdename_, iterationCounter, incrementalErr,
			    incrementalErr);
    }
    
    // boolean variable, holds condition if another iteration step
    //   is necessary
    performOneMoreStep =    (incrementalErr > incStopCrit_);
    
  } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
  
  
}


void SolveStepPiezo::AddPolarizationToRHS() {
  
  ENTER_FCN( "SolveStepPiezo::AddPolarizationToRHS" );

  //we assume, that the actual solution is stored in sol_!
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

  GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
								  *solhelp, ELEC_POTENTIAL, 
								  actlevel_, isaxi_);

  Vector<Double> LCoord, Efield;
  Double Ecomp, Pval;
  Integer comp = Dir_ - 1;
  Integer pdeElem=1;

  Matrix<Double>     ptCoord;
  Vector<Double>     sol, solderiv1, solderiv2, rhs;
  BaseFE             *ptElem;
  StdVector<Integer> connect, connect_PDE;
  
  BaseForm * rhsInt = new  PiezoPolarizationInt(Dir_, PDE_.getPDE_dofspernode(), isaxi_);

  for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {
    StdVector<Elem*> elemssd;
    ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
    
    for (Integer iel=0; iel < elemssd.GetSize(); iel++) {

      //compute the electric field intensity
      elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
      FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);
      
      //get correct component of electric field for scalar Preisach model
      Ecomp = Efield[comp]; 

      //compute polarization
      hyst_-> computeValue(Ecomp,pdeElem);
     
      ptElem  = elemssd[iel]->ptElem;
      connect = elemssd[iel]->connect;
      GetElemCoords(connect, ptCoord, actlevel_);
      
      rhsInt->SetElemPtr(ptElem);
      rhsInt->SetFactor(Pval);

      rhsInt->CalcElemVector(ptCoord, rhs);
      
      //get equation numbers 
      eqnData_->Node2EQN(connect, connect_PDE);
      
      //assemble
      algsys_->SetElementRHS(&rhs[0], connect_PDE.GetPointer(), connect_PDE.GetSize());

      pdeElem++;
    }  
  }

}


void SolveStepPiezo::DoUpdateHyst() {
  
  ENTER_FCN( "SolveStepPiezo::DoUpdateHyst" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

  GradientFieldOp<Double> * FieldOp = new GradientFieldOp<Double>(ptgrid_, &PDE_, eqnData_,
								  *solhelp, ELEC_POTENTIAL, 
								  actlevel_, isaxi_);

  Vector<Double> LCoord, Efield;
  StdVector<Elem*> elemssd;
  Integer pdeElem=1;
  Double Ecomp;
  Integer comp = Dir_-1;
  
  // loop over all subdomains
  for (Integer isd=0; isd<subdoms_.GetSize(); isd++) {
    // get vector of Elem of subdomain with color: subdoms[isd]
    ptgrid_->GetElemSD(elemssd,subdoms_[isd],actlevel_);
      
    // loop over elements of subdomain
      for (Integer iel=0; iel< elemssd.GetSize(); iel++)	{
	elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);

	//compute electric field
	FieldOp->CalcElemGradField( Efield, elemssd[iel], LCoord, 1);

	//get correct component of electric field for scalar Preisach model
	//and invoke the update MinMaxList method
	Ecomp = Efield[comp]; 
	hyst_->updateMinMaxList(Ecomp, pdeElem);

	pdeElem++;
      }
  }

  delete FieldOp;

}

} // end of namespace

