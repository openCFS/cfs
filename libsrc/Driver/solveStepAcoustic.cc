#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcoustic.hh"
#include "Forms/forms_header.hh"

namespace CoupledField {

SolveStepAcoustic::SolveStepAcoustic(BasePDE& apde) : BaseSolveStep(apde) {
  ENTER_FCN( "SolveStepAcoustic::SolveStepAcoustic" );
}

SolveStepAcoustic::~SolveStepAcoustic() {
  ENTER_FCN( "SolveStepAcoustic::~SolveStepAcoustic" );
}
 
// ======================================================
// Solve Step Transient SECTION  
// ======================================================

void SolveStepAcoustic::StepTransNonLin(const Integer kstep, const Double asteptime,
					const Integer level, const Boolean reset) {

  ENTER_FCN( "SolveStepAcoustic::StepTransNonLin" );

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
		   <<   "Nonlinear Acoustics: Perform internal loop no. "
		   << iterationCounter << std::endl;      
#endif
        
	// set solution of previous iteration
	oldSol = newSol;
        
	if ( iterationCounter>1 ) 
	  job = 3;
        
	// Set linear part of RHS
	algsys_->InitRHS(RhsLinVal_.GetPointer());
        
	// put nonlinear part to RHS
	// decide depending on first region
	if ( nonLinPDEName_[0] == WESTERVELT )
	  AddWesterveltRHS();
	if ( nonLinPDEName_[0] == KUZNETSOV )
	  AddKuznetsovRHS();

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

void SolveStepAcoustic::AddWesterveltRHS() {
  
  ENTER_FCN( "SolveStepAcoustic::AddWesterveltRHS" );

  Matrix<Double>     ptCoord;
  Vector<Double>     sol, solderiv1, solderiv2, rhs;
  BaseFE             * ptElem;
  StdVector<Integer> connect, connect_PDE;
  
  Double coeff;
  Double density, compressibility, c0, BoverA;

  BaseForm * rhsInt = new  nLinWesterveltRHSInt(1.0, isaxi_);

  for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {
	StdVector<Elem*> elemssd;
	ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
        
	// get material data
	density         = materialData_[actSD].GetDensity();
	compressibility = materialData_[actSD].GetCompressibility();
	c0 = sqrt(compressibility/density);
	BoverA = materialData_[actSD].GetBoverA();

	// set correct factors for bilinear forms
	coeff = (1+0.5*BoverA) / pow(c0,4);
	rhsInt->SetFactor(coeff);
        
	for (Integer j=0; j < elemssd.GetSize(); j++) {

	  ptElem  = elemssd[j]->ptElem;
	  connect = elemssd[j]->connect;
 	  GetElemCoords(connect, ptCoord, actlevel_);

	  GetSolVecOfElement(sol, connect);
	  GetDerivSolVecOfElement(solderiv1, connect);
	  GetDeriv2SolVecOfElement(solderiv2, connect);

	  rhsInt->SetElemPtr(ptElem);

	  rhsInt->SetActElemSol(sol);
	  rhsInt->SetActElemSolDeriv1(solderiv1);
	  rhsInt->SetActElemSolDeriv2(solderiv2);
	  rhsInt->CalcElemVector(ptCoord, rhs);

	  //get equation numbers 
	  eqnData_->Node2EQN(connect, connect_PDE);

	  //assemble
	  algsys_->SetElementRHS(&rhs[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
	}  
  }

}

void SolveStepAcoustic::AddKuznetsovRHS() {
  
  ENTER_FCN( "SolveStepAcoustic::AddKuznetsovRHS" );

  Matrix<Double>     ptCoord;
  Vector<Double>     sol, solderiv1, solderiv2, rhs;
  BaseFE             * ptElem;
  StdVector<Integer> connect, connect_PDE;
  
  Double coeff1, coeff2;
  Double density, compressibility, c0, BoverA;

  BaseForm * rhsInt = new  nLinKuznetsovRHSInt(1.0, isaxi_);
  //   BaseForm * N1 = new nLinAcoustic1(coeffN1, isaxi_);
  //   BaseForm * N2 = new nLinAcoustic2(coeffN2, isaxi_);

  for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {
	StdVector<Elem*> elemssd;
	ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
        
	// get material data
	density         = materialData_[actSD].GetDensity();
	compressibility = materialData_[actSD].GetCompressibility();
	c0 = sqrt(compressibility/density);
	BoverA = materialData_[actSD].GetBoverA();

	// set correct factors for bilinear forms
	coeff1 = density * BoverA / pow(c0,4);
	rhsInt->SetFactor(coeff1);
	coeff2 = density * 2.0 / (c0*c0);
	rhsInt->SetSecondFactor(coeff2);
        
	for (Integer j=0; j < elemssd.GetSize(); j++) {

	  ptElem  = elemssd[j]->ptElem;
	  connect = elemssd[j]->connect;
 	  GetElemCoords(connect, ptCoord, actlevel_);

	  GetSolVecOfElement(sol, connect);
	  GetDerivSolVecOfElement(solderiv1, connect);
	  GetDeriv2SolVecOfElement(solderiv2, connect);

	  rhsInt->SetElemPtr(ptElem);

	  rhsInt->SetActElemSol(sol);
	  rhsInt->SetActElemSolDeriv1(solderiv1);
	  rhsInt->SetActElemSolDeriv2(solderiv2);
	  rhsInt->CalcElemVector(ptCoord, rhs);
          
	  //        N1->SetElemPtr(ptElem);
	  //        N1->SetActElemSolDeriv1(solderiv1);
	  //        N1->CalcElementMatrix(ptCoord, elemmat);
	  //        rhs = elemmat * solderiv2;
	  //        N2->SetElemPtr(ptElem);
	  //        N2->SetActElemSol(sol);
	  //        N2->CalcElementMatrix(ptCoord, elemmat);
	  //        rhs += elemmat * solderiv1;

	  //get equation numbers 
	  eqnData_->Node2EQN(connect, connect_PDE);

	  //assemble
	  algsys_->SetElementRHS(&rhs[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
	}  
  }

}


} // end of namespace

