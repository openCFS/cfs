#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"
#include "Utils/mathfunctions.hh"

namespace CoupledField {

// =====================================================================
// set solution information
// =====================================================================
AcousticPDE::AcousticPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
						 FileType *aptFileType, WriteResults *aptOut)
  :BasePDE(aptgrid,aptbcs,aptFileType,aptOut,aptTimeFunc) {

  ENTER_FCN( "AcousticPDE::AcousticPDE" );

  dofspernode_ = 1;
  solTypes_ = ACOU_POTENTIAL;
  solDofs_ = 1;
  pdename_          = "acoustic";
  pdematerialclass_ = "fluid";

  coarsealpha_ = 0.01;

  laststepcalc_ = 0;
    
  isaxi_ = params->HasValue( "type", "axi", "geometry" );

  // *********************************************
  //   Check what type of damping should be used
  // *********************************************

  dampingType_ = NONE;
  params->GetList( "type", dampingList_, pdename_, "damping");

  if ( dampingList_.GetSize() == 0 )
	Info->PrintF( pdename_, "\nNo Damping information detected!\n" );

  Boolean sorted_ = TRUE;
  for ( Integer k = 0; k < dampingList_.GetSize(); k++) {
	// Rayleigh damping
	if (dampingList_[k] == "rayleigh") {
	  dampingType_ = RAYLEIGH;
	  Info->PrintF( pdename_, "\nUsing RAYLEIGH damping for region: %d\n", k );
	}
	// Thermoviscous damping
	else if (dampingList_[k] == "thermoViscous") {
	  dampingType_ = THERMOVISCOUS;
	  Info->PrintF(pdename_,"\nUsing THERMOVISCOUS damping for region: %d\n", k );
	}
	// Fractional damping
	else if (dampingList_[k] == "fractional") {
	  dampingType_ = FRACTIONAL;
	  Info->PrintF( pdename_, "\nUsing FRACTIONAL damping for region: %d\n", k );
	}
	else if (dampingList_[k] == "no") {
	  dampingType_ = NONE;
	  Info->PrintF( pdename_, "\nUsing NO damping at all for region: %d\n", k );
	}
	if ( k > 1 )
	  if ( dampingList_[k] != dampingList_[k-1] )
		sorted_ = FALSE;
  }
  if ( sorted_ == FALSE )
	Error("Specify same type of damping for all regions!",__FILE__,__LINE__);

  // get additional information for fractional damping model
  if ( dampingType_ == FRACTIONAL ) {

	StdVector<std::string> fracAlgList_;
	params->GetList( "fracAlg", fracAlgList_, pdename_, "damping" );
	StdVector<Integer> fracMemoryList_;
	params->GetList( "fracMemory", fracMemoryList_, pdename_, "damping" );
	StdVector<std::string> interpolationList_;
	params->GetList( "interpolation", interpolationList_, pdename_, "damping" );

	if( fracAlgList_.GetSize() == 0 || fracMemoryList_.GetSize() == 0
		|| interpolationList_.GetSize() == 0 )
	  Error("Specify attributes fracAlg, fracMemory and interpolation!",__FILE__,__LINE__);
	// up to now take values from first subdomain
	else {
	  if ( fracAlgList_[0] == "gl" )
		Info->PrintF( pdename_, "  with Gruenwald-Letnikov algorithm,\n");
	  else if (fracAlgList_[0] == "blank")
		Info->PrintF( pdename_, "  with Blanks algorithm,\n");

	  fracMemory_ = fracMemoryList_[0];
	  Info->PrintF( pdename_, "  memory size is: %d,\n", fracMemory_);

	  if ( interpolationList_[0] == "lin1pt")
		inType_ = LIN1PT;
	  else
		inType_ = NOTUSED;
	  Info->PrintF( pdename_, "  %s interpolation of past values\n\n"
					, interpolationList_[0].c_str() );
	}
	// modify dampingList, so that fracAlg is included
	for ( Integer k = 0; k < dampingList_.GetSize(); k++) {
	  if ( fracAlgList_[k] == "gl" )
		dampingList_[k] = "fractional_gl";
	  else if (fracAlgList_[k] == "blank")
		dampingList_[k] = "fractional_blank";
	}
  }

  // *************************************************************
  //   Check what type of nonlinear PDE formulation should be used
  // *************************************************************

  nonLin_ = FALSE; //declaration in basePDE.hh
  params->GetList( "nonLinear", nonLinPDEType_, pdename_, "region" );

  for ( Integer k = 0; k < nonLinPDEType_.GetSize(); k++ ) {
	if ( nonLinPDEType_[k] != "no" ) {
	  nonLin_ = TRUE;
	  break;
	}
  }

  if( nonLin_ ) {
	// solution method
	params->Get("method", nonLinMethod_, pdename_, "nonLinear" );
	// perform logging?
	nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );
	// type of line search
	params->Get("type", lineSearch_, pdename_, "lineSearch" );
	// incremental stopping criterion
	params->Get("incStopCrit", incStopCrit_, pdename_, "nonLinear" );
	// residual stopping criterion
	params->Get("resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
	// maximal number of NL-iterations
	params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
  }

  // ***************************************************************
  //   If no other damping type is specified and we have absorbing
  //   boundary conditions, then use ABCDAMP
  // ***************************************************************
  absorbingBCs_ = FALSE;
  params->GetList( "name", absBCs_, pdename_, "absorbingBCs" );
  if ( absBCs_.GetSize() ) {
	absorbingBCs_ = TRUE;
	Info->PrintF( pdename_, " Apply Absorbing Boundary Conditions\n" );
	surfdoms_ = absBCs_;
  }

  needsDampingMatrix_ = FALSE;
  // two matrices formulation
  // if ( absorbingBCs_ == TRUE || (dampingType_ != NONE  && dampingType_ != FRACTIONAL) )
  if ( absorbingBCs_ == TRUE || dampingType_ != NONE )
	needsDampingMatrix_ = TRUE;  
}


void AcousticPDE::DefineIntegrators(const Integer level) {

  ENTER_FCN( "AcousticPDE::DefineIntegerators" );

  Boolean nonLin = FALSE;
  Double density, compressibility, c0, alpha0, y, BoverA;
  Double coeffmass, coeffdamp;

  for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

	density         = materialData_[actSD].GetDensity();
	compressibility = materialData_[actSD].GetCompressibility();
	c0 = sqrt(compressibility/density);

	alpha0 = materialData_[actSD].GetDampingAlfa();
	y      = materialData_[actSD].GetDampingBeta();

	// **********************************************************************
	//   Linear wave equation
	// **********************************************************************

	// stiffness integrator
	BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);	  
	IntegratorDescriptor * stiffIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);
	assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);

	// mass integrator
	coeffmass = density / (c0*c0);
	//coeffmass  = density*density/compressibility;
	BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
	IntegratorDescriptor * massIntDescr = new IntegratorDescriptor(bilinearMass, MASS);
	assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]);

	// **********************************************************************
	//   Additional terms for damping in linear wave equation
	// **********************************************************************

	if ( nonLinPDEType_[actSD]=="no" && dampingList_.GetSize() != 0) {

	  if (dampingList_[actSD] == "rayleigh") {
		// This works even after assemble_->AddIntegrator() is executed because of the pointers...
		// stiffness part
		stiffIntDescr->SetSecondaryMat(DAMPING, materialData_[actSD].GetDampingBeta(), analysistype_);
		
		// mass part
		massIntDescr->SetSecondaryMat(DAMPING, materialData_[actSD].GetDampingAlfa(), analysistype_);
	  }
	  
	  else if ( dampingList_[actSD] == "thermoViscous" ) {
		coeffdamp  =  density * 2.0 * alpha0 * c0;
		BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);	
		IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearStiff, DAMPING);
		
		assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	  }

	  else if ( dampingList_[actSD] == "fractional_gl" ) {
 		coeffdamp = - density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
		BaseForm * bilinearDamp  = new MassInt(coeffdamp, dofspernode_, isaxi_);
		bilinearDamp->SetFracDamping();
		// conservative formulation
		// adapt NewmarkFracDamp::Init and BasePDE::GetFracDampMatrixCoeff
		// IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, DAMPING);

		// two matrices formulation
		IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, STIFFNESS);

		assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	  }

	  else if  ( dampingList_[actSD] == "fractional_blank" ) {
 		coeffdamp =  - density * 2.0 * alpha0 / c0 / sin((y-1.0)*PI/2.0);
		coeffdamp *= exp(-gammaln(1.0- (y- 1.0)) ); // prefactor of blank alg
		coeffdamp *= 1.0/(1.0- (y- 1.0));           // weight factor of index 0
		BaseForm * bilinearDamp  = new MassInt(coeffdamp, dofspernode_, isaxi_);
		bilinearDamp->SetFracDamping();
		// conservative formulation
		// adapt NewmarkFracDamp::Init and BasePDE::GetFracDampMatrixCoeff		
		// IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, DAMPING);

		// two matrices formulation
		IntegratorDescriptor * dampIntDescr = new IntegratorDescriptor(bilinearDamp, STIFFNESS);

		assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
	  }

	}

	// **********************************************************************
	//   Additional terms for nonlinear wave equation
	// **********************************************************************
      
	else if ( nonLinPDEType_[actSD] == "kuznetsov" ) {
	
	  BoverA = materialData_[actSD].GetBoverA();
	  Double coeffN1 = BoverA * pow(compressibility/density,2);
	  BaseForm * N1 = new nLinAcoustic1(isaxi_);
	  assemble_->AddRhsIntegrator(N1, subdoms_[actSD], nonLin_);
	
	  Double coeffN2 = 2 * compressibility/density;
	  BaseForm * N2 = new nLinAcoustic2(isaxi_);
	  assemble_->AddRhsIntegrator(N2, subdoms_[actSD], nonLin_);
	
	  if ( dampingType_ == THERMOVISCOUS ) {
		alpha0 = materialData_[actSD].GetDampingAlfa();
		coeffdamp = 1.0;
		stiffIntDescr->SetSecondaryMat(DAMPING, coeffdamp, analysistype_);
		assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);
	  }
	
	}
	else if ( nonLinPDEType_[actSD] == "westervelt" ) {
	  alpha0 = materialData_[actSD].GetDampingAlfa();
	  coeffdamp = 1.0;
	  stiffIntDescr->SetSecondaryMat(DAMPING, coeffdamp, analysistype_);
	  assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);
	
	  BoverA = materialData_[actSD].GetBoverA();
	  Double coeffN1 = (1.0+0.5*BoverA)*pow(compressibility/density,2)/density;
	  BaseForm * N1 = new nLinAcoustic1(isaxi_);
	  assemble_->AddRhsIntegrator(N1, subdoms_[actSD], nonLin_);  
	}
  }

  // **********************************************************************
  //   surface-integration: Absorbing boundaries
  // **********************************************************************

  if ( absorbingBCs_ == TRUE && analysistype_ != HARMONIC ) {
	for (Integer actSD = 0; actSD < absBCs_.GetSize(); actSD++) {
	  //currently hard-coded!!
	  Double density = materialData_[0].GetDensity();
	  Double compressibility = materialData_[0].GetCompressibility();
	  Double coeffdamp = density/((sqrt(compressibility/density)));
	
	  BaseForm * bilinear_damp = new MassInt(coeffdamp,dofspernode_, isaxi_);
	  assemble_->AddSurfIntegrator(bilinear_damp,  absBCs_[actSD],
								   DAMPING, nonLin);
	}
  }
 
}

// ======================================================
// SOLVING SECTION
// ======================================================

void AcousticPDE::InitTimeStepping() {
  ENTER_FCN( "AcousticPDE::InitTimeStepping" );
    
  if ( dampingType_ == FRACTIONAL ) {
	TS_alg_ = new NewmarkFracDamp(pdename_, algsys_, eqnData_, ptgrid_, this,
								  subdoms_, dampingList_, 
								  fracMemory_, inType_, isaxi_);
  }
  else {
	// this includes rayleigh and thermoviscous damping
	TS_alg_ = new Newmark(pdename_, algsys_, eqnData_, needsDampingMatrix_);
  }

  // Needed for fractional damping model.
  // See Assemble::AssembleMatrices
  assemble_->SetPDEPointer(this);

}

void AcousticPDE::StepTransNonLin(const Integer kstep, const Double asteptime,
 								  const Integer level, const Boolean reset) {

  ENTER_FCN( "AcousticPDE::StepTransNonLin" );

//    lasttimecalc_ = asteptime;
//    laststepcalc_ = kstep;


//    const Integer job = 1;
  
//    static Integer timeStepCounter=1;
//    Boolean performOneMoreStep;
//    Integer iterationCounter=0;

//    Vector<Double> actSol(numPDENodes_);
//    Vector<Double> solInc(numPDENodes_);
  
//    // Cast BaseStoreSol into StoreSol<Double>,
//    // since this function is only called
//    // in the transient case
//    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

//    //set actual solution  
//    actSol = solhelp->GetAlgSysVector();

//    //compute predictors
//    TS_alg_->Predictor(solhelp->GetAlgSysVector());

//    // inner forces due to nonlin formulation
//    assemble_->AssembleNLRHS(level, lasttimecalc_);  

//    //Update RHS (mass matrix on right hand side)
//    TS_alg_->UpdateRHS(solhelp->GetAlgSysVector());
  
//    timeStepCounter++;
//    do {
// 	  iterationCounter++;
// 	  std::cout << std::endl
// 				<< "Nonlinear Acoustics: Perform internal loop no. " 
// 				<< iterationCounter << std::endl;

// #ifdef DEBUG
// 	  *debug << std::endl
// 			 << "====================================================== "
// 			 << std::endl
// 			 <<	"Nonlinear Acoustics: Perform internal loop no. "
// 			 << iterationCounter << std::endl;      
// #endif

// 	  // setup and solve new system (rhs is already set) ====================
// 	  assemble_->InitNonLinMatrices();
// 	  assemble_->AssembleMatrices(level);
// 	  algsys_->ConstructEffectiveMatrix(matrix_factor_);

// 	  SetBCs(level, lasttimecalc_);

// #ifdef USE_OLAS
//       algsys_->BuildInDirichlet();
//       algsys_->SetupPrecond(job);
// #else
//       algsys_->CalcPrecond(job);
// #endif

// 	  algsys_->Solve();

// 	  // new solution is only an increment of the full solution =============
// 	  StoreAlgsysToVec(solInc, algsys_->GetSolutionVal() );

// 	  Double residualL2Norm;
//       Double etaLineSearch = 0;
      
//       if ( lineSearch_ != "no" ) {
// 		actSol += solInc;
//       }
//       else {
		
// 		// TRUE is for transient simulation
// 		residualL2Norm = LineSearch(solInc, actSol, etaLineSearch, level,TRUE);
//       }

// 	  //store A_/n+1) in the solution-object sol_
// 	  sol_->SetAlgSysVector(actSol);

// 	  // recalculate RHS with new values to get new residual (f^(k+1))=======
// #ifndef USE_OLAS    
// 	  algsys_->InitRHS(RhsLinVal_.GetPointer());
// #endif

// 	  //Update RHS (mass matrix on right hand side)
// 	  TS_alg_->UpdateRHS(actSol);

// 	  // inner forces due to nonlin formulation
// 	  assemble_->AssembleNLRHS(level, lasttimecalc_);
 

// 	  // ====================================================================
// 	  // calculation of error norms
// 	  // ====================================================================

//       if ( lineSearch_ != "no" ) {
		
// 		Vector<Double> actRHS;
// 		StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
		
// 		// ------------------------------------------------------------------
// 		// calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
// 		// ------------------------------------------------------------------
// 		residualL2Norm = RhsL2Norm(actRHS);
//       }
	  
// 	  Double residualErr;
// 	  // Please adapt!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
// 	  if (1)
// 		residualErr    = residualL2Norm;
// 	  else
// 		residualErr    = residualL2Norm;
	  
// 	  // --------------------------------------------------------------------
// 	  // calculate incremental error
// 	  // --------------------------------------------------------------------
// 	  Double solIncrL2Norm = solInc.NormL2();
// 	  Double actSolL2Norm = actSol.NormL2();
// 	  Double incrementalErr;
      
// 	  if (actSolL2Norm > 1)
// 		incrementalErr = solIncrL2Norm / actSolL2Norm;
// 	  else
// 		incrementalErr = solIncrL2Norm;
	  
// 	  // --------------------------------------------------------------------
// 	  // output of norms and data
// 	  // --------------------------------------------------------------------
// 	  if ( nonLinLogging_ == TRUE ) {
// 		Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
// 							  incrementalErr);
// 	  }
	  
// 	  // boolean variable, holds condition if another iteration step
// 	  // is necessary
// 	  performOneMoreStep = 
// 		(incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
      
//    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

//    //perform corrector step  
//    TS_alg_->Corrector(actSol);
}


// calculates L2-norm of RHS regarding dirichlet entries due to penalty
// formulation by setting them 0
Double AcousticPDE::RhsL2Norm(Vector<Double>& actRHS) {
  ENTER_FCN( "AcousticPDE::RhsL2Norm" );
    
  Integer node, dof;
  StdVector<Integer> eqn(1);
  std::list<Integer> nodes;
    
  // Eliminate dirichlet node from RHS (due to penalty formulation)
  for (Integer i=0; i< bcs_hd_.GetSize(); i++) {

	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end();
		 p++) {
	  node=*p;
	  eqnData_->Node2EQN(node,eqn);
	  if (eqn[0] != 0){
		actRHS[(eqn[0]-1)] = 0;
	  }
	}
  }
  return actRHS.NormL2();
}

Double AcousticPDE::LineSearch( Vector<Double>& solInc, Vector<Double>& actSol, 
								Double& etaLineSearch, Integer level,
								Boolean trans ) {
  ENTER_FCN( "AcousticPDE::LineSearch" );

  //  Vector<Double> solOld(actSol);
  const Integer nrEtas = 4;
  const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125};
  Double etaOpt;
  Double residualL2NormOpt = 1e15;
  
  actSol += solInc;
  
}

// ======================================================
// COUPLING SECTION
// ======================================================

void AcousticPDE::InitCoupling(PDECoupling * Coupling) {

  ENTER_FCN( "AcousticPDE::InitCoupling" );

  pdeIsCoupled_ = TRUE;
  ptCoupling_   = Coupling;
  
  // Intialize the memory of the coupling values
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
	if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE)	{
	  ptCoupling_->CreateCouplingVector(i,isComplex_);

	  //now since we need a incremental formulation, initialize some necessary vectors
	  isIncrFormulation_ = TRUE;
	  solIncr_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
	  actSol_.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());
	  solIncr_.Init(0);
	  actSol_.Init(0);
	}
  }
}


void AcousticPDE::CalcOutputCoupling() {

  ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

  Integer dof;
  SolutionType quantity;
  StdVector<Elem*> * couplingElems = NULL;
  StdVector<Elem*> * interfaceVolElems = NULL;
  StdVector<Integer> * couplingNodes = NULL;
  StdVector<MaterialData*> * couplingMaterials = NULL;
  CFSVector * temp_values = NULL;
  
  // loop over all output coupling quantities
  for (Integer i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
	quantity = ptCoupling_->GetOutputQuantity(i);
	ptCoupling_->GetOutputValues(i, temp_values);
      
	Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

	switch(ptCoupling_->GetOutputType(i)) {

	case NODE:
	  if (quantity == ACOU_FORCE) {
		ptCoupling_->GetOutputElements(i, couplingElems);
		ptCoupling_->GetOutputNodes(i, couplingNodes);
		ptCoupling_->GetOwnMaterials(i, couplingMaterials);
		ptCoupling_->GetInputNeighbourElems(i, interfaceVolElems);
		dof = ptCoupling_->GetOutputDof(i);
	    
		CalcMechCouplingRHS(couplingElems, *couplingNodes,
							couplingMaterials, *values, dof,
							interfaceVolElems);
	  }
	  break;

	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
  }
}


void AcousticPDE::CalcMechCouplingRHS(StdVector<Elem*> * couplingElems, 
									  StdVector<Integer> & couplingNodes,
									  StdVector<MaterialData*> * couplingMaterials,
									  Vector<Double>& elemCouplingSols,
									  Integer couplingdof,
									  StdVector<Elem*> * interfaceVolElems)
{
  ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );

  Double density=0;
   
  elemCouplingSols.Init(0.0);

  for (Integer actElem=0; actElem<couplingElems->GetSize(); actElem++) {
	Elem * actCoupleElem = (*couplingElems)[actElem];
	BaseFE * ptElem          = actCoupleElem->ptElem;

	StdVector<Integer> connecth = actCoupleElem->connect;

	Matrix<Double> ptCoord; 
	GetElemCoords(connecth, ptCoord, actlevel_);

	Boolean found = FALSE;
      
	// get correct density belonging to the neighbouring element
	// of the interface
	for (Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++)	{  
	  if ((*interfaceVolElems)[actElem]->namesd ==  subdoms_[actSD]) {
		density = materialData_[actSD].GetDensity();
		found = TRUE;
	  }
	}
      
	if (found ==FALSE) {
	  std::string msg = "Cannot find correct density to compute acoustic ";
	  msg += "pressure forces! Using density of acoustic subdomain 1";
	  Info->Warning( msg );
	  density = materialData_[0].GetDensity();
	}
          
	BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
	Matrix<Double> elemmat;
	bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
	delete bilinear_mass;	  
	Vector<Double> sol;
	GetDerivSolVecOfElement(sol, connecth);
      
	Vector<Double> forceOnElem = elemmat * sol;
	// force has to be added on RHS with negative sign
	forceOnElem *= -1;

	// the normal vector points outwards of the MECHANICAL domain
	// (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
	Vector<Double> n;

	// points outward own domain
	ptgrid_->CalcSurfNormalOutOfVol(n, *actCoupleElem, *(*interfaceVolElems)[actElem]); 
	n *=-1;
	//std::cerr << "acoustic Normal =\n" << n << std::endl;

	for (Integer actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
	  Integer nodePos = 0;
	  
	  while(connecth[actNode] != couplingNodes[nodePos] &&
			nodePos < couplingNodes.GetSize()) {
		nodePos++;
	  }
	  
	  for (Integer actDof=0; actDof < couplingdof ; actDof++) {
		elemCouplingSols[nodePos*couplingdof +actDof] += 
		  forceOnElem[actNode] * n[actDof];
		//std::cerr << "acouforceonElem +=" << forceOnElem[actNode] * n[actDof] << std::endl;
	  }
	}
  }
}


Boolean AcousticPDE::HasOutput(SolutionType output) {
  ENTER_FCN( "AcousticPDE::HasOutput" );
  if (output == ACOU_FORCE) {
	return TRUE;
  }
  return FALSE;
}


// ======================================================
// POSTPROCESSING SECTION
// ======================================================

void AcousticPDE::WriteResultsInFile(Integer stepOffset,
									 Double timeOffset) {
  ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
  int commrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
  if (!commrank) {
#endif
	NodeStoreSol<Double> solIm_mesh;
	NodeStoreSol<Double> * solTransient;
	NodeStoreSol<Complex> * solHarmonic;
       
	Double actTime = lasttimecalc_ + timeOffset;
	Integer actStep = laststepcalc_ + stepOffset;
      
	if (analysistype_==HARMONIC) {
	  solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

	  if (saveSol_)
		outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actFreqStep_, 
											actFrequency_, complexFormat_);
	  if (saveSolHist_)
		outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actFreqStep_, 
										   actFrequency_, complexFormat_);
	}
	else {  
	  solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

	  if (saveSol_){
		outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
	  
		if (saveSolHist_)
		  outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
	  }
	
	  if (saveDeriv1_) {
		solDeriv1_.SetAlgSysVector(getS1()); 
		outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
	
		if (saveDeriv1Hist_)
		  outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
	  }

	  if (saveDeriv2_) {
		solDeriv2_.SetAlgSysVector(getS2());
		outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
	  
		if (saveDeriv2Hist_)
		  outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
	  }
	}
#ifdef PARALLEL
  }//!commrank
#endif
}


// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void AcousticPDE::ReadStoreResults() {  
  ENTER_FCN( "AcousticPDE::ReadStoreResults" );
    
  // Construct vectors for restricted parameter search
  StdVector<std::string> keyVec;
  StdVector<std::string> attrVec;
  StdVector<std::string> valVec;
  std::string quantity;
    
  // *****************************
  // Determine nodal results
  // ***************************** 
  StdVector<std::string> nodeValues;
  keyVec  = pdename_, "storeResults", "nodeResults", "region";
  attrVec = "", "", "type";  
    
  // --- Acou Potential, 1. Deriv ---
  Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
  valVec = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, nodeValues);
  if (nodeValues.GetSize() > 0) {
	saveDeriv1_ = TRUE;
      
	// intialize corresponding storesolution object
	solDeriv1_.SetNumSolutions(1);
	solDeriv1_.SetNumNodes(numPDENodes_);
	solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
	solDeriv1_.SetNumDofs(1);
	solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
	solDeriv1_.Init(0);
  }

  // --- Acou Potential, 2. Deriv ---
  Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
  valVec = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, nodeValues);
  if (nodeValues.GetSize() > 0) {
	saveDeriv2_ = TRUE;
      
	// intialize corresponding storesolution object
	solDeriv2_.SetNumSolutions(1);
	solDeriv2_.SetNumNodes(numPDENodes_);
	solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
	solDeriv2_.SetNumDofs(1);
	solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_); 
	solDeriv2_.Init(0);
  }
    
  // *****************************
  // Determine element results
  // *****************************
    
  // --- nothding to do here ---

  // *****************************
  // Determine nodal history
  // *****************************
  StdVector<std::string> saveNodeHist; 
  keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
  attrVec = "", "", "type";
    
  // --- Acoustic Potential  ---
  Enum2String(ACOU_POTENTIAL, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
  if (saveNodeHist.GetSize() > 0) {
	if (saveSol_ == FALSE) {
	  std::string errMsg = pdename_;
	  errMsg += ": History of ";
	  errMsg += quantity + " can only be written, if it is activated ";
	  errMsg += "in section 'nodalResults', too.";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
	saveSolHist_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
  }
    
  // --- Acoustic Potential, 1. Deriv ---
  Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
  if (saveNodeHist.GetSize() > 0) {
	if (saveDeriv1_ == FALSE) {
	  std::string errMsg = pdename_;
	  errMsg += ": History of ";
	  errMsg += quantity + " can only be written, if it is activated ";
	  errMsg += "in section 'nodalResults', too.";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
	saveDeriv1Hist_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotentialD1 for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
  }
  // --- Acoustic Potential, 1. Deriv ---
  Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
  valVec  = "", "", quantity;
  params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
  if (saveNodeHist.GetSize() > 0) {
	if (saveDeriv2_ == FALSE) {
	  std::string errMsg = pdename_;
	  errMsg += ": History of ";
	  errMsg += quantity + " can only be written, if it is activated ";
	  errMsg += "in section 'nodalResults', too.";
	  Error(errMsg.c_str(), __FILE__, __LINE__);
	}
	saveDeriv1Hist_ = TRUE;
	Info->PrintF( pdename_, "Saving acouPotetentialD2 for Nodes:\n" );
	for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	  Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
	}
  }
    
}

} // end of namespace
