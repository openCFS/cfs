#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "Forms/forms_header.hh"
#include "Forms/curlfieldop.hh"
#include "Estimator/spaceerror.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "trapezoidal.hh"
#include "Utils/Coil.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Utils/SmoothSpline.hh"
#include "PDE/scalarnodeEQN.hh"
#include "magneticPDE.hh"

namespace CoupledField {

  // **************
  //  Constructor
  // **************
  MagPDE::MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		 FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) {

    ENTER_FCN( "MagPDE::MagPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================
    dofspernode_ = 1;
    solTypes_ = MAG_POTENTIAL;
    solDofs_ = 1;
    pdename_          = "magnetic";
    pdematerialclass_ = "magnetic";

    // ---------------------------------------------------------
    //   Get special geometry
    // ---------------------------------------------------------
    isaxi_ = params->HasValue( "type", "axi", "geometry" )
      == TRUE ? TRUE : FALSE;

    // ---------------------------
    //   Set coupling parameters
    // ---------------------------
    deltCoords_.Resize( dim_, numPDENodes_ );
    
  }


  // *************
  //  Destructor
  // *************
  MagPDE::~MagPDE() {
  ENTER_FCN( "MagPDE::~MagPDE" );
  for ( UInt k = 0; k < coilDef_.GetSize(); k++ ) {
    delete coilDef_[k];
  }
}


  // *********************************
  //  Read special boundary conitions
  // *********************************

  void MagPDE::ReadSpecialBCs()
  {
    ENTER_FCN( "MagPDE::ReadSpecialBCs" );

    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
    
    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();
  }
  

  // ****************************
  //  Initialize Nonlinearities
  // ****************************
  void MagPDE::InitNonLin()
  {
    ENTER_FCN( "MagPDE::InitNonLin" );

    // ------------------------------------
    //   Get information on non-linearity
    // ------------------------------------
    nonLin_ = FALSE;
    params->GetList( "nonLinear", nonLinType_, pdename_, "region" );
    
    for ( Integer k = 0; k < nonLinType_.GetSize(); k++ ) {
      if ( nonLinType_[k] != "no" ) {
	nonLin_ = TRUE;
	break;
      }
    }

    if( nonLin_ ) {

      // solution method
      params->Get( "method", nonLinMethod_, pdename_, "nonLinear" );

      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );

      // incremental stopping criterion
      params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      
      // residual stopping criterion
      params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );

      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }
  }


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagPDE::DefineIntegrators(const Integer level) {

    ENTER_FCN( "MagPDE::DefineIntegerators" );

    // Loop over all regions this PDE lives on
    for ( Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // Get reluctivity for this domain and perform consistency check
      Double reluctivity;
      materialData_[actSD].GetPermeability(2,2,reluctivity);

      if ( reluctivity == 0 ) {
	Info->Error( "Region '" + subdoms_[actSD] + "' has zero" +
		     " permeability!", __FILE__, __LINE__ );
      }
      
      reluctivity = 1/reluctivity;
 
      if ( nonLinType_[actSD] != "no" ) {

	//read in the BH-curve data and compute the approximation
	std::string nlfnc = materialData_[actSD].GetBHCurveFileName();
	ApproxData *nlinFnc = new SmoothSpline(nlfnc);
	nlinFnc->CalcBestParameter();
	nlinFnc->CalcApproximation();

	BaseForm *curlcurl2D = new nLinCurlCurlNode2DInt(nlinFnc,reluctivity,
							 isaxi_);
	curlcurl2D->SetNonLinMethod(nonLinMethod_);      
	assemble_->AddIntegrator(curlcurl2D, subdoms_[actSD], STIFFNESS, TRUE);

	// nonlinear RHS linearform!!
	BaseForm * rhsSource = new nLinMagNode2D_linFormInt(nlinFnc,
							    reluctivity,
							    isaxi_);
	assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], TRUE);
      }
      else {
	BaseForm *curlcurl2D = new CurlCurlNode2DInt(reluctivity, isaxi_);
	assemble_->AddIntegrator(curlcurl2D, subdoms_[actSD], STIFFNESS,FALSE);
	if (nonLin_==TRUE) {
	  // for nonlinear RHS linearform we need linear and nonlinear
	  // subdomains
	  BaseForm * rhsSource = new nLinMagNode2D_linFormInt(reluctivity,
							      isaxi_);
	  assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], TRUE);
	}
      }

      // Mass matrix
      Double conductivity;
      materialData_[actSD].GetConductivity(2,2,conductivity);
      BaseForm *bilinear_mass = new MassInt(conductivity, dofspernode_,isaxi_);
      assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, FALSE );

      // If this subdomain is a coil we have to do special things
      for ( Integer coil = 0; coil < coilDef_.GetSize(); coil++ ) {
	if ( subdoms_[actSD] == coilName_[coil] ) {
	  Double factor = coilDef_[coil]->value_ /
	    coilDef_[coil]->windingCrossSection_;
	  BaseForm *coilSource = new VolumeSrcInt( factor, isaxi_ );	  
	  assemble_->AddRhsSrcIntegrator( coilSource,
					  subdoms_[actSD],
					  coilDef_[coil]->dynamicsFile_,
					  nonLin_ );
	}
      }

      // check, if this subdomain is a permanent magnet
      for ( Integer perm = 0; perm < magnetsDomain_.GetSize(); perm++ ) {
	if ( subdoms_[actSD] == magnetsDomain_[perm] ) {

	  if (dim_ ==3)
	    Error("Permanent magnets for 3D computation not implementes");

	  //change direction of magnetization, so that we can use the
	  //standard GlobalDerivatives and obtain: (curl w) . M
	  Vector<Double> magnetization(dim_);
	  if (isaxi_) {
	    magnetization[0] = magnetsOriY_[perm];
	    magnetization[1] = -magnetsOriX_[perm];
	  }
	  else {
	    magnetization[0] = -magnetsOriY_[perm];
	    magnetization[1] = magnetsOriX_[perm];
	  }
	  
	  std::string fncname = "none";
	  Boolean nlin = FALSE;
	  BaseForm *permSource = new MagPerm2DInt(magnetization, 
						  reluctivity, isaxi_ );	  
	  assemble_->AddRhsSrcIntegrator( permSource,
					  subdoms_[actSD],
					  fncname,
					  nlin );
	}
      }

    }
  }


  // ======================================================
  // SOLVING SECTION
  // ======================================================

  void MagPDE :: InitTimeStepping() {
    ENTER_FCN( "MagPDE::InitTimeStepping" );
    
    TS_alg_ = new Trapezoidal(pdename_, algsys_, eqnData_);
  }


  void MagPDE:: PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) {
    ENTER_FCN( "MagPDE::PreStepStatic" );
    if (pdeIsCoupled_) 
      algsys_->InitSol();

  }


  void MagPDE::PostStepStatic(const Integer level) {
    ENTER_FCN( "MagPDE::PostStepStatic" );
    if (pdeIsCoupled_) iterCoupledCounter_++;
  }


  void MagPDE::StepStaticNonLin(const Integer kstep, const Double aTime,
				const Integer level, const Boolean reset)
  {
    ENTER_FCN( "MagPDE::SolveStepStaticNonLin" );

    const Integer job = 1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;
  
    Vector<Double> solInc(eqnData_->GetNumEQNs());
    Vector<Double> actSol(eqnData_->GetNumEQNs());

    sol_->GetAlgSysVector(actSol);

    SetBCs(level, 0);

    // setup right hand side ==============================================
    Double RhsLinL2Norm = SetLinRHS(level); 

    //add nonlinear part to RHS 
    assemble_->AssembleNLRHS(level);

    do {
      iterationCounter++;

      std::cout << std::endl << "Nonlinear Magnetics: Perform internal loop "
		<< "nr. " << iterationCounter << std::endl;      

#ifdef DEBUG
      *debug << std::endl << "=============================================="
	     << "========\n"
	     <<	"Nonlinear Mechanics: Perform internal loop nr. "
	     << iterationCounter << std::endl;      
#endif

      // setup and solve new system (rhs is already set) =====================
      assemble_->InitNonLinMatrices();
      assemble_->AssembleMatrices(level);
      
#ifdef USE_OLAS
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond(job);
#else
      algsys_->CalcPrecond(job);
#endif

      algsys_->Solve();

      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solInc, algsys_->GetSolutionVal() );
      
      actSol += solInc;
      sol_->SetAlgSysVector(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))========
#ifndef USE_OLAS    
      algsys_->InitRHS(RhsLinVal_.GetPointer());
#endif

      assemble_->AssembleNLRHS(level);  // nonlinear part of RHS


      // calculation of residual error (takes care for Dirichlet BCs========
      Vector<Double> actRHS;
      StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
	  
      Double residualL2Norm;
      residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )


      // calculation of residual error =======================================
      Double residualErr;
      if (RhsLinL2Norm > 1.0)
	residualErr = residualL2Norm / RhsLinL2Norm;
      else
	residualErr = residualL2Norm;

      // calculate incremental error ========================================
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm  = actSol.NormL2();
      
      Double incrementalErr;      
      if (actSolL2Norm > 1.0)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	incrementalErr = solIncrL2Norm;

      // =====================================================================
      // output of norms and data
      // =====================================================================
      Double etaLineSearch = 1;
      if ( nonLinLogging_ == TRUE ) {
	Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
			      incrementalErr, etaLineSearch);
      }
      

      // boolean variable, holds condition if another
      // iteration step is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
      
      if (!(performOneMoreStep && iterationCounter < nonLinMaxIter_))
	mycout << "incrementalErr " << incrementalErr << myendl 
	       << "incStopCrit_ " << incStopCrit_ << myendl
	       << "residualErr " << residualErr  << myendl
	       << "residualStopCrit_ " << residualStopCrit_ << myendl;

    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  }


  // **************
  //   LineSearch
  // **************
  Double MagPDE::LineSearch( Vector<Double>& solInc, Vector<Double>& actSol, 
			     Double& etaLineSearch, Integer level,
			     Boolean trans ) {
    ENTER_FCN( "MagPDE::LineSearch" );

    //  Vector<Double> solOld(actSol);
    const Integer nrEtas = 4;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;
  
    actSol += solInc;
  
  }
  

  void MagPDE::StepTransNonLin(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) {

    ENTER_FCN( "MagPDE::StepTransNonLin" );

    lasttimecalc_ = asteptime;
    laststepcalc_ = kstep;


    const Integer job = 1;
  
    static Integer timeStepCounter=1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;

    Vector<Double> actSol(numPDENodes_);
    Vector<Double> solInc(numPDENodes_);
  
    // Cast BaseStoreSol into StoreSol<Double>,
    // since this function is only called
    // in the transient case
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //set actual solution  
    actSol = solhelp->GetAlgSysVector();

    //compute predictors
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    //now set up RHS: all linear source terms
    Double RhsLinL2Norm = SetLinRHS(level); 

    // inner forces due to nonlin formulation
    assemble_->AssembleNLRHS(level, lasttimecalc_);  

    //Update RHS (mass matrix on right hand side)
    TS_alg_->UpdateRHS(solhelp->GetAlgSysVector());
  
    timeStepCounter++;
    do {
      iterationCounter++;
      std::cout << std::endl
		<< "Nonlinear Magnetics: Perform internal loop nr. " 
		<< iterationCounter << std::endl;

#ifdef DEBUG
      *debug << std::endl
	     << "====================================================== "
	     << std::endl
	     <<	"Nonlinear Magnetics: Perform internal loop nr. "
	     << iterationCounter << std::endl;      
#endif

      // setup and solve new system (rhs is already set) ====================
      assemble_->InitNonLinMatrices();
      assemble_->AssembleMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      SetBCs(level, lasttimecalc_);

#ifdef USE_OLAS
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond(job);
#else
      algsys_->CalcPrecond(job);
#endif

      algsys_->Solve();

      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solInc, algsys_->GetSolutionVal() );

      Double residualL2Norm;
      Double etaLineSearch = 0;
      
      if ( lineSearch_ != "no" ) {
	actSol += solInc;
      }
      else {

	// TRUE is for transient simulation
	residualL2Norm = LineSearch(solInc, actSol, etaLineSearch, level,TRUE);
      }

      //store A_/n+1) in the solution-object sol_
      sol_->SetAlgSysVector(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))=======
#ifndef USE_OLAS    
      algsys_->InitRHS(RhsLinVal_.GetPointer());
#endif

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(actSol);

      // inner forces due to nonlin formulation
      assemble_->AssembleNLRHS(level, lasttimecalc_);
 

      // ====================================================================
      // calculation of error norms
      // ====================================================================

      if ( lineSearch_ != "no" ) {

	Vector<Double> actRHS;
	StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
	  
	// ------------------------------------------------------------------
	// calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
	// ------------------------------------------------------------------
	residualL2Norm = RhsL2Norm(actRHS);
      }
      
      Double residualErr;
      if ( RhsLinL2Norm > 1)
	residualErr    = residualL2Norm /  RhsLinL2Norm;
      else
	residualErr    = residualL2Norm;

      // --------------------------------------------------------------------
      // calculate incremental error
      // --------------------------------------------------------------------
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm = actSol.NormL2();
      Double incrementalErr;
      
      if (actSolL2Norm > 1)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	incrementalErr = solIncrL2Norm;

      // --------------------------------------------------------------------
      // output of norms and data
      // --------------------------------------------------------------------
      if ( nonLinLogging_ == TRUE ) {
	Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
			      incrementalErr, etaLineSearch);
      }

      // boolean variable, holds condition if another iteration step
      // is necessary
      performOneMoreStep = 
	(incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }


  // sets excitation coil and returns L2Norm of them
  Double MagPDE::SetLinRHS(const Integer level)
  {
    ENTER_FCN( "MagPDE::SetCoilExcitations" );

    Double RhsLinL2Norm;  


   // to incorporate loads
    assemble_->AssembleSrcRHS(level, lasttimecalc_); 

    // Stores rhs vector into extForces and returns that L2-norm
    StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

    RhsLinL2Norm = RhsLinVal_.NormL2();
 
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      Warning("Zero external force vector!! ", __FILE__,__LINE__);
    }

    return RhsLinL2Norm;
  }

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty
  // formulation by setting them 0
  Double MagPDE::RhsL2Norm(Vector<Double>& actRHS) {
    ENTER_FCN( "MagPDE::RhsL2Norm" );

    Integer node;
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


  //   // calculates L2-norm of RHS regarding dirichlet entries due
  //   // to penalty formulation by setting them to zero


  //   Double MagPDE::RhsL2Norm(Vector<Double>& actRHS) {
  //     ENTER_FCN( "MagPDE::RhsL2Norm" );

  //     Integer node;
  //     std::list<Integer> nodes;

  //     // Eliminate dirichlet node from RHS (due to penalty formulation)
  //     for (Integer i=0; i< bcs_hd_.GetSize(); i++) {
  //       nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
  //       for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end();
  // 	   p++) {
  // 	node=*p;
  // 	actRHS[mesh2PDENode_[node-1]-1] = 0;
  //       }
  //     }
  //     return actRHS.NormL2();
  //   }

  // stores an algsys_ vector into a StdVector and returns that L2-norm
  void MagPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {
    ENTER_FCN( "MagPDE::StoreAlgsysToVec" );
    //const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    vec.Resize(numElems);
    for (Integer i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void MagPDE::WriteResultsInFile(Integer stepOffset,
				  Double timeOffset) {

    ENTER_FCN( "MagPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = lasttimecalc_ + timeOffset;
    Integer actStep = laststepcalc_ + stepOffset;
    
    if (analysistype_ == STATIC ||
	analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

      if (saveSol_)
	  outFile_->WriteNodeSolutionTransient(*solTransient, 
					       actStep, actTime);
      if (saveSolHist_)
	outFile_->WriteNodeHistoryTransient(*solTransient, 
					    actStep, actTime);
      
      if (calcBfield_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(B_, actStep, actTime);
      }
      
      if (calcEddy_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(Jeddy_, actStep, actTime);
      }
    } else  
      if (analysistype_ == HARMONIC) {
	
	solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
	
	if (saveSol_)
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actFreqStep_, 
					      actFrequency_, complexFormat_);

	if (saveSolHist_)
	  outFile_->WriteNodeHistoryHarmonic(*solHarmonic, actFreqStep_, 
					      actFrequency_, complexFormat_);
      }
      else
	Error("MagPDE: Only static, transient and harmonic results can be written",
	      __FILE__, __LINE__);
  }


  // ***************
  //   PostProcess
  // ***************
  void MagPDE::PostProcess(const Integer level) {

    ENTER_FCN( "MagPDE::PostProcess" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  
    if (calcEnergy_.GetSize() !=0 )
      CalcEnergy();

    if (calcBfield_.GetSize() !=0 ) {

      CurlNodeOp * FieldOp = new CurlNodeOp(ptgrid_, this, eqnData_,
					    *solhelp, level);
      FieldOp->Set2DType(isaxi_);
 
      // ------ Calculation of the electric field ------

      Vector<Double> LCoord;
      LCoord.Resize(dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      StdVector<Elem*> elemssd;
      Integer counterElems=0;
      Vector<Double> TempE;
      Integer pdeElem;
      
      // Resize solution arrays
      B_.SetNumSolutions(1);
      B_.SetSolutionType(MAG_FLUX_DENSITY);
      B_.SetNumElems(numElems_);
      B_.SetNumDofs(dim_);
      B_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      B_.Init(0);
      
      // loop over all subdomains
      for (Integer isd=0; isd<calcBfield_.GetSize(); isd++) {

	// get vector of Elem of subdomain with color: subdoms[isd]
	ptgrid_->GetElemSD(elemssd,calcBfield_[isd],level);
	  
	// loop over elements of subdomain
	for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++) {
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	  FieldOp->CalcElemCurlNode( TempE, elemssd[iel], LCoord); 
	  // B_.SetNodalResult(mesh2PDEElem_[elemssd[iel]->elemNum - 1]-1,
	  // TempE);
	  B_.SetElemResult(pdeElem-1, TempE);
	}
      }
      delete FieldOp;
    }

    if (calcEddy_.GetSize() !=0 ) {

      Vector<Double> LCoord;
      LCoord.Resize(dim_);
      LCoord[0] = 0;
      LCoord[1] = 0;
      
      StdVector<Elem*> elemssd;
      Vector<Double> ShpFnc, tmp;
      Vector<Double> magVecDeriv1Elem;
      StdVector<Integer> connect, connect_PDE;
      Double conductivity = 0.0;

      Integer counterElems=0;
      Integer pdeElem;

      // dimension hard coded for .unv file!
      Vector<Double> JeddyElem(3);
      JeddyElem.Init();
      
      // Resize solution arrays
      Jeddy_.SetNumSolutions(1);
      Jeddy_.SetSolutionType(MAG_EDDY_CURRENT);
      Jeddy_.SetNumElems(numElems_);

      // dimension hard coded for .unv file!
      Jeddy_.SetNumDofs(3);  
      Jeddy_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      Jeddy_.Init(0);

      // loop over all subdomains
      for (Integer actSD=0; actSD<calcEddy_.GetSize(); actSD++) {

	// get vector of Elem of subdomain with color: subdoms[isd]
	ptgrid_->GetElemSD(elemssd,calcEddy_[actSD],level);
	  
	// Get the right material parameter for actual subdomain
	for (Integer iSD=0; iSD<subdoms_.GetSize(); iSD++)
	  if (subdoms_[iSD] == calcEddy_[actSD])
            materialData_[iSD].GetConductivity(2,2,conductivity);

	// loop over elements of subdomain
	for ( Integer actEl=0; actEl< elemssd.GetSize();
	      actEl++,counterElems++ ) {
	  BaseFE * ptEl = elemssd[actEl]->ptElem;
	  ptEl->GetShFnc(ShpFnc,LCoord);
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[actEl]->elemNum);

	  connect = elemssd[actEl]->connect;
	  
	  GetDerivSolVecOfElement(magVecDeriv1Elem,connect);
	  JeddyElem[0] = magVecDeriv1Elem * ShpFnc;
	  JeddyElem[0] *= -conductivity;
	  // Jeddy_.SetNodalResult(mesh2PDEElem_[elemssd[actEl]->elemNum-1]-1,
	  // JeddyElem);
	  Jeddy_.SetElemResult(pdeElem-1, JeddyElem);
	}
      }
    }

    if (UIfilename_.size() > 0 && analysistype_ == TRANSIENT) {
      Vector<Double> uiSD;
      ComputeUI(uiSD);
      WriteUI2File(uiSD);
    }
  }


  // **************
  //   CalcEnergy
  // **************
  void MagPDE::CalcEnergy() {

    ENTER_FCN( "MagPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    StdVector<Integer> connecth, Eqns;  
    Vector<double> help;

    Integer i, j;
    Vector<Double> energy(calcEnergy_.GetSize());

    for (i=0; i<calcEnergy_.GetSize(); i++) {

      //reads eps33 (matrix notation starts with 0)
      Double eps33 = materialData_[i].GetPermittivity(2,2);

      StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,calcEnergy_[i],actlevel_);

      energy[i] = 0;
      for (j=0; j < elemssd.GetSize(); j++) {

	ptElem=elemssd[j]->ptElem;
	BaseForm * bilinear_stiff = new LaplaceInt(ptElem, eps33, isaxi_);

	connecth=elemssd[j]->connect;
	GetElemCoords(connecth, ptCoord, actlevel_);
	bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

	// 	  EqnData_->Mesh2Eqn(Eqns,connecth);
	// 	  (*debug) << "Nodes:" << connecth << std::endl;
	// 	  (*debug) << "Eqns :" << Eqns << std::endl;


	Vector<Double> magvecpot;
	
	sol_->GetElemSolution(magvecpot,connecth);
	
	help = elemmat * magvecpot;
	energy[i] += help * magvecpot;
	    
	delete bilinear_stiff;	  
      }  
    }

    std::string resulttype = "Electric Energy";
    Info->WriteResult(pdename_,  resulttype, calcEnergy_ , energy);
  }


  // *************
  //   ComputeUI
  // *************
  void MagPDE::ComputeUI(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::ComputeUI" );

    uiSD.Resize(coilDef_.GetSize());
  
    // loop over all subdomains
    for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {

      for (Integer dom=0; dom<coilDef_.GetSize(); dom++) {
	if (subdoms_[actSD] == coilName_[dom]) {
	   
	  StdVector<Elem*> elemssd;		
	  // get vector of Elem of subdomain with color: subdoms[isd]
	  ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
	    

	  // loop over elements of subdomain	    
	  for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++) {
	    BaseFE * ptEl = elemssd[actEl]->ptElem;
		
	    const Integer nrIntPts= ptEl->GetNumIntPoints();
	    const Vector<Double> & intWeights = ptEl->GetIntWeights();  
	    Double jacDet;
		
	    StdVector<Integer> connect;
	    connect = elemssd[actEl]->connect;

	    Matrix<Double> ptCoord;
	    GetElemCoords(connect, ptCoord, actlevel_);

	    Vector<Double> magVecDeriv1Elem;
	    GetDerivSolVecOfElement(magVecDeriv1Elem,connect);
	    Double uiElem=0;
		
	    for (Integer actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {
	      Vector<Double> shapeFnc;
	      jacDet = ptEl->CalcJacobianDetAtIp(actIntPt, ptCoord);	
	      ptEl -> GetShFncAtIp(shapeFnc, actIntPt);

	      uiElem += shapeFnc * magVecDeriv1Elem;
		    
	      if (isaxi_) {
		Vector<Double> coordAtIP = ptCoord * shapeFnc;
		uiElem += shapeFnc * magVecDeriv1Elem * 2 * PI * coordAtIP[0]
		  * intWeights[actIntPt-1];
	      }
	      else {
		uiElem += shapeFnc * magVecDeriv1Elem * intWeights[actIntPt-1];
	      }
	    }
	    
	    uiSD[dom] += uiElem;
	  }    
	}
      }
    }
  }


  // ***************************************************
  //   Store currents/voltages and inductivity in file
  // ***************************************************
  void MagPDE::WriteUI2File(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::WriteUI2File" );

    Vector<Integer> coilIDs;   // just positive ids
    coilIDs.Push_back(abs(coilDef_[0]->id_));

    Integer maxID = coilDef_[0]->id_;

    for (Integer dom=1; dom < coilDef_.GetSize(); dom++) {

      Boolean isInVec = FALSE;
      
      for (Integer dom2=0; dom2 < coilIDs.GetSize(); dom2++)	
	if (abs(coilDef_[dom]->id_) == coilIDs[dom2])
	  isInVec = TRUE;
      
      if (!isInVec)
	coilIDs.Push_back(abs(coilDef_[dom]->id_));
      
      if (abs(coilDef_[dom]->id_) > maxID)
	maxID = coilDef_[dom]->id_;
    }

    Vector<Double> uiID(maxID);
    uiID.Init();

    for (Integer dom=0; dom < coilDef_.GetSize(); dom++) {
      Integer actCoilID = coilDef_[dom]->id_;
      uiID[abs(actCoilID)-1] += uiSD[dom] * actCoilID/abs(actCoilID);
    }

    *UIfile_ << lasttimecalc_ << " \t";
    for (Integer actID=0; actID < coilIDs.GetSize(); actID++) {
      if ( coilDef_[coilIDs[actID]-1]->coilType_ == Coil::MEASUREMENT2D )
	*UIfile_ << uiID[coilIDs[actID]-1] *
	  coilDef_[coilIDs[actID]-1]->windingCrossSection_ << " \t";
    }
  
    *UIfile_ << myEndl;
  
  }


  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagPDE::ReadCoils() {

    ENTER_FCN( "MagEdgePDE::ReadCoils" );

    StdVector<std::string> coilNamesAux, helper;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Get list of coils for magnetic PDE
    params->GetCoilList( coilName_, pdename_, bcSequenceTag_);

//     keyVec = pdename_, "coils", "name";
//     attrVec= "", "tag";
//     valVec = "", bcSequenceTag_;
//     params->GetList(keyVec, attrVec, valVec, helper);

//     std::cerr << "helper = " << helper << std::endl;
//     for (Integer i=0; i<coilNamesAux.GetSize(); i++) {
//       for (Integer j=0; j<helper.GetSize(); j++)
// 	if (helper[j] == coilNamesAux[i])
// 	  coilName_.Push_back(coilNamesAux[i]);
//     }

//     std::cerr << "coilName_ = " << coilName_ << std::endl;
    // Read parameters for individual coils and log to info file
    UInt nrCoils = coilName_.GetSize();
    if ( nrCoils > 0 ) {
      Info->PrintF( pdename_, " Using the following coils:\n" );
      coilDef_.Reserve( nrCoils );
      for ( UInt k = 0; k < nrCoils; k++ ) {
	coilDef_.Push_back( new Coil( coilName_[k], pdename_ ) );
	Info->PrintCoil( (*coilDef_[k]), analysistype_ );
      }
    }
  }


  // ********************************************************
  //   Query parameter object for information about magnets
  // ********************************************************
  void MagPDE::ReadMagnets() {

    ENTER_FCN( "MagEdgePDE::ReadMagnets" );

    // get domain names of magnets
    //    params->GetList( "name", magnetsDomain_, pdename_, "magnets" );
  // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    keyVec = pdename_, "magnets", "magnet", "name";
    attrVec = "", "", "tag";
    valVec  = "", "", bcSequenceTag_;

    params->GetList(keyVec, attrVec, valVec, magnetsDomain_);

    if ( magnetsDomain_.GetSize() > 0 ) {

      Info->PrintF( pdename_,
		    " Found permanent magnets in the following regions:" );

      Double tmpDir;

      // Construct vectors for restricted search parameter
      StdVector<std::string> keyVec;
      StdVector<std::string> attrVec;
      StdVector<std::string> valVec;
      attrVec = "", "", "name";

      // for each magnet ...
      for ( UInt k = 0; k < magnetsDomain_.GetSize(); k++ ) {

	// ... read direction of magnetisation
	valVec = "", "", magnetsDomain_[k];

	keyVec  = pdename_, "magnets", "magnet", "orientX";
	params->Get( keyVec, attrVec, valVec, tmpDir);
	magnetsOriX_.Push_back( tmpDir);

	keyVec  = pdename_, "magnets", "magnet", "orientY";
	params->Get( keyVec, attrVec, valVec, tmpDir );
	magnetsOriY_.Push_back( tmpDir );

	keyVec  = pdename_, "magnets", "magnet", "orientZ";
	params->Get( keyVec, attrVec, valVec, tmpDir );
	magnetsOriZ_.Push_back( tmpDir );

	// ... report name to logfile
	Info->PrintF( pdename_, "%s", magnetsDomain_[k].c_str());
      }
    }
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MagPDE::ReadStoreResults() {

    ENTER_FCN( "MagPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string quantity;

    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    Enum2String(MAG_POTENTIAL, quantity);
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
   	  hasOutput_ = TRUE;
    }

    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elemResults;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";  
    
    // --- Magnetic Flux Density ---
    Enum2String(MAG_FLUX_DENSITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, calcBfield_ );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcBfield_.GetSize() == 1 && calcBfield_[0] == "all" ) {
      calcBfield_ = subdoms_;
    }
    
    // Log to info file
    if ( calcBfield_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
		    " Computing magFluxDensity for regions:");
      for ( Integer k = 0; k < calcBfield_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcBfield_[k].c_str() );
      }
    }
    
    // --- Magnetic Energy ---
    Enum2String(MAG_ENERGY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, calcEnergy_ );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == "all" ) {
      calcEnergy_ = subdoms_;
    }
    
    // Log to info file
    if ( calcEnergy_.GetSize() > 0 ) {
	  hasOutput_ = TRUE;	
      Info->PrintF( pdename_,
		    " Computing magEnergy for regions:");
      for ( Integer k = 0; k < calcEnergy_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcEnergy_[k].c_str() );
      }
    }

   // --- Magnetic Eddy Current ---
    Enum2String(MAG_EDDY_CURRENT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, calcEddy_ );
    
    // If the symbolic name is "all" compute magnetic Fieldfor all regions
    if ( calcEddy_.GetSize() == 1 && calcEddy_[0] == "all" ) {
      calcEddy_ = subdoms_;
    }
    
    // Log to info file
    if ( calcEddy_.GetSize() > 0 ) {
      hasOutput_ =TRUE;
      Info->PrintF( pdename_,
		    " Computing magEddyCurrent for regions:");
      for ( Integer k = 0; k < calcEddy_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcEddy_[k].c_str() );
      }
    }
   
    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist;
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";
    
    // --- Magnetic Potential ---
    Enum2String(MAG_POTENTIAL, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving magPotential for Nodes:" );
      for ( Integer k = 0; k < saveNodeHist.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", saveNodeHist[k].c_str() );
      }
    }

    // *****************************
    // Determine element history
    // *****************************
    StdVector<std::string> saveElemHist;
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "";
    valVec = "", "", "";
    params->GetList(keyVec, attrVec, valVec, saveElemHist);
    
    if (saveElemHist.GetSize() > 0) {
      std::string errMsg = pdename_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }
    
    
  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================


void MagPDE::InitCoupling(PDECoupling * Coupling)
{
  ENTER_FCN( "MagPDE::InitCoupling" );
  
  pdeIsCoupled_ = TRUE;
  ptCoupling_   = Coupling;


  // Enable update of geometry
  const Integer numCouplings = ptCoupling_->GetNumOutputCouplings();  

  StdVector<StdVector<Integer> > elemNodeToCouplingNode_tmp;
  elemNodeToCouplingNode_.Resize(numCouplings);


 for (Integer actCoupling=0; actCoupling<numCouplings; actCoupling++)
    {
      if (ptCoupling_->GetOutputQuantity(actCoupling) == MAG_FORCE_LORENTZ)
	{
	  // Intialize the memory of the coupling values
	  ptCoupling_->CreateCouplingVector(actCoupling,isComplex_);

	  //get the element-node to coupling node matching
	  StdVector<std::string> couplRegions;
	  ptCoupling_->GetOutputRegions(actCoupling, couplRegions);

	  //Get total number of coupling elements
	  Integer totalCouplingElems = ptgrid_->GetMaxnumElem(actlevel_, couplRegions);;

	  elemNodeToCouplingNode_tmp.Clear();
	  elemNodeToCouplingNode_tmp.Resize(totalCouplingElems);

	  Integer offset = 0;
	  for (Integer reg=0; reg<couplRegions.GetSize(); reg++) 
	    {
	      //find subdomain index
	      Integer SDidx=-1; for (Integer sd=0; sd<subdoms_.GetSize(); sd++) 
		{
		  if (couplRegions[reg] == subdoms_[sd]) 
		    {
		      SDidx = sd;
		      break;
		    }
		}
	      
	      if (SDidx==-1)
		Error("magnetics: Coupling Region is not within the subdomains of the PDE");
	    
	      StdVector<Elem*> elemssd;
	      ptgrid_->GetElemSD(elemssd, subdoms_[SDidx], actlevel_);

	      StdVector<Integer> * couplingnodes;
	      ptCoupling_->GetOutputNodes(actCoupling, couplingnodes);
	      if (couplingnodes == 0)
		Error("magnetics: Couplingnodes = 0!!!!");

	      for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++) 
		{
		  StdVector<Integer> connecth = elemssd[actEl]->connect;
		  elemNodeToCouplingNode_tmp[offset+actEl].Resize(connecth.GetSize());

		  for (Integer ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++)
		    for (Integer cnode=0; cnode<(*couplingnodes).GetSize(); cnode++)
		      if (connecth[ielemnode] == (*couplingnodes)[cnode] ) 
			{
			  elemNodeToCouplingNode_tmp[offset+actEl][ielemnode] = cnode;
			  break;
			} // end if
		}

	      //in the case, that we have more than one coupling region!
	      offset = elemssd.GetSize();
	    }

	  elemNodeToCouplingNode_[actCoupling]  = elemNodeToCouplingNode_tmp;
	}
    }
}

void MagPDE::CalcOutputCoupling()
{
  ENTER_FCN( "MagPDE::CalcOutputCoupling" );

  SolutionType quantity;
  StdVector<Integer> * couplingNodes     = NULL;
  CFSVector * values = NULL;
  Integer forcesCount = 0;

  // loop over all output coupling quantities
  for (Integer actCoupling=0; actCoupling<ptCoupling_->GetNumOutputCouplings(); actCoupling++)
    {
      quantity = ptCoupling_->GetOutputQuantity(actCoupling);
      ptCoupling_->GetOutputValues(actCoupling, values);

      Vector<Double> *  temp = dynamic_cast<Vector<Double> *>(values);
      
      switch(ptCoupling_->GetOutputType(actCoupling))
	{
	  
	case NODE:	  
	  ptCoupling_->GetOutputNodes(actCoupling, couplingNodes);
	  
	  if (quantity == MAG_FORCE_LORENTZ) 
	    {
	      CalcNodeForceLorentz(*temp, elemNodeToCouplingNode_[forcesCount],
				   actCoupling, couplingNodes->GetSize());
	      forcesCount++;;
	    }

	  break;
	  
	case ELEM:
	  Error("No Element coupling output", __FILE__,__LINE__);
	}
    }
}

void MagPDE::CalcNodeForceLorentz(Vector<Double> & force, 
				  StdVector<StdVector<Integer> > & elemNodeToCouplingNode,
				  Integer actCoupling, Integer numCouplingNodes)
{
  ENTER_FCN( "MagPDE::CalcNodeForceLorentz" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double> *>(sol_);

  //get the coupling regions
  StdVector<std::string> couplRegions;
  ptCoupling_->GetOutputRegions(actCoupling, couplRegions);

  
  MagLorentzForceOp * ForceOp = new MagLorentzForceOp(ptgrid_, this, eqnData_, *solhelp, actlevel_, isaxi_);
   
  force.Init(0.0);

  Vector<Double> Jeddy;

  Integer offset = 0;
  for (Integer reg=0; reg<couplRegions.GetSize(); reg++) 
    {
      //find subdomain index
      Integer SDidx=-1; for (Integer sd=0; sd<subdoms_.GetSize(); sd++) 
	{
	  if (couplRegions[reg] == subdoms_[sd]) 
	    {
	      SDidx = sd;
	      break;
	    }
	}

      Double conductivity;
      materialData_[SDidx].GetConductivity(2,2,conductivity);      
      
      StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd, subdoms_[SDidx], actlevel_);
  
      for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++) 
	{
	  StdVector<Integer> connecth = elemssd[actEl]->connect;
	  Matrix<Double> elemForce;
	  GetDerivSolVecOfElement(Jeddy,connecth);
	  Jeddy *= -conductivity;

	  ForceOp->CalcElemMagLorentzForce(elemForce, Jeddy, elemssd[actEl]);

	  // Add the element force to the according coupling node
	  for (Integer ielemnode=0; ielemnode<connecth.GetSize(); ielemnode++)
	    for( Integer idim=0; idim<dim_; idim++)
	      force[elemNodeToCouplingNode[actEl+offset][ielemnode]*dim_+idim] += 
		elemForce[ielemnode][idim];
	}

      //in the case, that we have more than one coupling region!
      offset = elemssd.GetSize();
	  
    }
}

Boolean MagPDE::HasOutput(SolutionType output)
{
  ENTER_FCN( "MagPDE::HasOutput" );

  if (output == MAG_FORCE_LORENTZ)
    return TRUE;

  return FALSE;
}



} // end of namespace

