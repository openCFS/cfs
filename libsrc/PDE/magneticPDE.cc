#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "Forms/forms_header.hh"
#include "Forms/elecfieldop.hh"
#include "Forms/elecforceop.hh"
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

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Old version without XMLPARAMS and with old way of treating coil parameters
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef XMLPARAMS

  MagPDE::MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		 FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) {

    ENTER_FCN( "MagPDE::MagPDE" );

    dofspernode_ = 1;
  
    pdename_          = "magnetic";
    pdematerialclass_ = "magnetic";
  
    conf->getsubdompde(subdoms_,pdename_);
    ReadBCs(pdename_);

    //check, if problem is axisymmetric
    isaxi_ = FALSE;
    std::string subtype;
    conf->ifget("subtype",subtype,pdename_);
    if (subtype == "axi")
      isaxi_ = TRUE;

    //check for coils
    ReadCoils();

    conf->ifget("calc_UI",UIfilename_,pdename_);
    if (UIfilename_.size() > 0)
      {    
	std::string ui_info = "File name for saving voltages/currents is: "
	  + UIfilename_ + "\n";
	Info->PrintF(pdename_,"%s",ui_info.c_str());
	UIfile_ = new std::ofstream(UIfilename_.c_str());

	if (!UIfile_) 
	  Error("Can't open UI-file");

      }
  

    //check for permanent magnets
    conf->ifgetliststr("list_magnets", magnetsDomain_, pdename_);


    //check for nonlinearity
    nonLin_ = FALSE;
    nonLin_ = conf->get_option( "nonlin",  pdename_ );

    if( nonLin_ == TRUE ) {

      // incremental stopping criterion
      incStopCrit_ = 1e-3;
      conf->ifget("incStopCrit", incStopCrit_, pdename_);

      // residual stopping criterion
      residualStopCrit_ = 1e-3;
      conf->ifget("residualStopCrit", residualStopCrit_, pdename_);

      // maximal number of NL-iterations
      nonLinMaxIter_ = 100;
      conf->ifget("nonlinMaxIter", nonLinMaxIter_, pdename_);
    }

    //check for postprocessing
    conf->ifgetliststr("calc_BField",calcBfield_,pdename_); 
    conf->ifgetliststr("calc_Energy",calcEnergy_,pdename_); 
    conf->ifgetliststr("calc_Eddy",calcEddy_,pdename_); 

    // initialize eqation data object
    eqnData_  = new ScalarNodeEQN(ptgrid_, ptBCs_, subdoms_, actlevel_, dofspernode_);    
    eqnData_->SetHomoDirichletBCs(bcs_hd_, homDirichDof_);
    eqnData_->CalcMapping();
    //eqnData_->Print(std::cerr);
    
    numPDENodes_ = eqnData_->GetNumLocalNodes();
    numElems_ = eqnData_->GetNumLocalElems();

    deltCoords_.Resize(dim_,numPDENodes_);

  
    // set analysis parameters
    assemble_->SetPtr2EQNData(eqnData_);
    assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_,
				surfdoms_);
    assemble_->SetGraphType(NODEGRAPH);

#ifdef USE_OLAS
    assemble_->SetMatrixEntryType(OLAS::DOUBLE);
    assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
#else
    assemble_->SetMatrixType(RSCALAR);
#endif  
  
    assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    

    // Initalize solution class
    sol_->SetNumSolutions(1);
    sol_->SetSolutionType(MAG_POTENTIAL);
    sol_->SetNumNodes(numPDENodes_);
    sol_->SetNumDofs(dofspernode_);
    sol_->SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    sol_->Init();
    
    ReadMaterialData();
   
    DefineIntegrators(actlevel_);  

    ReadSavings();
  }



  void MagPDE::DefineIntegrators(const Integer level)
  {
    ENTER_FCN( "MagPDE::DefineIntegerators" );

    Boolean nonLin = FALSE;
    Boolean iscoil;
    Integer idxcoil;

    for (int actSD = 0; actSD < subdoms_.GetSize(); actSD++)
      {
	Double reluctivity  = materialData_[actSD].GetPermeability();
	if ( reluctivity == 0)
	  Error("Permeability can not be zero!");

	reluctivity = 1/reluctivity;

	nonLin = FALSE;

	BaseForm * curlcurl2D;

	if (actSD == 1 && nonLin_==TRUE)
	  {
	    //just for testing

	    //read in the BH-curve data and compute the approximation
	    std::string nlfnc = "bh1.fnc";
	    ApproxData *nlinFnc = new SmoothSpline(nlfnc);
	    nlinFnc->CalcBestParameter();
	    nlinFnc->CalcApproximation();

	    curlcurl2D = new nLinCurlCurlNode2DInt(nlinFnc,reluctivity, isaxi_);
	    nonLin = TRUE;
	    assemble_->AddIntegrator(curlcurl2D, subdoms_[actSD], STIFFNESS,
				     nonLin);

	    // do RHS!!
	    BaseForm * rhsSource = new nLinMagNode2D_linFormInt(nlinFnc, reluctivity, isaxi_);
	    assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
	  }
	else
	  {
	    curlcurl2D = new CurlCurlNode2DInt(reluctivity, isaxi_);
	    assemble_->AddIntegrator(curlcurl2D, subdoms_[actSD], STIFFNESS,
				     nonLin);
	    if (nonLin_==TRUE) {
	      // do RHS!!
	      BaseForm * rhsSource = new nLinMagNode2D_linFormInt(reluctivity, isaxi_);
	      assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
	    }
	  }

	Double conductivity = materialData_[actSD].GetConductivity();      
	BaseForm * bilinear_mass  = new MassInt(conductivity, dofspernode_,
						isaxi_);
	assemble_->AddIntegrator(bilinear_mass, subdoms_[actSD], MASS, nonLin);


	//=========================== RHS =====================================
	iscoil = FALSE;
	for (Integer dom=0; dom<coilDef_.GetSize(); dom++)
	  if (subdoms_[actSD] == coilName_[dom]) 
	    {
	      iscoil = TRUE;
	      idxcoil = dom;
	    }

	if (iscoil)
	  {
	    Double factor = coilDef_[idxcoil].current /
	      coilDef_[idxcoil].coilArea;
	    BaseForm * coilSource = new VolumeSrcInt(factor,isaxi_);
	    assemble_->AddRhsSrcIntegrator(coilSource, subdoms_[actSD],
					   coilDef_[idxcoil].timefnc, nonLin_);
	  }
	 
      }
  }


  // ======================================================
  // SOLVING SECTION
  // ======================================================
  Double MagPDE::LineSearch(Vector<Double>& solInc, Vector<Double>& actSol, 
			    Double& etaLineSearch, Integer level, Boolean trans)
  {
    ENTER_FCN( "MagPDE::LineSearch" );

    //  Vector<Double> solOld(actSol);
    const Integer nrEtas = 4;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;
  
    actSol += solInc;
  
  }

  void MagPDE::StepTransNonLin(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset)
  {
    ENTER_FCN( "MagPDE::StepTransNonLin" );

    lasttimecalc_ = asteptime;
    laststepcalc_ = kstep;


    const Integer job = 1;
    const Integer update = 0;  
  
    static Integer timeStepCounter=1;
    Double * ptsol;
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
    do
      {
	iterationCounter++;
	std::cout << std::endl << "Nonlinear Magnetics: Perform internal loop nr. " 
		  << iterationCounter << std::endl;

#ifdef DEBUG
	*debug << std::endl << "====================================================== " << std::endl
	       <<	"Nonlinear Magnetics: Perform internal loop nr. " << iterationCounter << std::endl;      
#endif



	// setup and solve new system (rhs is already set) =====================
	assemble_->InitNonLinMatrices();
	assemble_->AssembleMatrices(level);
	algsys_->ConstructEffectiveMatrix(matrix_factor_);

	SetBCs(level, update, lasttimecalc_);

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
      
	if (!lineSearch_)
	  actSol += solInc;
	else
	  // TRUE is for transient simulation
	  residualL2Norm = LineSearch(solInc, actSol, etaLineSearch, level, TRUE);


	//store A_/n+1) in the solution-object sol_
	sol_->SetAlgSysVector(actSol);


	// recalculate RHS with new values to get new residual (f^(k+1))========
#ifndef USE_OLAS    
	algsys_->InitRHS(RhsLinVal_.GetPointer());
#endif

	//Update RHS (mass matrix on right hand side)
	TS_alg_->UpdateRHS(actSol);

	assemble_->AssembleNLRHS(level, lasttimecalc_);  // inner forces due to nonlin formulation
 

	// =====================================================================
	// calculation of error norms
	// =====================================================================

	if (!lineSearch_)
	  {
	    Vector<Double> actRHS;
	    StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
	  
	    // calculation of residual error =======================================
	      residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
	    }
      
	Double residualErr;
	if ( RhsLinL2Norm > 1)
	  residualErr    = residualL2Norm /  RhsLinL2Norm;
	else
	  residualErr    = residualL2Norm;

	// calculate incremental error ========================================
	Double solIncrL2Norm = solInc.NormL2();
	Double actSolL2Norm = actSol.NormL2();
	Double incrementalErr;

	if (actSolL2Norm > 1)
	  incrementalErr = solIncrL2Norm / actSolL2Norm;
	else
	  incrementalErr = solIncrL2Norm;

	// =====================================================================
	// output of norms and data
	// =====================================================================


	if ( nonLinLogging_ == TRUE ) {
	  Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
				incrementalErr, etaLineSearch);
	}

	// boolean variable, holds condition if another iteration step is necessary
	performOneMoreStep = 
	  (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < maxnumiter_);  

  
    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }

  void MagPDE :: InitTimeStepping(const Double dt) {
    ENTER_FCN( "MagPDE::InitTimeStepping" );
    TS_alg_ = new Trapezoidal(pdename_, algsys_, eqnData_);
    TS_alg_->Init(matrix_factor_, dt);
  }


  void MagPDE:: PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) {
    ENTER_FCN( "MagPDE::PreStepStatic" );
    if (pdeIsCoupled_) 
      algsys_->InitSol();
  }


  void MagPDE::StepStaticNonLin(const Integer kstep, const Double aTime,
				const Integer level, const Boolean reset)
  {
    ENTER_FCN( "MagPDE::SolveStepStaticNonLin" );

    const Integer job = 1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;
  
    Vector<Double> solInc(numPDENodes_);
    Vector<Double> actSol(numPDENodes_);

    sol_->GetAlgSysVector(actSol);

    SetBCs(level, updateBCs_, 0);

    // setup right hand side ==============================================      
    Double RhsLinL2Norm = SetLinRHS(level); 

    //add nonlinear part to RHS 
    assemble_->AssembleNLRHS(level);

    do
      {
	iterationCounter++;

	std::cout << std::endl << "Nonlinear Magnetics: Perform internal loop nr. " 
		  << iterationCounter << std::endl;      

#ifdef DEBUG
	*debug << std::endl << "====================================================== " << std::endl
	       <<	"Nonlinear Mechanics: Perform internal loop nr. " << iterationCounter << std::endl;      
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
	Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);


	// boolean variable, holds condition if another iteration step is necessary
	performOneMoreStep = 
	  (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
      
	if (!(performOneMoreStep && iterationCounter < maxnumiter_))
	  mycout << "incrementalErr " << incrementalErr << myendl 
		 << "incStopCrit_ " << incStopCrit_ << myendl
		 << "residualErr " << residualErr  << myendl
		 << "residualStopCrit_ " << residualStopCrit_ << myendl;

      }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  }


  // sets excitation coil and returns L2Norm of them
  Double MagPDE::SetLinRHS(const Integer level)
  {
    ENTER_FCN( "MagPDE::SetCoilExcitations" );

    Double RhsLinL2Norm;  

    // to incorporate loads
    assemble_->AssembleSrcRHS(level, lasttimecalc_); 

    // stores rhs vector into extForces and returns that L2-norm
    StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

    RhsLinL2Norm = RhsLinVal_.NormL2();
 
    //  if extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm)
      Warning("Zero external force vector!! ", __FILE__,__LINE__);
  
    return RhsLinL2Norm;
  }

  Double MagPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "MagPDE::RhsL2Norm" );

    Integer node, eqnNr, eqnDof;

  
    std::list<Integer> nodes;
  
    // Eliminate dirichlet node from RHS (due to penalty formulation)
    for (Integer i=0; i< bcs_hd_.GetSize(); i++)
      {
	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	  {
	    node=*p;
	    eqnData_->Node2EQN(node,1,eqnNr,eqnDof);
	    if (eqnNr != 0){
	      actRHS[(eqnNr-1)] = 0;
	    }
	  }
      }
    return actRHS.NormL2();
  }

  //   // calculates L2-norm of RHS regarding dirichlet entries due to penalty formulation by setting them 0
  //   Double MagPDE::RhsL2Norm(Vector<Double>& actRHS)
  //   {
  //     ENTER_FCN( "MagPDE::RhsL2Norm" );

  //     Integer node;
  
  //     std::list<Integer> nodes;

  //     // Eliminate dirichlet node from RHS (due to penalty formulation)
  //     for (Integer i=0; i< bcs_hd_.GetSize(); i++)
  //       {
  // 	nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
  // 	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
  // 	  {
  // 	    node=*p;
  // 	    actRHS[mesh2PDENode_[node-1]-1] = 0;
  // 	  }
  //       }

  //     return actRHS.NormL2();
  //   }

  void MagPDE::PostStepStatic(const Integer level) {
    ENTER_FCN( "MagPDE::PostStepStatic" );
    if (pdeIsCoupled_) iterCoupledCounter_++;
  }



  // stores an algsys_ vector into a StdVector and returns that L2-norm
  void MagPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt)
  {
    ENTER_FCN( "MagPDE::StoreAlgsysToVec" );

    //    const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    vec.Resize(numElems);

    for (Integer i=0; i<numElems; i++)   
      vec[i] = pt[i];
  }


  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================

  void MagPDE::WriteResultsInFile() {

    ENTER_FCN( "MagPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    if (analysistype_ == STATIC ||
	analysistype_ == TRANSIENT) {

      // transient/static case
      if (savesol_)
	{
	  solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
	  outFile_->WriteNodeSolutionTransient(*solTransient, 
					       laststepcalc_, lasttimecalc_);
	}
      
      if (calcBfield_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(B_, laststepcalc_, lasttimecalc_);
      }
      
      if (calcEddy_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(Jeddy_, laststepcalc_, lasttimecalc_);
      }
    } else  if (analysistype_ == HARMONIC) {
      
      // harmonic case
      if (savesol_)
	{
	  solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actFreqStep_, 
					      actFrequency_, complexFormat_);
	}
    }
    else
      Error("MagPDE: Only static, transient and harmonic results can be written",
	    __FILE__, __LINE__);
  }

	
  

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
      B_.SetSolutionType(MAG_FIELD);
      B_.SetNumElems(numElems_);
      B_.SetNumDofs(dim_);
      B_.SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
      B_.Init(0);
      
      // loop over all subdomains
      for (Integer isd=0; isd<calcBfield_.GetSize(); isd++) {

	// get vector of Elem of subdomain with color: subdoms[isd]
	ptgrid_->GetElemSD(elemssd,calcBfield_[isd],level);
	  
	// loop over elements of subdomain
	for (Integer iel=0; iel< elemssd.GetSize(); iel++,counterElems++) 
	  {
	    pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
	    FieldOp->CalcElemCurlNode( TempE, elemssd[iel], LCoord); 
	    // B_.SetNodalResult(mesh2PDEElem_[elemssd[iel]->elemNum - 1]-1, TempE);
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
	    conductivity = materialData_[iSD].GetConductivity(); 	  

	// loop over elements of subdomain
	for (Integer actEl=0; actEl< elemssd.GetSize(); actEl++,counterElems++) {
	  BaseFE * ptEl = elemssd[actEl]->ptElem;
	  ptEl->GetShFnc(ShpFnc,LCoord);
	  pdeElem = eqnData_->Mesh2PDEElem(elemssd[actEl]->elemNum);

	  connect = elemssd[actEl]->connect;
	  
	  GetDerivSolVecOfElement(magVecDeriv1Elem,connect);
	  JeddyElem[0] = magVecDeriv1Elem * ShpFnc;
	  JeddyElem[0] *= -conductivity;
	  // 	  Jeddy_.SetNodalResult(mesh2PDEElem_[elemssd[actEl]->elemNum - 1]-1,
	  // 				JeddyElem);
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
	    const Integer nrNodes = ptEl->GetNumNodes();
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
	      ptEl->GetShFncAtIp(shapeFnc, actIntPt);

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


  void MagPDE::WriteUI2File(Vector<Double>& uiSD) {

    ENTER_FCN( "MagPDE::WriteUI2File" );

    Vector<Integer> coilIDs;   // just positive ids
    coilIDs.Push_back(abs(coilDef_[0].ID));

    Integer maxID = coilDef_[0].ID;

    for (Integer dom=1; dom < coilDef_.GetSize(); dom++) {

      Boolean isInVec = FALSE;
      
      for (Integer dom2=0; dom2 < coilIDs.GetSize(); dom2++)	
	if (abs(coilDef_[dom].ID) == coilIDs[dom2])
	  isInVec = TRUE;
      
      if (!isInVec)
	coilIDs.Push_back(abs(coilDef_[dom].ID));
      
      if (abs(coilDef_[dom].ID) > maxID)
	maxID = coilDef_[dom].ID;
    }

    Vector<Double> uiID(maxID);
    uiID.Init();

    for (Integer dom=0; dom < coilDef_.GetSize(); dom++) {
      Integer actCoilID = coilDef_[dom].ID;
      uiID[abs(actCoilID)-1] += uiSD[dom] * actCoilID/abs(actCoilID);
    }

    *UIfile_ << lasttimecalc_ << " \t";
    for (Integer actID=0; actID < coilIDs.GetSize(); actID++) {
      if ( coilDef_[coilIDs[actID]-1].type == MEASUREMENT)   
	*UIfile_ << uiID[coilIDs[actID]-1] *
	  coilDef_[coilIDs[actID]-1].coilArea << " \t";
    }
  
    *UIfile_ << myEndl;
  
  }


  // reads all information in the config file concerning coils 
  void MagPDE::ReadCoils() {

    ENTER_FCN( "MagEdgePDE::ReadCoils" );

    conf->ifgetliststr("list_coils", coilName_, pdename_);
    Integer nrCoils = coilName_.GetSize();
    if (nrCoils) {
      coilDef_.Resize(nrCoils);
      for (Integer i=0; i < nrCoils; i++) {
	conf->getCoilData(coilDef_[i],pdename_,coilName_[i]+":");
	Info->PrintCoil(coilName_[i], coilDef_[i], analysistype_);
      }
    }
  }

  // ======================================================
  // COUPLING SECTION
  // ======================================================

#else // for #ifndef XMLPARAMS


  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // New version with XMLPARAMS and different way of treating coil parameters
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  // **************
  //  Constructor
  // **************
  MagPDE::MagPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc,
		 FileType *aptFileType, WriteResults *aptOut)
    :BasePDE(aptgrid, aptbcs, aptFileType, aptOut, aptTimeFunc) {

    ENTER_FCN( "MagPDE::MagPDE" );

    // ---------------------------
    //   Set PDE characteristics
    // ---------------------------
    dofspernode_ = 1;
    pdename_          = "magnetic";
    pdematerialclass_ = "magnetic";

    // ---------------------------------------------------------
    //   Get regions, boundary conditions and special geometry
    // ---------------------------------------------------------
    params->GetList( "name", subdoms_, pdename_, "region" );
    ReadBCs(pdename_);
    // once isaxi_ is bool:
    //isaxi_ = params->HasValue( "type", "axi", "geometry" );
    // so long:
    isaxi_ = params->HasValue( "type", "axi", "geometry" )
      == TRUE ? TRUE : FALSE;

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


    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();

    // -----------------------------
    // Check for permanent magnets
    // -----------------------------
    ReadMagnets();

    // ----------------------------
    //   Initalize equation data class
    // ----------------------------
    eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				  dofspernode_ );
    eqnData_->SetHomoDirichletBCs( bcs_hd_, homDirichDof_ );
    eqnData_->CalcMapping();
    //eqnData_->Print(std::cerr);
    
    numPDENodes_ = eqnData_->GetNumLocalNodes();
    numElems_ = eqnData_->GetNumLocalElems();

    // ---------------------------
    //   Set coupling parameters
    // ---------------------------
    deltCoords_.Resize( dim_, numPDENodes_ );

    // ---------------------------
    //   Set analysis parameters
    // ---------------------------
    assemble_->SetPtr2EQNData(eqnData_); 
    assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_,
				surfdoms_);
    assemble_->SetGraphType(NODEGRAPH);
   
#ifdef USE_OLAS
    assemble_->SetMatrixEntryType(OLAS::DOUBLE);
    assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
#else
    assemble_->SetMatrixType(RSCALAR);
#endif  
  
    assemble_->SetNumDirichlet(GetNumRestraints(actlevel_));
    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    // ----------------------------
    //   Initalize solution class
    // ----------------------------
    sol_->SetNumSolutions(1);
    sol_->SetSolutionType(MAG_POTENTIAL);
    sol_->SetNumNodes(numPDENodes_);
    sol_->SetNumDofs(dofspernode_);
    sol_->SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    sol_->Init();
    ReadMaterialData();
   
    DefineIntegrators(actlevel_);


    // ----------------------------
    //   Get post-processing info
    // ----------------------------
    ReadStoreResults();

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


  // *****************************
  //   Definition of Integrators
  // *****************************
  void MagPDE::DefineIntegrators(const Integer level) {

    ENTER_FCN( "MagPDE::DefineIntegerators" );

    // Loop over all regions this PDE lives on
    for ( Integer actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // Get reluctivity for this domain and perform consistency check
      Double reluctivity = materialData_[actSD].GetPermeability();
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
      Double conductivity = materialData_[actSD].GetConductivity();      
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
    }
  }


  // ======================================================
  // SOLVING SECTION
  // ======================================================

  void MagPDE :: InitTimeStepping(const Double dt) {
    ENTER_FCN( "MagPDE::InitTimeStepping" );
    TS_alg_ = new Trapezoidal(pdename_, algsys_, eqnData_);
    TS_alg_->Init(matrix_factor_, dt);
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

    SetBCs(level, updateBCs_, 0);

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
      
      if (!(performOneMoreStep && iterationCounter < maxnumiter_))
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
    const Integer update = 0;  
  
    static Integer timeStepCounter=1;
    Double * ptsol;
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

      SetBCs(level, update, lasttimecalc_);

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
      
    } while(performOneMoreStep && iterationCounter < maxnumiter_);  

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

  void MagPDE::WriteResultsInFile() {

    ENTER_FCN( "MagPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    if (analysistype_ == STATIC ||
	analysistype_ == TRANSIENT) {

      // transient/static case
      if (savesol_)
	{
	  solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
	  outFile_->WriteNodeSolutionTransient(*solTransient, 
					       laststepcalc_, lasttimecalc_);
	}
      
      if (calcBfield_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(B_, laststepcalc_, lasttimecalc_);
      }
      
      if (calcEddy_.GetSize() !=0 ) {
	outFile_->WriteElemSolutionTransient(Jeddy_, laststepcalc_, lasttimecalc_);
      }
    } else  if (analysistype_ == HARMONIC) {
      
      // harmonic case
      if (savesol_)
	{
	  solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
	  outFile_->WriteNodeSolutionHarmonic(*solHarmonic, actFreqStep_, 
					      actFrequency_, complexFormat_);
	}
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
      B_.SetSolutionType(MAG_FIELD);
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
	    conductivity = materialData_[iSD].GetConductivity(); 	  

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
	    const Integer nrNodes = ptEl->GetNumNodes();
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

    // Get list of coils for magnetic PDE
    params->GetCoilList( coilName_, pdename_ );

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
    params->GetList( "name", magnetsDomain_, pdename_, "magnets" );

    if ( magnetsDomain_.GetSize() > 0 ) {

      Info->PrintF( pdename_,
		    " Found permanent magnets in the following regions:" );

      Double tmpDir[3];

      // for each magnet ...
      for ( UInt k = 0; k < magnetsDomain_.GetSize(); k++ ) {

	// ... report name to logfile
	Info->PrintF( pdename_, " %s", magnetsDomain_[k].c_str() );

	// ... read direction of magnetisation
	params->CGet( "orientX", tmpDir[0], "name", magnetsDomain_[k], true,
		      pdename_, "magnets" );
	params->CGet( "orientX", tmpDir[1], "name", magnetsDomain_[k], true,
		      pdename_, "magnets" );
	params->CGet( "orientX", tmpDir[2], "name", magnetsDomain_[k], true,
		      pdename_, "magnets" );
	magnetsOriX_.Push_back( tmpDir[0] );
	magnetsOriY_.Push_back( tmpDir[1] );
	magnetsOriZ_.Push_back( tmpDir[2] );
      }
    }
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MagPDE::ReadStoreResults() {

    ENTER_FCN( "MagPDE::ReadStoreResults" );

    // By default we only save the solution at nodal values
    savesol_    = TRUE;
    savederiv1_ = FALSE;
    savederiv2_ = FALSE;

    // -----------
    //   B-FIELD
    // -----------

    // Determine regions for which B-field must be computed
    params->CGetList( "region", calcBfield_, "type", "bfield", 0, pdename_,
		      "elemResults" );

    // If the symbolic name is "all" compute electric field for all regions
    if ( calcBfield_.GetSize() == 1 && calcBfield_[0] == "all" ) {
      calcBfield_ = subdoms_;
    }

    // Log to info file
    if ( calcBfield_.GetSize() > 0 ) {
      Info->PrintF( pdename_,
		    " Computing magnetic flux density / B-field for regions:");
      for ( Integer k = 0; k < calcBfield_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcBfield_[k].c_str() );
      }
    }

    // ----------
    //   ENERGY
    // ----------

    // Determine regions for which energy must be computed
    params->CGetList( "region", calcEnergy_, "type", "energy", 0, pdename_,
		      "elemResults" );

    // If the the symbolic name is "all" compute energy for all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == "all" ) {
      calcEnergy_ = subdoms_;
    }

    // Log to info file
    if ( calcEnergy_.GetSize() > 0 ) {
      Info->PrintF( pdename_, " Computing energy for regions:" );
      for ( Integer k = 0; k < calcEnergy_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcEnergy_[k].c_str() );
      }
    }

    // --------
    //   EDDY
    // --------

    // Determine regions for which eddy current density must be computed
    params->CGetList( "region", calcEddy_, "type", "eddy", 0, pdename_,
		      "elemResults" );

    // If the the symbolic name is "all" compute energy for all regions
    if ( calcEddy_.GetSize() == 1 && calcEddy_[0] == "all" ) {
      calcEddy_ = subdoms_;
    }

    // Log to info file
    if ( calcEddy_.GetSize() > 0 ) {
      Info->PrintF( pdename_,
		    " Computing eddy current densities for regions:" );
      for ( Integer k = 0; k < calcEddy_.GetSize(); k++ ) {
	Info->PrintF( pdename_, " %s", calcEddy_[k].c_str() );
      }
    }

  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================


#endif // for #ifndef XMLPARAMS


} // end of namespace

