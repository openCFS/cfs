#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/elem.hh"
#include "Estimator/spaceerror.hh"
#include "basePDE.hh"
#include "scalarblockEQN.hh"
#include "scalarnodeEQN.hh"
#include "blocknodeEQN.hh"
#include "superblockEQN.hh"


namespace CoupledField
{

BasePDE::BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
		 WriteResults * aOutFile, TimeFunc *aptTimeFunc)
  :nonLin_(FALSE),
   incStopCrit_(1e-2), 
   residualStopCrit_(1e-3),
   firstTimeStepStatic_(TRUE),
   isaxi_(FALSE),
   isComplex_(FALSE),
   numPDENodes_(0),
   numElems_(0),
   sol_(NULL),
   isAlwaysStatic_(FALSE),
   dampingType_(NONE),
   needsDampingMatrix_(FALSE)
{

  ENTER_FCN( "BasePDE::BasePDE" );

  nonLin_ = FALSE;
  incStopCrit_ = 1e-2;
  residualStopCrit_ = 1e-3;


  // =====================================================================
  // set file pointers
  // =====================================================================
  inFile_     = aInFile;
  outFile_    = aOutFile;
  ptTimeFunc_ = aptTimeFunc;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;
  assemble_   = NULL;
  algsys_     = NULL;
  eqnData_    = NULL;

  // =====================================================================
  // set analysis parameters
  // =====================================================================
  actlevel_ = 0;
  actFrequency_ = 0;
  complexFormat_ = AMPLITUDE_PHASE;
  couplingBCsCounter_ = 0;
  numDirichletBCs_ = 0;
  pdeIsCoupled_ = FALSE;
  updateCouplingBCs_ = FALSE;
  updateBCs_ = 0;
  dim_ = ptgrid_->GetDim();
  initMatrices_ = FALSE;
  geoUpdate_ = FALSE;


  // =====================================================================
  // set solver parameters
  // =====================================================================
  eps_         = 1.0e-8;
  dampiter_    = 1.0;
  maxnumiter_  = 500;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;


  // =====================================================================
  // set postprocessing parameters
  // =====================================================================
  saveSol_        = TRUE;
  saveDeriv1_     = FALSE;
  saveDeriv2_     = FALSE;
  saveSolHist_    = FALSE;
  saveDeriv1Hist_ = FALSE;
  saveDeriv2Hist_ = FALSE;

  
  //standard parameter for solver
#ifdef USE_OLAS
  
  precondtype_ = JACOBI;
  solvertype_  = CG;
#else
  precondtype_ = ID;
  solvertype_  = RealDirect;
#endif

 

}


  
  void BasePDE::Init(Integer bcSequenceIndex,
		     std::string  bcSequenceTag)
  {
    ENTER_FCN( "BasePDE::Init()" );
    
    bcSequenceIndex_ = bcSequenceIndex;
    bcSequenceTag_ = bcSequenceTag;

    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================
#ifndef XMLPARAMS
    conf->getsubdompde(subdoms_,pdename_);
#else
    params->GetList( "name", subdoms_, pdename_, "region" );
    Info->PrintF( pdename_, " %s lives on regions:", pdename_.c_str());
    for ( Integer k = 0; k < subdoms_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s", subdoms_[k].c_str() );
    }
#endif

 //    //allocate according algebraic system
    algsys_ = new StandardSystem();
    
#ifdef USE_OLAS
    // Get parameter and report object of OLAS
    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport();
    
    std::string parallel = "no";
#ifndef XMLPARAMS
    conf->ifget("parallel",parallel);
#else
    Info->Warning( "Parameter 'parallel' unknown to XML" );
#endif
    
    if (parallel == "yes")
      olasParams_->SetValue( "Parallel", true);
    else
      olasParams_->SetValue( "Parallel", false);
    
#endif
    
    // =====================================================================
    // Get type of analysis
    // =====================================================================
    
    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string stepString;
    std::string analysis;

#ifndef XMLPARAMS
    conf->get("analysis", analysis);
#else
    params->Get( "type", analysis, "analysis" );
#endif

    AnalysisType analysisHelp;
    StdVector<std::string> tags, analysisTypes, pdenames;
    String2Enum(analysis,analysisHelp);
    
    if (analysisHelp == STATIC ||
	isAlwaysStatic_ == TRUE) {
      isComplex_ = FALSE;
      assemble_ = new StaticAssemble(algsys_, ptgrid_);
      analysistype_ = STATIC;
    }
    else if (analysisHelp == TRANSIENT) {
      isComplex_ = FALSE;
      assemble_ = new TransientAssemble(algsys_, ptgrid_);
      analysistype_ = TRANSIENT;
      laststepcalc_ = 1;
    }    
    else if (analysisHelp == HARMONIC) {
      isComplex_ = TRUE;
      assemble_ = new HarmonicAssemble(algsys_, ptgrid_);
      analysistype_ = HARMONIC;
      //overwrite defualt solver
#ifdef USE_OLAS
	Info->Error( "Implement this branch!", __FILE__, __LINE__ );
#else
	solvertype_ = ComplexDirectSolver;
#endif
      }
    else if (analysisHelp == MULTI_SEQUENCE) {

      //std::cerr << "BasePDE::Init: In Step for multisequence" << std::endl;
      
    
      stepString = Info->GenStr(bcSequenceIndex_);
      attrVec = "", "index", "type";
      valVec = "", stepString, pdename_;
      //std::cerr << "pdename_ = " << pdename_ << std::endl;
      //std::cerr << "stepString = " << stepString << std::endl;

      keyVec = "multiSequence", "step", "pde", "analysis";
      params->Get(keyVec, attrVec, valVec, analysis);
     
      //std::cerr << "analysis = " << analysis;

      String2Enum(analysis, analysistype_);

      if (analysistype_ == STATIC){
	assemble_ = new StaticAssemble(algsys_, ptgrid_);
	isComplex_ = FALSE;
      }
      else if (analysistype_ == TRANSIENT) {
	isComplex_ = FALSE;
	assemble_ = new TransientAssemble(algsys_, ptgrid_);
	laststepcalc_ = 1;	
      }
      else {
	Error("BasePDE::Init: AnalysisType not supported", __FILE__, __LINE__);
      }
  }
    else
      Error("Analysis Type not supported",__FILE__,__LINE__);
    
    // Determine if solution is of complex type or not
    
    if (analysistype_ == HARMONIC)
      sol_ = new NodeStoreSol<Complex>;
    else
      sol_ = new NodeStoreSol<Double>;
    
    
    // =====================================================================
    // initialize adaptivity
    // =====================================================================
    
#ifndef XMLPARAMS
    if (conf->get_option("adaptspace"))
      conf->get("tolerance_space_error", tolSpaceErr_);
    else
      tolSpaceErr_ = .0;
#else
    if( params->IsSet( "adaptspace" ) )
      {
	params->Get( "tolerance_space_error", tolSpaceErr_ );
      }
    else
      {
	tolSpaceErr_ = 0;
      }
#endif 
    
    // =====================================================================
    // read in boundary conditions
    // =====================================================================
    ReadBCs();
    ReadSpecialBCs();
    numDirichletBCs_ += GetNumRestraints(actlevel_);
    
    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    InitNonLin();
    
    // =====================================================================
    // initialize EQN-object and Storeresults class
    // =====================================================================


#ifdef XMLPARAMS
    // What type of equation numbering does the user want?
    keyVec = pdename_, "solver", "matrix", "eqnNumbering";
    std::string typeOfNumbering;
    params->Get( keyVec, typeOfNumbering );

    // Assemble a system matrix with scalar complex or double entries
    if ( typeOfNumbering == "scalar" ) {
      if ( dofspernode_ == 1 ) {
	eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				      dofspernode_ );
      }
      else {
	eqnData_ = new ScalarBlockEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				       dofspernode_ );
      }
      Info->PrintF( pdename_, "Using scalar equation numbering" );
    }

    // Treat all dofs of a node together and assemble a system matrix with
    // small square matrices as entries
    else if ( typeOfNumbering == "block" ) {
      if ( dofspernode_ == 1 ) {
	Warning("dopspernode = 1, so 'block' numbering identical to 'scalar'");
	eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				      dofspernode_ );
      }
      else {
	eqnData_ = new BlockNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				     dofspernode_ );
	Info->PrintF( pdename_, "Using block equation numbering" );
      }
    }

    // Only sensible for piezoPDE, will be phased out, once direct
    // coupling is possible
    else if ( typeOfNumbering == "superBlock" ) {
      if ( pdename_ != "piezo" ) {
	Error( "superBlock numbering only implemented for Piezo PDEs!",
	       __FILE__, __LINE__ );
      }
      else {
	eqnData_ = new SuperBlockEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
				      dofspernode_);
	Info->PrintF( pdename_, "Using super-block equation numbering" );
      }
    }

#else
	
    // For conf-file always use scalar entries in system matrix
    if (dofspernode_ == 1) {
      eqnData_  = new ScalarNodeEQN(ptgrid_, ptBCs_, subdoms_, 
				    actlevel_, dofspernode_);
    } else {
      eqnData_ = new ScalarBlockEQN(ptgrid_, ptBCs_, subdoms_, 
				      actlevel_, dofspernode_);
    }

#endif

    eqnData_->SetHomoDirichletBCs(bcs_hd_, homDirichDof_);
    eqnData_->CalcMapping();
    // eqnData_->Print(*debug);
    numPDENodes_ = eqnData_->GetNumLocalNodes();
    numElems_ = eqnData_->GetNumLocalElems();
    
    // Initialize Storesolution class
    sol_->SetNumSolutions(solTypes_.GetSize());
    sol_->SetNumNodes(numPDENodes_);
    for (Integer iSol=0; iSol<solTypes_.GetSize(); iSol++) {
      sol_->SetSolutionType(solTypes_[iSol],iSol);
      sol_->SetNumDofs(solDofs_[iSol], solTypes_[iSol]);
    }
    sol_->SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    sol_->Init(); 

    
    
    // =====================================================================
    // initialize assemble object
    // =====================================================================
    assemble_->SetPtr2EQNData(eqnData_); 
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);
    //assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, subdoms_, pressSurf_);
    assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, 
				subdoms_, surfdoms_, bcSequenceTag_);
    assemble_->SetGraphType(NODEGRAPH);
#ifdef USE_OLAS
    if (isComplex_)
      assemble_->SetMatrixEntryType(OLAS::COMPLEX);
    else
      assemble_->SetMatrixEntryType(OLAS::DOUBLE);
    assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
#else
    if (eqnData_->IsBlockMapped()) {
      if (isComplex_)
	assemble_->SetMatrixType(CBLOCK);
      else
	assemble_->SetMatrixType(RBLOCK);
    }
    else {
      if (isComplex_)
	assemble_->SetMatrixType(CSCALAR);
      else
	assemble_->SetMatrixType(RSCALAR);
    }
#endif 
    assemble_->SetNumDirichlet(numDirichletBCs_);
    
    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    if (needsDampingMatrix_) 
      assemble_->NeedDampingMatrix();

//     SetAlgSys();
    
    // =====================================================================
    // read in material data
    // =====================================================================
    ReadMaterialData();
    
    // =====================================================================
    // define the integratos for PDE
    // =====================================================================
    DefineIntegrators(actlevel_);
    
    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
#ifndef XMLPARAMS
    ReadSavings();
#else
    ReadStoreResults();
#endif
  }
  
  // For XML parameter handling we have replaced this method by the pure
  // virtual ReadStoreResults() method which must be implemented by each
  // PDE according to its demands.
  
#ifndef XMLPARAMS
  
  void BasePDE::ReadSavings() {

    ENTER_FCN( "BasePDE::ReadSavings" );

    //set saving of solution to yes, if user has not used the nodalsave-command
    savesol_ = TRUE;

    //check for node saving
    StdVector<std::string> savings;
  
    //reset saving of solution, if user has used the nodalsave-command
    if (conf->ifgetliststr("nodalsave", savings, pdename_))
      savesol_ = FALSE;

    for (Integer i=0; i<savings.GetSize(); i++) {
      if (savings[i] == "dof")  savesol_ = TRUE;
      if (savings[i] == "deriv1") savederiv1_ = TRUE;
      if (savings[i] == "deriv2") savederiv2_ = TRUE;
    }
  }

#endif
  

  Integer BasePDE::GetBCDof(const std::string dofStartString)
  {
    ENTER_FCN( "BasePDEPDE::GetBCDof" );

    Integer nrActDof = 0;

    if ( dofStartString == "ux" )
      nrActDof = 1;
    if ( dofStartString == "uy" )
      nrActDof = 2;
    if ( dofStartString == "uz" )
      nrActDof = 3;
    // HARD coded for piezo case
    if ( dofStartString == "ep" )
      nrActDof = dofspernode_;
    if ( nrActDof == 0 )
      Error("Unknown dof-type in homog. BC; substring must start with ux, uy, uz or ep!!",
	    __FILE__, __LINE__);

    return nrActDof;
  }


void BasePDE::WriteGeneralPDEdefines()
{

  ENTER_FCN( "BasePDE::WriteGeneralPDEdefines" );

  //BC-section
  for (Integer i=0; i< bcs_hd_.GetSize(); i++) 
    {
      Integer dof;
      std::string doftype = bcs_hd_[i];
      if (dofspernode_ > 1) 
	dof = GetBCDof(homDirichDof_[i]);
      else
	dof = 1;

      Info->WriteHomBC(pdename_, bcs_hd_[i], dof);	
    }

  for (Integer i=0; i< bcs_id_.GetSize(); i++) 
    {
      Integer dof;
      std::string doftype = bcs_id_[i];
      if (dofspernode_ > 1) 
	dof = GetBCDof(inhomDirichDof_[i]);
      else
	dof = 1;

      Info->WriteInHomBC( pdename_, bcs_id_[i], val_id_[i], fncnames_id_[i],
			  dof );	
    }

  // Loads
  StdVector<std::string> loadDom = GetLoadDom();
  StdVector<std::string> loadDof = GetLoadDof();
  StdVector<Double> loadVals = GetLoadVals();
  StdVector<std::string> loadfncs = GetLoadFncs();

  for(int i=0; i < loadDom.GetSize(); i++)
    {
      Integer dof;
      std::string doftype = loadDom[i];
      if (dofspernode_ > 1) 
	dof = GetBCDof(loadDof[i]);
      else
	dof = 1;

      Info->WriteLoad(pdename_, loadDom[i], loadVals[i], loadfncs[i], dof);
    }
}


  // ======================================================
  // Solve Step SECTION  
  // ======================================================

// time is used for a series of static calculations
// don't get confused with REAL transient simulations!
void BasePDE::SolveStepStatic(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset)
{
  ENTER_FCN( "BasePDE::SolveStepStatic" );

  lasttimecalc_ = asteptime;
  laststepcalc_ = kstep;


  if (nonLin_)
    StepStaticNonLin(kstep,asteptime,level,reset);
  else
    StepStaticLin(kstep,asteptime,level,reset);

}


void BasePDE::StepStaticLin(const Integer kstep, const Double aTime,
			    const Integer level, const Boolean reset)
{
  ENTER_FCN( "BasePDE::StepStaticLin" );

  Integer job;
  Double * ptsol;
  lasttimecalc_ = aTime; // for correct output in unv-file
  

  if ( pdeIsCoupled_ == FALSE || iterCoupledCounter_ == 0 && firstTimeStepStatic_)
    assemble_->AssembleMatrices(level);


  //this has to be done each time!
  assemble_->AssembleSrcRHS(level, aTime);

  const Integer numElems = numPDENodes_ * dofspernode_;

  if ( pdeIsCoupled_ == FALSE || iterCoupledCounter_ == 0)
    {
      //account for bcs
      SetBCs(level, updateBCs_, aTime);
      updateBCs_ = 0;

      if (firstTimeStepStatic_)
	job = 1; // calc new preconditioner
      else 
	job = 3;  // just set new BCs    

#ifdef USE_OLAS
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond(job);
#else
      algsys_->CalcPrecond(job);
#endif

    }
  else
    job = 3;

#ifdef USE_OLAS
  algsys_->BuildInDirichlet();
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif  

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;
  
  //sol_->SetDataPointer(ptsol);
  sol_->CopyFromAlgSysDataPointer(ptsol);

  firstTimeStepStatic_ = FALSE;
}


void BasePDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
			     const Integer level, const Boolean reset)
{

  ENTER_FCN( "BasePDE::SolveStepTrans" );

  lasttimecalc_= asteptime;
  recalc_ = FALSE;

  if (laststepcalc_ == kstep && kstep != 1) 
    recalc_=TRUE;
  else 
    laststepcalc_= kstep;

  if (nonLin_)
    StepTransNonLin(kstep, asteptime, level, reset);
  else
    StepTransLin(kstep, asteptime, level, reset);
}


void BasePDE::PreStepTrans(const Integer kstep, const Double asteptime,
			   const Integer level, const Boolean reset)
{

  ENTER_FCN( "BasePDE::PreStepTrans" );

  lasttimecalc_= asteptime;

  // due to coupling-pdes, the RHS has to be initialized BEFORE 
  // the coupling forces are assembled to the RHS
  algsys_->InitRHS();

  if (geoUpdate_)
    {
      assemble_->SetNonlinGeo();

      algsys_->InitRHS();
      algsys_->InitSol();
      assemble_->InitMatrices();

      assemble_->SetReassemble();   
    }
  
  
}



void BasePDE::PostStepTrans(const Integer kstep, const Double asteptime, const Integer level)
{
  ENTER_FCN( "BasePDE::PostStepTrans" );

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
  if (pdeIsCoupled_)
    {
      //save solution
      Vector<Double> solvector= solhelp->GetAlgSysVector();

      //perform corrector step
      TS_alg_->Corrector(solvector); 
    }
  
  if (pdeIsCoupled_)
    iterCoupledCounter_++;
  
}


// initialize PDEs before iteration (done for each time step)
 void BasePDE::InitStepTransCoupled(Double asteptime)
  {
    ENTER_FCN( "BasePDE::InitStepTransCoupled" );
    lasttimecalc_= asteptime;
    iterCoupledCounter_ = 0;
  }
 


void BasePDE::StepTransLin(const Integer kstep, const Double asteptime,
			   const Integer level, const Boolean reset)
{

  ENTER_FCN( "BasePDE::StepTransLin" );

  Double * ptsol;
  Integer update,job;
   

  //account for RHS
  assemble_->AssembleSrcRHS(level,lasttimecalc_);

  Integer k=0;
  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

  if ( pdeIsCoupled_ == FALSE || iterCoupledCounter_ == 0)
    {    
      
      Vector<Double> solvector= solhelp->GetAlgSysVector();
      TS_alg_->Predictor(solvector);
      
    }
  

  if (laststepcalc_ == 1)
    {
      update = 0;
      job = 3;

      // why is the first statement checking for 'pdeIsCoupled'?
      if ( pdeIsCoupled_ == FALSE 
	   || iterCoupledCounter_ == 0 
	   || geoUpdate_ == TRUE)
	{
	  job = 1;
	  assemble_->AssembleMatrices(level);
	  algsys_->ConstructEffectiveMatrix(matrix_factor_);
	}  
    }
  else if (reset)
    {
      update = 1;
      job    = 1;

      algsys_->InitMatrix(SYSTEM);
      algsys_->InitMatrix(STIFFNESS);
      algsys_->InitMatrix(MASS);
      if (dampingType_)
        algsys_->InitMatrix(DAMPING);

      algsys_->ConstructEffectiveMatrix(matrix_factor_); 
    }
  else
    {
      update = 1;
      job    = 3; 

      // The following section is only an experiment up to now
      if (geoUpdate_ == TRUE && pdeIsCoupled_ == TRUE){
	job = 1;
	assemble_->AssembleMatrices(level);
	algsys_->ConstructEffectiveMatrix(matrix_factor_);	
      }
    }


  TS_alg_->UpdateRHS();

  SetBCs(level,update,lasttimecalc_);
#ifdef USE_OLAS
  algsys_->BuildInDirichlet();
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif

  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  //sol_->SetDataPointer(ptsol);
  sol_->CopyFromAlgSysDataPointer(ptsol);
  

  Vector<Double> & solvector =\
    dynamic_cast<NodeStoreSol<Double>*>(sol_)->GetAlgSysVector();

  if (!pdeIsCoupled_)
    TS_alg_->Corrector(solvector);
}


void BasePDE::PreStepHarmonic(const Integer freqStep, const Double frequency,  
				const Integer level, const Boolean reset)
{

  ENTER_FCN( "BasePDE::PreStepHarmonic" );

  actFrequency_ = frequency;
  actFreqStep_  = freqStep;

  assemble_->SetFrequency(frequency);

  algsys_->InitRHS();

  if (reset)
      assemble_->InitMatrices();  
}


void BasePDE::SolveStepHarmonic(const Integer freqStep, const Double frequency,  
				const Integer level, const Boolean reset)
{
  ENTER_FCN( "BasePDE::SolveStepHarmonic" );

  if (nonLin_)
    StepHarmonicNonLin(freqStep, frequency,level,reset);
  else
    StepHarmonicLin(freqStep, frequency,level,reset);

}


void BasePDE::StepHarmonicLin(const Integer freqStep, const Double frequency, 
			      const Integer level, const Boolean reset)
{
  ENTER_FCN( "BasePDE::StepHarmonicLin" );

  Integer job;
  Double * ptsol;

  if (reset)
    assemble_->AssembleMatrices(level);


  //this has to be done each time!
  assemble_->AssembleSrcRHS(level, frequency);

  const Integer numElems = numPDENodes_ * dofspernode_;

  if (reset)
    {
      //account for bcs
      SetBCs(level, updateBCs_, frequency);
    
      job = 1; // calc new preconditioner

    }
  else
    job = 3;

#ifdef USE_OLAS
  algsys_->BuildInDirichlet();
  algsys_->SetupPrecond(job);
#else
  algsys_->CalcPrecond(job);
#endif  

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  //Vector<Complex> tmp(numPDENodes_*dofspernode_);

  //  if (dofspernode_ > 1)
  //    Error("Currenrly just dofpernode=1 supported in StepHarmonicLin");

  //Integer k=0;
  
  //for (Integer node=0; node<numPDENodes_*dofspernode_; node++)
  //{
  //  Complex val(ptsol[k],ptsol[k+1]);
  //  tmp[node] = val;
  //  k+=2;
    //}

  //sol_->SetAlgSysVector(tmp);

  sol_->CopyFromAlgSysDataPointer(ptsol);

}


void  BasePDE::SetBCs(const Integer level, const Integer update, const Double time)
{

  ENTER_FCN( "BasePDE::SetBCs" );

  //std::cerr << "BasePDE::SetBCs" << std::endl;
  //std::cerr << "===============" << std::endl;
  
  
  Integer node, dof;
  Double val, val_tfunc;

  std::list<Integer> nodes;

  Integer i;
  Integer j;
  Integer eqnNr, eqnDof;

  if (pdeIsCoupled_)
    j = couplingBCsCounter_;
  else
    j=0;

  val = 0;
  for (i=0; i<bcs_hd_.GetSize(); i++)
    {  
      dof = 1;
      if (dofspernode_ > 1)
	{
	  std::string doftype = bcs_hd_[i];
	  dof = GetBCDof(homDirichDof_[i]);
	}

      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	{
	  node=*p;
	  eqnData_->Node2EQN(node, dof, eqnNr, eqnDof);
	  if (eqnNr > 0)
	    {
	      
	      // Increment counter for BCs
	      //	      std::cerr << "   -> SET!" << std::endl;
#ifdef USE_OLAS
	      algsys_->SetDirichlet(j+1, eqnNr, val, eqnDof);
#else
	      if (analysistype_ == HARMONIC) 
		{
		  algsys_->SetDirichlet(j*2+1, eqnNr, val, eqnDof,SYSTEM);
		  // set imag part 
		  algsys_->SetDirichlet(j*2+2, eqnNr, val, eqnDof+1, SYSTEM);
		}
	      else
		algsys_->SetDirichlet(j+1, eqnNr ,val, eqnDof, SYSTEM);
#endif
	      j++;
	    }
	}
    }

  for (i=0; i<bcs_id_.GetSize(); i++)
    {
      dof = 1;
      if (dofspernode_ > 1)
	{
	  std::string doftype = bcs_id_[i]; 
	  dof = GetBCDof(inhomDirichDof_[i]);
	}
      
      nodes=ptBCs_->GetNodesLevel(bcs_id_[i]); 

      //get the correct time function value
      val_tfunc = 1.0;
      if (ptTimeFunc_->GetmaxTimeFnc() > 0 )
	  val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncnames_id_[i]);
      
      val=val_id_[i]*val_tfunc;

      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  eqnData_->Node2EQN(node, dof, eqnNr, eqnDof);
#ifdef USE_OLAS
	  algsys_->SetDirichlet(j+1, eqnNr, val, eqnDof);
#else
	  if (analysistype_ == HARMONIC) 
	    {
	      algsys_->SetDirichlet(j*2+1, eqnNr, val, eqnDof, SYSTEM);
	      // set imag part 
	      algsys_->SetDirichlet(j*2+2, eqnNr, 0.0, eqnDof+1, SYSTEM);
	    }
	  else
	    algsys_->SetDirichlet(j+1, eqnNr, val, eqnDof, SYSTEM);
#endif
	}
    }
}



void BasePDE::SetSolution(CFSVector &sol)
{
  ENTER_FCN( "BasePDE::SetSolution" );
  
  sol_->SetAlgSysVector(sol);

}


void BasePDE::ReadBCs()
{

  ENTER_FCN( " entering BasePDE::ReadBCs " );


#ifndef XMLPARAMS
  conf->ifgetliststr( "homogeneous_dirichlet"  , bcs_hd_, pdename_ ); 
  conf->ifgetliststr( "inhomogeneous_dirichlet", bcs_id_, pdename_ );
  
  val_id_.Resize(bcs_id_.GetSize());
  fncnames_id_.Resize(bcs_id_.GetSize());

  for( Integer i = 0; i < bcs_id_.GetSize(); i++ )
    {
      conf->get2( bcs_id_[i], val_id_[i], fncnames_id_[i], pdename_,
		  "bc_conditions", "inhomogeneous_dirichlet" );
    }

  if (dofspernode_ > 1)
    {
      conf->ifgetliststr("homogenBCDof", homDirichDof_, pdename_);  
      conf->ifgetliststr("inhomogenBCDof", inhomDirichDof_, pdename_);
      
      // just for consistency with old script
      conf->ifgetliststr("homoBCDof", homDirichDof_, pdename_);
      conf->ifgetliststr("homoBCdof", homDirichDof_, pdename_);
      conf->ifgetliststr("inhomoBCDof", inhomDirichDof_, pdename_);
      conf->ifgetliststr("inhomoBCdof", inhomDirichDof_, pdename_);
    }
#else

  // vectors for parameter handling
  StdVector<std::string> keyVec;
  StdVector<std::string> attrVec;
  StdVector<std::string> valVec;

  
  
   // Get names of node sets for homogeneous Dirichlet boundary conditions
  //   params->GetList( "name", bcs_hd_, pdename_, "dirichletHom"   );
  keyVec = pdename_, "bcsAndLoads", "dirichletHom", "name";
  attrVec = "", "tag", "";
  valVec = "", bcSequenceTag_, "";
  //std::cerr << "ReadBCs: bcSequenceTag = " << bcSequenceTag_ << std::endl;
  params->GetList(keyVec, attrVec, valVec, bcs_hd_);

  // Get names of node sets, values and filenames for inhomogenous
  // Dirichlet boundary conditions
  keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "name";
  params->GetList(keyVec, attrVec, valVec, bcs_id_);

  keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "value";
  params->GetList(keyVec, attrVec, valVec, val_id_);

  keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dynamics";
  params->GetList(keyVec, attrVec, valVec, fncnames_id_);

//   std::cerr << "===============================" << std::endl;
//   std::cerr << "PDENAME : " << pdename_ << std::endl;
//    std::cerr << "dirichletInhom = " << bcs_id_ << std::endl;
//    std::cerr << "dirichletHom = " << bcs_hd_ << std::endl;
//    std::cerr << "dirichletInhom-value = " << val_id_ << std::endl;
//    std::cerr << "dirichletInhom-dynamics= " << fncnames_id_ << std::endl;
//    std::cerr << std::endl;

  // Check consistency
  if ( bcs_id_.GetSize() != val_id_.GetSize() ||
       fncnames_id_.GetSize() > bcs_id_.GetSize() )
    {
      std::string errmsg = "dirichletInhom: ";
      errmsg += "#names of node sets = " + Info->GenStr(bcs_id_.GetSize());
      errmsg += ", #values = " + Info->GenStr(val_id_.GetSize());
      errmsg += ", #dynamics = " + fncnames_id_.GetSize() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

  // We need not have as many function/filenames as boundary conditions!
  for ( Integer k = fncnames_id_.GetSize(); k < bcs_id_.GetSize(); k++ )
    {
      fncnames_id_.Push_back( "none" );
    }
  if (dofspernode_ > 1)
    {
      keyVec = pdename_, "bcsAndLoads", "dirichletHom", "dof";
      params->GetList(keyVec, attrVec, valVec, homDirichDof_);
		      
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dof";
      params->GetList(keyVec, attrVec, valVec, inhomDirichDof_);
//       params->GetList( "dof", homDirichDof_  , pdename_, "dirichletHom" );  
//       params->GetList( "dof", inhomDirichDof_, pdename_, "dirichletInhom" );
    }
#endif

  // =====================================================================
  // if pde has more than one dof, initialize dof of boundary
  // conditions
  // =====================================================================
  
  if (dofspernode_ > 1)
    {
      if (bcs_hd_.GetSize() != homDirichDof_.GetSize())
	{
	  std::string errmsg = "Inconsistent definition of homogeneous ";
	  errmsg += "Dirichlet Boundary Conditions\n";
	  errmsg += " bcs_hd_.GetSize() = " + Info->GenStr( bcs_hd_.GetSize() );
	  errmsg += "\n homDirichDof_.GetSize() = "
	  + Info->GenStr( homDirichDof_.GetSize() ) + '\n';
	  Info->Error( errmsg, __FILE__, __LINE__ );
	}
      if (bcs_id_.GetSize() != inhomDirichDof_.GetSize()) 
	{
	  std::string errmsg = "Inconsistent definition of inhomogeneous ";
	  errmsg += "Dirichlet Boundary Conditions";
	  Info->Error( errmsg, __FILE__, __LINE__ );
	}
    }
  

  

}


void BasePDE::ReadMaterialData()
{
  ENTER_FCN( "BasePDE::ReadMaterialData" );

#ifndef XMLPARAMS

  // read material-file name from config-file
  std::string matFileName;
  conf->get( "material_file", matFileName );
  LoadMaterialData loadMaterial( matFileName.c_str() );

  //read material data for each subdomain
  materialData_ = new MaterialData[subdoms_.GetSize()];

  std::string matName;
  for (Integer i=0; i<subdoms_.GetSize(); i++)
    {
      // load material data into array "materialData_"
      conf->getstr(subdoms_[i], matName, "list_subdomains");
      loadMaterial.GetMaterial(materialData_[i], matName, pdematerialclass_);
    }

#else

  // Query name of file with material data
  std::string matFileName;
  params->Get( "file", matFileName, "materialData" );

  // Generate new material reader
  LoadMaterialData loadMaterial( matFileName.c_str() );

  // Allocate space to hold material data for each subdomain of this PDE
  materialData_ = new MaterialData[subdoms_.GetSize()];

  // Get list of subdomains and materials
  StdVector< std::string > subdomName;
  StdVector< std::string > subdomMaterial;
  params->GetList( "name", subdomName, "domain", "region" );
  params->GetList( "material", subdomMaterial, "domain", "region" );

  // Load material data for subdomains on which this PDE lives
  // from data file
  for( Integer i = 0; i < subdoms_.GetSize(); i++ )
    {
      for( Integer k = 0; k <= subdomName.GetSize(); k++ )
	{
	  if( subdoms_[i] == subdomName[k] )
	    {
	      loadMaterial.GetMaterial( materialData_[i], subdomMaterial[k],
					pdematerialclass_ );
	      break;
	    }
	}
    }

#endif

}


Integer BasePDE::GetNumRestraints(const Integer level)
{

  ENTER_FCN( "BasePDE::GetNumRestraints" );
    
  Integer res=0;
  Integer i;

  for (i=0; i<bcs_hd_.GetSize(); i++)
    {
      res+=ptBCs_->GetNumNodesLevel(bcs_hd_[i]);
    }

   for (i=0; i<bcs_id_.GetSize(); i++)
    {
      res+=ptBCs_->GetNumNodesLevel(bcs_id_[i]);
    }

  return res;
}





// **********************
//   Default Destructor
// **********************
BasePDE::~BasePDE()
{

  ENTER_FCN( " entering BasePDE::~BasePDE" );
  // ATTENTION: Dummy value for as_id!!!!!!!!!!!!!!!!!!!!!!!!!!
  DeleteAlgSys(0);
  
  if (assemble_)
    delete assemble_;
  
  // ATTENTION: Dummy value for as_id!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  if (sol_)
    delete sol_;

  

  if (eqnData_)
    delete eqnData_;
}






  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  void BasePDE::SetAlgSys()
  {

    ENTER_FCN( " Analysis::SetAlgSys" );

    // Set parameter for solver and preconditioner

    (*cla) <<  "--- PDE: " << pdename_ << " ---" << std::endl;

#if defined(USE_OLAS) && defined(XMLPARAMS)
    CFSOLASParams::SetParams( pdename_, params, olasParams_ );
#else
    SetSolverParameters();
#endif

    //set the graph type used for the system matrices
    //std::cerr << "BasePDE: We have " << eqnData_->GetNumEQNs();
    //std::cerr << " Equations" << std::endl;
    assemble_->SetupMatrixGraph(eqnData_->GetNumEQNs());

    //allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }



void BasePDE::SetSolverParameters()
{

  ENTER_FCN( " BasePDE::SetSolverParameters" );

  // Get Solver and Precond parameters
#ifdef USE_OLAS
#ifdef XMLPARAMS
  Error( "Use CFSOLASParams::SetSolverParams/SetPrecondParams instead",
	 __FILE__, __LINE__ );

  // The following is compatibility code for interfacing with old
  // .conf file. This is to be phased out!!!

#else

  //if values are defined in conf-file, take these

  // relative accuracy in the precond. energy
  conf->ifget("eps",eps_,pdename_);
  // damping parameter for Jacobi, SSOR
  conf->ifget("dampiter",dampiter_,pdename_);
  // maximal number of iterations
  conf->ifget("maxnumit",maxnumiter_,pdename_);
  // number of equation for coarsing
  conf->ifget("numeqcoarse",numeqcoarse_,pdename_);
  // coarsing parameter for AMG
  conf->ifget("coarsealpha",coarsealpha_,pdename_);

  olasParams_->SetValue( "eps", eps_ );
  olasParams_->SetValue( "MaxIter", maxnumiter_);
  olasParams_->SetValue( "epsmach", 1e-30 );
  olasParams_->SetValue( "Solver", solvertype_ );
  olasParams_->SetValue( "Precond", precondtype_);
  // The following parameters are not passed yet
  // -> Contact Uwe regarding multrigrid parameters
  // olasParams_->SetValue( "dampiter", dampiter_);
  // olasParams_->SetValue( "numeqcoarse", numeqcoarse_);
  // olasParams_->SetValue( "coarseAlpha", coarsealpha_);

  Integer solverIntegerVal = 0;
  Integer precondIntegerVal = 0;
  Integer matrixStorageTypeVal = 0;

  conf->ifget("solvertype",solverIntegerVal,pdename_);
  conf->ifget("precondtype", precondIntegerVal,pdename_);
  conf->ifget("matrixstoragetype",matrixStorageTypeVal,pdename_);

  // Assign correct solver
  switch(solverIntegerVal) {

  case 0:
    // nothing specified in conf-file use default value
    // set in the according specialized PDE
    break;
  case 1:
    solvertype_ = OLAS::RICHARDSON;
    break;
  case 2:
    solvertype_ = OLAS::RICHARDSON;
    break;
  case 3:
    solvertype_ = OLAS::CG;
    break;
  case 4:
    solvertype_ = OLAS::CG;
    break;
  case 5:
    solvertype_ = OLAS::LANCZOS;
    break;
  case 6:
    solvertype_ = OLAS::QMR;
    break;
  case 7:
    solvertype_ = OLAS::QMR;
    break;
  case 8:
    solvertype_ = OLAS::DIRECT;
    break;
  case 9:
    Error("The specified solver in the config file is not implemented in OLAS",
	  __FILE__,__LINE__);
    break;
  case 10:
    Error("The specified solver in the config file is not implemented in OLAS",
	  __FILE__,__LINE__);
    break;
  case 12:
    Error("The specified solver in the config file is not implemented in OLAS",
	  __FILE__,__LINE__);
    break;
  case 13:
    solvertype_ = OLAS::HYPRE_PCG;
    break;
  case 14:
    solvertype_ = OLAS::LAPACK_LU;
    break;
  default:
    Error("The specified solver in the config file is not known",__FILE__,
	  __LINE__);
    break;
  }

  // Assign correct preconditioner
  switch(precondIntegerVal) {

  case 0:
    // nothing specified in conf-file  
    // use default value set in the according specialized PDE
    break;
  case 1:
    precondtype_ = OLAS::ID;
    break;
  case 2:
    precondtype_ = OLAS::MG;
    break;
  case 3:
    precondtype_ = OLAS::JACOBI;
    break;
  case 4:
    precondtype_ = OLAS::ILU;
    break;
  case 5:
    precondtype_ = OLAS::SSOR;
    break;
  case 6:
    precondtype_ = OLAS::HYPRE_JACOBI;
    break;
  case 7:
    precondtype_ = OLAS::HYPRE_AMG;
    break;
  case 8:
    precondtype_ = OLAS::HYPRE_ILU;
  case 9:
    precondtype_ = OLAS::HYPRE_SPAI;
    break;
  default:
    Error("The specified preconditioner in the config file is not known",
	  __FILE__,__LINE__);
    break;
  }

  // Assign matrixStorageType
  switch(matrixStorageTypeVal) {
  case 0:
    // nothing specified in conf-file.
    // use default value set in the according specialized PDE
    break;
  case 1:
    assemble_->SetMatrixStorageType(OLAS::SPARSE_SYM);      
    break;
  case 2:
    assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);      
    break;
  case 3:
    assemble_->SetMatrixStorageType(OLAS::SKYLINE_SYM);      
    break;
  case 4:
    assemble_->SetMatrixStorageType(OLAS::SKYLINE_NONSYM);      
    break;
  case 5:
    assemble_->SetMatrixStorageType(OLAS::HYPRE_MATRIX);      
    break;
  case 6:
    assemble_->SetMatrixStorageType(OLAS::LAPACK_GBMATRIX);      
    break;
  case 7:
    assemble_->SetMatrixStorageType(OLAS::LAPACK_PBMATRIX);      
    break;
  default:
    Error("The matrix storage format in the config file is not known",
	  __FILE__,__LINE__);
    break;
  }

  // Adapt preconditioner in case of a direct solver
  if ( solvertype_ == DIRECT && precondtype_ != ID ) {
    precondtype_ = ID;
  }
#endif //XMLPARAMS

  // Here starts the branch for LAS
#else

#ifndef XMLPARAMS

  //if values are defined in conf-file, take these

  conf->ifget( "solvertype",  solvertype_,  pdename_ );
  conf->ifget( "precondtype", precondtype_, pdename_ );
  if ( solvertype_ == RealDirect && precondtype_ != ID ) {
    precondtype_ = ID;
  }

  // relative accuracy in the precond. energy
  conf->ifget("eps",eps_,pdename_);
  // damping parameter for Jacobi, SSOR
  conf->ifget("dampiter",dampiter_,pdename_);
  // maximal number of iterations
  conf->ifget("maxnumit",maxnumiter_,pdename_);
  // number of equation for coarsing
  conf->ifget("numeqcoarse",numeqcoarse_,pdename_);
  // coarsing parameter for AMG
  conf->ifget("coarsealpha",coarsealpha_,pdename_);
#else
  
 // Section for LAS and XMLPARAMS
  std::string solverType, errMsg;
  std::string precondType;
  StdVector<std::string> numIters;
  StdVector<Double> eps;
  
  params->Get( "type", solverType, pdename_, "solver");
  
  if (solverType == "Direct") {
    if (!isComplex_)
      solvertype_ = RealDirect;
    else
      solvertype_ = ComplexDirectSolver;
  }
  else if (solverType == "CG") {
    if (!isComplex_)
      solvertype_ = RealCG;
    else
      solvertype_ = ComplexCG;
  }
  else if (solverType == "Richardson") {
    if (!isComplex_)
      solvertype_ = RealRichardson;
    else
      solvertype_ = ComplexRichardson;
  }
  else if (solverType == "QMR") {
    if (!isComplex_)
      solvertype_ = RealQMR;
    else
      solvertype_ = ComplexQMR;
  }
  else {
    errMsg = "Solvertype '";
    errMsg += solverType;
    errMsg += "' is not known for PDE '";
    errMsg += pdename_;
    errMsg += "' !";
    Error(errMsg.c_str(), __FILE__, __LINE__);
  }
  
  
  
  params->Get( "precond", precondType, pdename_, "solver");
  if (precondType == "noPrecond")
    precondtype_ = ID;
  else if (precondType == "Id")
    precondtype_ = ID;
  else if (precondType == "Jacobi")
    precondtype_ = JACOBI;
  else if (precondType == "SSOR")
    precondtype_ = SSOR;
  else if (precondType == "ILU")
    precondtype_ = ILU;
  else if (precondType == "MG") {
    precondtype_ = MG;
    Warning("No addtional MG parameters available", __FILE__, __LINE__);
  }
  else {
    
    errMsg = "Precondtype '";
    errMsg += precondType;
    errMsg += "' is not known for PDE '";
    errMsg += pdename_;
    errMsg += "' !";
    Error(errMsg.c_str(), __FILE__, __LINE__);
  }
    
  params->GetList( "tol", eps, pdename_ );
  if (eps.GetSize() == 1) {
    eps_ = eps[0];
  }
  params->GetList( "maxIter", numIters, pdename_);
  if (numIters.GetSize() == 1) {
    maxnumiter_ = atoi(numIters[0].c_str());
  }
  
#endif
  //communicate with algebraic system
  algsys_->CreateParameter();
  algsys_->SetAccuracy(eps_);
  algsys_->SetMaxNumIter(maxnumiter_);
  algsys_->SetPrecond(precondtype_);
  algsys_->SetSolver(solvertype_);
  algsys_->SetDampIter(dampiter_);
  algsys_->SetCoarseSystem(numeqcoarse_);
  algsys_->SetAlpha(coarsealpha_);

#endif

} 




void BasePDE::CreateMatrices_Solver()
{

  ENTER_FCN( " BasePDE::CreateMatrices_Solver" );

  // create and initialize matrices 
  assemble_->CreateMatrices();  
  assemble_->InitMatrices();  

  //create solver and preconditioner
  algsys_->CreateSolver();

#ifdef USE_OLAS
  algsys_->CreatePrecond();
#else
  algsys_->CreatePrecond(assemble_->GetMatrixType());
#endif

  //now reset AlgebraicSystem 
  algsys_->InitRHS();
  algsys_->InitSol();
}


  // ======================================================
  // COUPLING SECTION
  // ======================================================


void BasePDE::CalcInputCoupling()
{

  ENTER_FCN( " BasePDE::CalcInputCoupling" );

  StdVector<Integer> * nodes;
  StdVector<Elem*> * elements;
  CFSVector * val;
  Double * help;
  Integer pdeNode, eqnNr,eqnDof, helpnode;
  Integer couplingDof;
  
  // Reset counter for boundary conditions
  couplingBCsCounter_ = 0;
  
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      couplingDof = ptCoupling_->GetInputDof(i);
    
      // Up to now, Coupling is only possible with
      // Real valued solutions
      Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);

      switch(ptCoupling_->GetInputType(i))
	{
	  
	case COORD:
	  initMatrices_ = TRUE;
	  ptCoupling_->GetInputNodes(i, nodes);
	  // -- OLD --
	  //deltCoords_.reshape(Dim_, numPDENodes_);

	  deltCoords_.Resize(dim_, numPDENodes_);
	  
	  // set ptr of deltCoords to assembly-object
	  assemble_->SetPtrDeltaCoordinates(&deltCoords_);
	  
	  for (Integer j=0; j<nodes->GetSize(); j++)
	    for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	      {
		//std::cerr << "processing dim = " << dim << ", j = " << j << std::endl;
		
		pdeNode = eqnData_->Mesh2PDENode((*nodes)[j]);
		if (pdeNode==-1)
		  Error("Coupling node not in my subdomains. See mesh- and config-file for errors",
			__FILE__,__LINE__);
		
		deltCoords_(dof,pdeNode-1) = help[dof + j*dim_];
	      }
	  
	  break;

	case RHS:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(RHS)" << std::endl;
	  ptCoupling_->GetInputNodes(i, nodes);
	  

	  //for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	  for (Integer dof=0; dof<couplingDof; dof++)
	    for (Integer j=0; j<nodes->GetSize(); j++)
	      {
		// std::cerr << "Getting input node " << (*nodes)[j] << std::endl;
		eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);
		// This warning is disabled, in multi-dof pdes
		// only one compomemt
		// if (eqnNr==0)
// 		  {
// 		    // std::cerr << "PDENODE: "  << PDEnode << "Node[" << (*nodes)[j] << "][" 
// 		    //      << dof+1 << "]= " << (*val)[dof][j] << std::endl; 
// 		    std::string msg = pdename_;
// 		    helpnode = (*nodes)[j];
// 		    msg += ": A node with node nr ";
// 		    //msg+= helpnode;
// 		    msg+= " for Right hand side coupling is a homogenous Dirichlet node.";
// 		    msg+= "Check your config/mesh file!";
// 		    Warning(msg.c_str()  , __FILE__,__LINE__);
// 		  }
// 		else {
// 		std::cerr << "RHS: Setting EQN " << pd << " value " << help[dof+couplingDof*j] << std::endl;
		  // PROBLEM !!!!
		  // SetNodeRHS erwartet Double* oder Complex*, aber Inhalt erst zu Laufzeit bekannt ...
		if (eqnNr != 0) 
		  {
		  algsys_->SetNodeRHS(help[dof+couplingDof*j], eqnNr, eqnDof);
		  }
	      }
	  
	  break;

	case ID_BC:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(ID_BC)" << std::endl;
 	  ptCoupling_->GetInputNodes(i, nodes);
	  //std::cerr << "found " << nodes->GetSize() << " IDBCs input nodes" << std::endl;
	  for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	    for (Integer j=0; j<nodes->GetSize(); j++, couplingBCsCounter_++)
	      {
		eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);
		if (pdeNode==0)
		  Error("The specified coupling node has no equation number"
			, __FILE__,__LINE__);
#ifdef USE_OLAS
		algsys_->SetDirichlet(couplingBCsCounter_+1, eqnNr, help[dof+j*couplingDof], eqnDof);
#else
		if (updateCouplingBCs_)
		  {
		    //		    std::cerr << "updating BC[" << dof << "][" << (*nodes)[j] << "] = " 
		    //		      << (*val)[dof][j] << std::endl;
		    algsys_->UpdateDirichlet(couplingBCsCounter_+1, help[dof+j*couplingDof], SYSTEM);
		  }
		else
		  {	
		    //  std::cerr << "BC[" << dim << "][" << (*nodes)[j] << "] = " << (*val)[dim][j] << std::endl;
		    algsys_->SetDirichlet(couplingBCsCounter_+1, eqnNr, help[dof+j*couplingDof], eqnDof, SYSTEM);
		  }
#endif
	      }
	  break;
	  
	case MAT:
	  Error("Not implemented yet",__FILE__,__LINE__);
	  break;

	}  // end switch
    } // end for

  updateCouplingBCs_ = TRUE;

}




  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================



void BasePDE::GetElemCoords(const StdVector<Integer> connect, 
			    Matrix<Double> &coordMat, 
			    const Integer level)
{

  ENTER_FCN( "BasePDE::GetElemCoords" );
  Integer pdeNode;

  ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
  if (deltCoords_.GetSizeRow() != 0 && geoUpdate_ == TRUE)
    {
      for (Integer i=0; i<coordMat.GetSizeRow(); i++)
	for (Integer j=0; j<coordMat.GetSizeCol(); j++) {
	  pdeNode = eqnData_->Mesh2PDENode(connect[j]);
	  //coordMat(i,j) += deltCoords_(i,mesh2PDENode_ [connect[j] - 1] - 1);
	  coordMat(i,j) += deltCoords_(i, pdeNode - 1);
	}
    }
}


  // ======================================================
  // ADAPTIVITY SECTION 
  // ======================================================


#ifdef ADAPTGRID
void BasePDE::RefineMesh(const Integer level)
{

  ENTER_FCN( "BasePDE::RefineMesh" );
  
  Integer                  numChilds;
  Integer                  numElems;
  StdVector<Elem*>       elemssd;
  Double                   theta_s;
  Double                   coarseFactor;
  Double                   tol4Elm;
  Integer                  numRefLoops=0;
  Integer                  numRefinements;
  Integer                  iem;

  // get max num elements for the domain,on which we have the equation
  numElems=ptgrid_->GetMaxnumElem(level,subdoms_);

  // get pointer to array with elements
  ptgrid_->GetElemSD(elemssd,subdoms_[0],level);

  // if element is marked, then value of the array's element is equal 1, else 0.
  markingElems_.Resize(numElems);

  if (!conf->ifget("safety_factor_for_space_adaptivity",theta_s))
    theta_s=1.0;

  coarseFactor = 0.7;
  tol4Elm = theta_s*tolSpaceErr_;
  
  std::cout << " tolerance space error: " << tolSpaceErr_ << std::endl;
 
  for (iem=0; iem<numElems; iem++) // loop over elements
    {
      elemssd[iem]->refinementFlag = 0;
      if (errorMap_[iem]>tol4Elm)
	{
	  elemssd[iem]->refinementFlag = 1; 

	  numChilds=elemssd[iem]->ptElem->getNumChilds();

	  numRefinements = defineRefinements(errorMap_[iem],tol4Elm,numChilds);

	  if (numRefinements>numRefLoops)
	    numRefLoops = numRefinements;

	    markingElems_[iem]=numRefinements;
	    elemssd[iem]->refinementNumber = numRefinements;
	}
      else
	if (errorMap_[iem] < coarseFactor*tol4Elm)
	    elemssd[iem]->refinementFlag = 0; /// !!!!
	else
	  elemssd[iem]->refinementFlag = 0;
    }
    
  // Fuehren die Verfeinerung durch
  std::cout << " number of refinement loops \n" << numRefLoops << std::endl;

  ptgrid_->Refine(numRefLoops);
  //    ptgrid_->RefineUniform();
  
}




Boolean BasePDE::TestError(const Integer level)
{

  ENTER_FCN( "BasePDE::TestError" );

  if (!ptError_) ConstructorError();
  
  // Berechnung der Fehlerkarte
  Double            totalErr;
  Vector<Double>    solVec;
  
  solVec.Resize(sol_.size());
  int i;
  for (i=0; i<sol_.size(); i++)
    solVec[i]=sol_[0][i];
  //      sol_.toVector(solVec,1);
  
  ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);
  
  std::cout << " space error: " << totalErr <<
    " tolerance: " << tolSpaceErr_ << std::endl;

  if (totalErr > tolSpaceErr_) return TRUE;
  else return FALSE;
  
}



//In this fnc we delete old pointer to Error-object & create new
void BasePDE::ConstructorError()
{

  ENTER_FCN( "BasePDE::ConstructorError" );

  if (ptError_) delete ptError_;
  
  ptError_=new SpaceErrorEstimator();
  ptError_->Init(this);

}




void BasePDE::WriteErrorInfo(WriteResults * ptmeshes)
{
  ptmeshes->WriteElemSolution(errorMap_,0,0,"ERR-errorMap");
  ptmeshes->WriteElemSolution(markingElems_,0,0,"ERR-markedElems");
}
#endif

void BasePDE::GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connecth)
{

  ENTER_FCN( "BasePDE::GetDerivSolVecOfElement" );

  // displacement of element nodes
  sol.Resize(dofspernode_ * connecth.GetSize());
  sol.Init(0);
  Integer eqnNr, eqnDof;
  Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();
  
  if (analysistype_ == TRANSIENT) {
    const Vector<Double> & sol_der1 = getS1();
    
    for(Integer actNode=0; actNode<connecth.GetSize(); actNode++)
      for(Integer actDof=0; actDof < dofspernode_; actDof++)
	{
	  eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
	  if (eqnNr!= 0)
	    //sol[actDof + actNode*dofspernode_] = sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];
	    sol[actDof + actNode*dofspernode_] = sol_der1[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
	  else
	    sol[actDof + actNode*dofspernode_] = 0.0;
	}
  }
}

void BasePDE::GetDerivSolOfElement(Matrix<Double>& sol, StdVector<Integer>& connect_PDE)
{

  ENTER_FCN( "BasePDE::GetDerivSolOfElement" );

  // displacement of element nodes
  sol.Resize(dofspernode_, connect_PDE.GetSize());

  const Vector<Double>& sol_der1 = getS1();  

  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof][actNode] =  sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];

}







void BasePDE::CalcLineNormalVec(Vector<Double>& n, 
				const Elem& interfaceElem, 
				const Elem& neighbour)
{

ENTER_FCN( "BasePDE::CalcLineNormalVec" );

  const Integer nrVecElem2d = 2;

  // get coords of interface elemmt ========================================
  BaseFE * ptElem = interfaceElem.ptElem;
  StdVector<Integer> connecth = interfaceElem.connect;
  
  Matrix<Double> ptCoord; 
  GetElemCoords(connecth, ptCoord, actlevel_);

  // calculates normal to element ==========================================
  CalcLineNormalVec(n, ptCoord);  // normal vector is a rotation of 90° in pos. math. direction
                                  // of the vector from node 1 to node 2 in ptCoord


  connecth = neighbour.connect;
  
  Integer indexNode1=-1;
  Integer indexNode2=-1;
  
  for(Integer actNode=0; actNode < connecth.GetSize(); actNode++)
    {
      if (connecth[actNode] == interfaceElem.connect[0])
	indexNode1 = actNode;
      if (connecth[actNode] == interfaceElem.connect[1])
	indexNode2 = actNode;
    }

  if (indexNode1==-1 || indexNode2==-1)
    Error("Nodes of neighbouring element not found!", __FILE__, __LINE__);


  // counterclockwise orientation of nodes (difference of node indizes is +1)
  if (indexNode2-indexNode1 == 1 || 
      (indexNode2-indexNode1)+connecth.GetSize() == 1 )
    n *= -1;
  
  else
    // if not clockwise orientation of nodes (difference of node indizes is -1)
    if (! (indexNode2-indexNode1 == -1 || 
	   (indexNode2-indexNode1)-connecth.GetSize()==-1) )
      Error("Nodes of interface don't lie beneath each other in neighbouring element!", __FILE__, __LINE__);  
}





// normal of line element: ATTENTION no defined sign!!
void BasePDE::CalcLineNormalVec(Vector<Double>& n, Matrix<Double>& ptCoord)
{

  ENTER_FCN( "BasePDE::CalcLineNormalVec" );

  const Integer nrVecElem2d = 2;
  
  if (ptCoord.GetSizeRow()!=nrVecElem2d)
    Error("Calc element normal: no line element! ", __FILE__,__LINE__);

  n.Resize(nrVecElem2d);
  
  // normal of a vector: interchange x-coord and y-coord and take the new x-coord as negative

  n[0] = -(ptCoord[1][1] - ptCoord[1][0]);
  n[1] =  (ptCoord[0][1] - ptCoord[0][0]);
  n /= n.NormL2();
}


void BasePDE::sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted){
  ENTER_FCN( "PiezoPDE::SortStresses" );

  //soring according to capa (unv) notation
  //our notation is Voit: Txx Tyy Tzz Tyz Txz Txy

  if (subType_ == "3d") {
    //Txx Txy Tyy Txz Tyz Tzz
    sorted[0] = unsorted[0];
    sorted[1] = unsorted[5];
    sorted[2] = unsorted[1];
    sorted[3] = unsorted[4];
    sorted[4] = unsorted[3];
    sorted[5] = unsorted[2];
  }
  else if (subType_ == "axi") {

    //Tphiphi 0 Trr 0 Trz Tzz
    sorted[0] = unsorted[3];
    sorted[1] = 0.0;
    sorted[2] = unsorted[0];
    sorted[3] = 0.0;
    sorted[4] = unsorted[2];
    sorted[5] = unsorted[1];
  }
  else {
     //0 0 Tyy 0 Tyz Tzz
    sorted[0] = 0.0;
    sorted[1] = 0.0;
    sorted[2] = unsorted[0];
    sorted[3] = 0.0;
    sorted[4] = unsorted[2];
    sorted[5] = unsorted[1];
  }

}


} // end of namespace

