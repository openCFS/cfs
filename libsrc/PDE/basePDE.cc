#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/elem.hh"
#include "Estimator/spaceerror.hh"
#include "basePDE.hh"


namespace CoupledField
{

BasePDE::BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
		 WriteResults * aOutFile, TimeFunc *aptTimeFunc)
  :nonLin_(FALSE),
   incStopCrit_(1e-2), 
   residualStopCrit_(1e-3),
   firstTimeStepStatic_(TRUE),
   isaxi_(FALSE),
   isComplex_(FALSE)
{

  ENTER_FCN( "BasePDE::BasePDE" );

  nonLin_ = FALSE;
  incStopCrit_ = 1e-2;
  residualStopCrit_ = 1e-3;

  inFile_     = aInFile;
  outFile_    = aOutFile;
  ptTimeFunc_ = aptTimeFunc;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;

  actlevel_ = 0;
  couplingBCsCounter_ = 0;
  numDirichletBCs_ = 0;
  pdeIsCoupled_ = FALSE;
  updateCouplingBCs_ = FALSE;
  updateBCs_ = 0;
  dim_ = ptgrid_->GetDim();
  initMatrices_ = FALSE;



  eps_         = 1.0e-8;
  dampiter_    = 1.0;
  maxnumiter_  = 500;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;

  savederiv1_ = FALSE;
  savederiv2_ = FALSE;
  

  //standard parameter for solver
#ifdef USE_OLAS
  
  precondtype_ = JACOBI;
  solvertype_  = CG;
#else
  precondtype_ = ID;
  solvertype_  = RealDirect;
#endif

  //get analysis type
  std::string analysis;
#ifndef XMLPARAMS
  conf->get("analysis", analysis);
#else
  params->Get( "type", analysis, "analysis" );
#endif

  //allocate according algebraic system
  algsys_ = new StandardSystem();


#ifdef USE_OLAS
  olasParams_ = algsys_->GetOLASParams();
  olasReport_ = algsys_->GetOLASReport();
  
  std::string parallel = "no";
  conf->ifget("parallel",parallel);
  
  if (parallel == "yes")
    olasParams_->SetValue( "Parallel", true);
  else
    olasParams_->SetValue( "Parallel", false);
  
#endif


  if (analysis=="static") 
    {
      assemble_ = new StaticAssemble(algsys_, ptgrid_);
      analysistype_ = STATIC;
    }

  else if (analysis=="transient")
    {
      assemble_ = new TransientAssemble(algsys_, ptgrid_);
      analysistype_ = TRANSIENT;
    }

  else if (analysis=="harmonic")
    analysistype_ = HARMONIC;

  else
    Error("Analysis Type not supported",__FILE__,__LINE__);

std::cout << "Analaysis: " << analysistype_ << std::endl;

  // Determine if solution is of complex type or not
  if (analysistype_ == HARMONIC)
    sol_ = new StoreSol<Complex>;
  else
    sol_ = new StoreSol<Double>;


  //for adaptivity

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
    std::vector<std::string> savings;
  
    //reset saving of solution, if user has used the nodalsave-command
    if (conf->ifgetliststr("nodalsave", savings, pdename_))
      savesol_ = FALSE;

    for (Integer i=0; i<savings.size(); i++) {
      if (savings[i] == "dof")  savesol_ = TRUE;
      if (savings[i] == "deriv1") savederiv1_ = TRUE;
      if (savings[i] == "deriv2") savederiv2_ = TRUE;
    }
  }

#endif
  


void BasePDE::WriteGeneralPDEdefines()
{

  ENTER_FCN( "BasePDE::WriteGeneralPDEdefines" );

  //BC-section
  for (Integer i=0; i< bcs_hd_.size(); i++) 
    {
      Integer dof;
      std::string doftype = bcs_hd_[i];
      if (dofspernode_ > 1) 
	dof = GetBCDof(homDirichDof_[i]);
      else
	dof = 1;

      Info->WriteHomBC(pdename_, bcs_hd_[i], dof);	
    }

  for (Integer i=0; i< bcs_id_.size(); i++) 
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
  std::vector<std::string> loadDom = GetLoadDom();
  std::vector<std::string> loadDof = GetLoadDof();
  std::vector<Double> loadVals = GetLoadVals();
  std::vector<std::string> loadfncs = GetLoadFncs();

  for(int i=0; i < loadDom.size(); i++)
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

  //  PreStepStatic(level);

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
  sol_->CopyFromDataPointer(ptsol);

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
  assemble_->AssembleSrcRHS(level,lasttimecalc_);
  
}



void BasePDE::PostStepTrans(const Integer kstep, const Double asteptime, const Integer level)
{
  ENTER_FCN( "BasePDE::PostStepTrans" );

  StoreSol<Double> * solhelp = dynamic_cast<StoreSol<Double>*>(sol_);
    
  if (pdeIsCoupled_)
    {
      //save solution
      Vector<Double> solvector= solhelp->GetCompleteVector();

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
   
  Integer k=0;
  StoreSol<Double> * solhelp = dynamic_cast<StoreSol<Double>*>(sol_);

  if ( pdeIsCoupled_ == FALSE || iterCoupledCounter_ == 0)
    {    
      
      Vector<Double> solvector= solhelp->GetCompleteVector();
      TS_alg_->Predictor(solvector);
      
    }
  

  if (laststepcalc_ == 1)
    {
      update = 0;
      job = 3;

      if ( pdeIsCoupled_ == FALSE || iterCoupledCounter_ == 0)
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
      if (damping_type_)
        algsys_->InitMatrix(DAMPING);

      algsys_->ConstructEffectiveMatrix(matrix_factor_); 
    }
  else
    {
      update = 1;
      job    = 3; 
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
  sol_->CopyFromDataPointer(ptsol);
  

  Vector<Double> & solvector =\
    dynamic_cast<StoreSol<Double>*>(sol_)->GetCompleteVector();

  if (!pdeIsCoupled_)
    TS_alg_->Corrector(solvector);
}

void  BasePDE::SetBCs(const Integer level, const Integer update, const Double time)
{

  ENTER_FCN( "BasePDE::SetBCs" );

  Integer node, dof;
  Double val, val_tfunc;

  std::list<Integer> nodes;

  Integer i;
  Integer j;
  if (pdeIsCoupled_)
    j = couplingBCsCounter_;
  else
    j=0;

  val = 0;
  for (i=0; i<bcs_hd_.size(); i++)
    {  
      dof = 1;
      if (dofspernode_ > 1)
	{
	  std::string doftype = bcs_hd_[i];
	  dof = GetBCDof(homDirichDof_[i]);
	}

      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
#ifdef USE_OLAS
	  algsys_->SetDirichlet(j+1, mesh2PDENode_[node-1],val, dof);
#else
	  algsys_->SetDirichlet(j+1, mesh2PDENode_[node-1],val, dof, SYSTEM);
#endif
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
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
	  // Mesh node numbers are mapped to PDE node numbers
#ifdef USE_OLAS
	  algsys_->SetDirichlet(j+1, mesh2PDENode_[node-1],
				val, dof);
#else
	  algsys_->SetDirichlet(j+1, mesh2PDENode_[node-1],
				val, dof, SYSTEM);
#endif
	}
    }
}




void BasePDE::ReadBCs( const std::string pde )
{

  ENTER_FCN( " entering BasePDE::ReadBCs " );


#ifndef XMLPARAMS
  conf->ifgetliststr( "homogeneous_dirichlet"  , bcs_hd_, pde ); 
  conf->ifgetliststr( "inhomogeneous_dirichlet", bcs_id_, pde );

  val_id_.resize(bcs_id_.size());
  fncnames_id_.resize(bcs_id_.size());

  for( Integer i = 0; i < bcs_id_.size(); i++ )
    {
      conf->get2( bcs_id_[i], val_id_[i], fncnames_id_[i], pde,
		  "bc_conditions", "inhomogeneous_dirichlet" );
    }

#else

  // Get names of node sets for homogeneous Dirichlet boundary conditions
  params->GetList( "name", bcs_hd_, pde, "dirichletHom"   );

  // Get names of node sets, values and filenames for inhomogenous
  // Dirichlet boundary conditions
  params->GetList( "name"    , bcs_id_     , pde, "dirichletInhom" );
  params->GetList( "value"   , val_id_     , pde, "dirichletInhom" );
  params->GetList( "dynamics", fncnames_id_, pde, "dirichletInhom" );

  // Check consistency
  if ( bcs_id_.size() != val_id_.size() ||
       fncnames_id_.size() > bcs_id_.size() )
    {
      std::string errmsg = "dirichletInhom: ";
      errmsg += "#names of node sets = " + Info->GenStr(bcs_id_.size());
      errmsg += ", #values = " + Info->GenStr(val_id_.size());
      errmsg += ", #dynamics = " + fncnames_id_.size() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

  // We need not have as many function/filenames as boundary conditions!
  for ( Integer k = fncnames_id_.size(); k < bcs_id_.size(); k++ )
    {
      fncnames_id_.push_back( "none" );
    }
#endif

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
  materialData_ = new MaterialData[subdoms_.size()];

  std::string matName;
  for (Integer i=0; i<subdoms_.size(); i++)
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
  materialData_ = new MaterialData[subdoms_.size()];

  // Get list of subdomains and materials
  std::vector< std::string > subdomName;
  std::vector< std::string > subdomMaterial;
  params->GetList( "name", subdomName, "domain", "region" );
  params->GetList( "material", subdomMaterial, "domain", "region" );

  // Load material data for subdomains on which this PDE lives
  // from data file
  for( Integer i = 0; i < subdoms_.size(); i++ )
    {
      for( Integer k = 0; k <= subdomName.size(); k++ )
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

  for (i=0; i<bcs_hd_.size(); i++)
    {
      res+=ptBCs_->GetNumNodesLevel(bcs_hd_[i]);
    }

  for (i=0; i<bcs_id_.size(); i++)
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

  if (sol_)
    delete sol_;
}






  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  void BasePDE::SetAlgSys(int sysid)
  {

    ENTER_FCN( " Analysis::SetAlgSys" );

    //set solver parameters  
    SetSolverParameters();

    //ask the PDE matrix settings
    //    assemble_->MatrixSettings();

    //set the graph type used for the system matrices
    assemble_->SetupMatrixGraph(numPDENodes_);
    //    assemble_->SetupMatrixGraph(numPDENodes_, NODEGRAPH);

    //allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }



void BasePDE::SetSolverParameters()
{

  ENTER_FCN( " BasePDE::SetSolverParameters" );

  // Get Solver and Precond parameters
#ifdef USE_OLAS
  
  Integer solverIntegerVal = 0;
  Integer precondIntegerVal = 0;
  Integer matrixStorageTypeVal = 0;
  

#ifndef XMLPARAMS
  conf->ifget("solvertype",solverIntegerVal,pdename_); // solver
  conf->ifget("precondtype", precondIntegerVal,pdename_); //preconditioner
  conf->ifget("matrixstoragetype",matrixStorageTypeVal,pdename_); // matrixStorageTypeVal
#else
#ifdef DEBUG
  Info->Warning("SetSolverParameters: Using defaults! No solvers yet in XML!");
#endif
#endif


  // Assign correct solver
  switch(solverIntegerVal) 
    {
    case 0:
      // nothing specified in conf-file  
      // use default value set in the according specialized PDE
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
       Error("The specified solver in the config file is not implemented in OLAS",__FILE__,__LINE__);
       break;
    case 10:
      Error("The specified solver in the config file is not implemented in OLAS",__FILE__,__LINE__);
      break;
    case 12:
      Error("The specified solver in the config file is not implemented in OLAS",__FILE__,__LINE__);
      break;
    case 13:
      solvertype_ = OLAS::HYPRE_PCG;
      break;
    case 14:
      solvertype_ = OLAS::LAPACK_LU;
      break;
    default:
      Error("The specified solver in the config file is not known",__FILE__,__LINE__);
      break;
 }


 // Assign correct preconditioner
  switch(precondIntegerVal) 
    {
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
      Error("The specified preconditioner in the config file is not known",__FILE__,__LINE__);
      break;
    }

  // Assign matrixStorageType
  switch(matrixStorageTypeVal)
    {
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
      Error("The specified matrix storage format in the config file is not known",__FILE__,__LINE__);
      break;
    }
#else
#ifndef XMLPARAMS
    conf->ifget("solvertype",solvertype_,pdename_); // solver
    conf->ifget("precondtype", precondtype_,pdename_); //preconditioner
#else
#ifdef DEBUG
  Info->Warning("SetSolverParameters: Using defaults! No solvers yet in XML!");
#endif
#endif
#endif


  //if values are defined in conf-file, take these
#ifndef XMLPARAMS
  conf->ifget("eps",eps_,pdename_); // relative accuracy in the precond. energy
  conf->ifget("dampiter",dampiter_,pdename_); // damping parameter for Jacobi, SSOR
  conf->ifget("maxnumit",maxnumiter_,pdename_); // maximal number of iterations
  conf->ifget("numeqcoarse",numeqcoarse_,pdename_); // number of equation for coarsing
  conf->ifget("coarsealpha",coarsealpha_,pdename_); // coarsing parameter for AMG
#else
#ifdef DEBUG
  Info->Warning("SetSolverParameters: Using defaults! No solvers yet in XML!");
#endif
#endif

#ifdef USE_OLAS
  if (solvertype_==DIRECT && precondtype_!=ID)  precondtype_=ID;
#else
  if (solvertype_==RealDirect && precondtype_!=ID)  precondtype_=ID;
#endif

  //communicate with algebraic system
 

#ifdef USE_OLAS
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
#else
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

  std::vector<Integer> * nodes;
  std::vector<Elem*> * elements;
  BaseStoreSol * val;
  Double * help;
  Integer PDEnode;
  Integer couplingDof;
  
  // Reset counter for boundary conditions
  couplingBCsCounter_ = 0;
  
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      couplingDof = ptCoupling_->GetInputDof(i);
    
      // Only for testing, until a common vectorclass
      // for CFS++ and OLAS is found
      help=val->GetDoublePointer();
      
      switch(ptCoupling_->GetInputType(i))
	{
	  
	case COORD:
	  initMatrices_ = TRUE;
	  ptCoupling_->GetInputNodes(i, nodes);
	  // -- OLD --
	  //deltCoords_.reshape(Dim_, numPDENodes_);

	  val->GetDataPointer(help);
	  deltCoords_.Resize(dim_, numPDENodes_);
	  
	  // set ptr of deltCoords to assembly-object
	  assemble_->SetPtrDeltaCoordinates(&deltCoords_);
	  
	  
	  for (Integer j=0; j<nodes->size(); j++)
	    for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	      {
		//std::cerr << "processing dim = " << dim << ", j = " << j << std::endl;
		PDEnode = mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  Error("Coupling node not in my subdomains. See mesh- and config-file for errors",
			__FILE__,__LINE__);
		
		deltCoords_(dof,PDEnode-1) = help[dof + j*dim_];
	      }
	  
	  break;

	case RHS:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(RHS)" << std::endl;
	  ptCoupling_->GetInputNodes(i, nodes);
	  
	  

	  //for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	  for (Integer dof=0; dof<couplingDof; dof++)
	    for (Integer j=0; j<nodes->size(); j++)
	      {
		PDEnode = mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  {
		    // std::cerr << "PDENODE: "  << PDEnode << "Node[" << (*nodes)[j] << "][" 
		    //      << dof+1 << "]= " << (*val)[dof][j] << std::endl; 
		    Error("Node not assigned to coupling domain: see mesh- and config-file"
			  , __FILE__,__LINE__);
		  }

		// PROBLEM !!!!
		// SetNodeRHS erwartet Double* oder Complex*, aber Inhalt erst zu Laufzeit bekannt ...
		algsys_->SetNodeRHS(help[dof+couplingDof*j], PDEnode, dof+1);
	      }
	  
	  break;

	case ID_BC:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(ID_BC)" << std::endl;
 	  ptCoupling_->GetInputNodes(i, nodes);
	  //std::cerr << "found " << nodes->size() << " IDBCs input nodes" << std::endl;
	  for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	    for (Integer j=0; j<nodes->size(); j++, couplingBCsCounter_++)
	      {
		PDEnode = mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);
#ifdef USE_OLAS
		algsys_->SetDirichlet(couplingBCsCounter_+1, PDEnode, help[dof+j*couplingDof], dof+1);
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
		    algsys_->SetDirichlet(couplingBCsCounter_+1, PDEnode, help[dof+j*couplingDof], dof+1, SYSTEM);
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


void BasePDE::Mesh2PDENode(Vector<Integer> & PDENodes, 
			   const Vector<Integer> & MeshNodes,
			   const std::vector<Integer> & Mesh2PDENode)
{

  ENTER_FCN( "BasePDE::Mesh2PDENode " );

  PDENodes.Resize(MeshNodes.GetSize());
  
  for (Integer i=0; i<MeshNodes.GetSize(); i++) 
     PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];
}




void BasePDE::PDE2MeshNode(Vector<Integer> & meshNodes, 
			   const Vector<Integer> & pdeNodes,
			   const std::vector<Integer> & pde2MeshNode)
{

  ENTER_FCN( "BasePDE::PDE2MeshNode " );

  meshNodes.Resize(pdeNodes.GetSize());

  for (Integer i=0; i<pdeNodes.GetSize(); i++)
    meshNodes[i] = pde2MeshNode[pdeNodes[i]-1];

#ifdef DEBUG
  (*debug) << "--------------------" << std::endl;
  (*debug) << " PDE2MeshNode()" << std::endl;
  for (Integer i=0; i<pdeNodes.GetSize(); i++)
    (*debug) << "in: " << pdeNodes[i] << " out: " << meshNodes[i] << std::endl;
#endif
}




void BasePDE::AssignPDENodeNumbers(std::vector<Integer> & mesh2PDENode,
			  std::vector<Integer> & pde2MeshNode,
			  const std::vector<std::string> & subdoms)
{

  ENTER_FCN( "BasePDE:AssignPDENodeNumbers:" );

  // Initialize Mesh2PDENode and PDE2MeshNode
  mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  std::vector<Elem*> SD;
  Integer nodecounter = 1;

#ifdef ADAPTGRID
  std::cout << "NO MAPPING OF NODES!! " << std::endl << std::endl;
  
  PDEMeshNode_.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  for (Integer i=0;i<ptgrid_->GetMaxnumnodes(actlevel_);i++)
    {
      mesh2PDENode_[i] = i+1;
      pde2MeshNode_[i] = i+1;
    }
  numPDENodes_ = pde2MeshNode_.size();
  
  return;
#endif  

  // Iterate over Subdomains
  for (Integer i=0; i<subdoms.size(); i++)
    {
      ptgrid_->GetElemSD(SD,subdoms[i],actlevel_);
      // Iterate over all elements in subdomain
      for (Integer j=0; j<SD.size(); j++)
	{
	  // Iterate over all element nodes
	  for (Integer numNodes=0; numNodes<SD[j]->connect.GetSize(); numNodes++)
	    {
	      // Check if node was already assigned
	      if (mesh2PDENode[SD[j]->connect[numNodes] - 1] == -1)
		{
		  mesh2PDENode[SD[j]->connect[numNodes] - 1] = nodecounter;
		  pde2MeshNode.push_back(SD[j]->connect[numNodes]);
		  nodecounter++;
		}
	    }
	}
    }

  numPDENodes_ = pde2MeshNode_.size();
  
}

void BasePDE::AssignPDEElemNumbers(std::vector<Integer> & mesh2PDEElem,
				   std::vector<Integer> & pde2MeshElem,
				   const std::vector<std::string> & subdoms)
{

  ENTER_FCN( "BasePDE::AssignPDEElemNumbers" );

  // Resize mesh2PDEElem and pde2MeshElem to correct size
  // and initialize with -1 = not assigned to PDE
  mesh2PDEElem.resize(ptgrid_->GetMaxnumElem(actlevel_),-1);
  pde2MeshElem.resize(ptgrid_->GetMaxnumElem(actlevel_,subdoms),-1);

  //std::cerr << "mesh2PDEElem.size = " << mesh2PDEElem.size() << std::endl;
  //std::cerr << "PDE2meshElem.size = " << pde2MeshElem.size() << std::endl;
  std::vector<Elem*> SD;
  Integer elemCounter = 1;


  // Iterate over Subdomains
  for (Integer iSD=0; iSD<subdoms.size(); iSD++)
    {
      ptgrid_->GetElemSD(SD,subdoms[iSD],actlevel_);

      // Iterate over all elements in subdomain
      for (Integer iElem=0; iElem<SD.size(); iElem++)
	{
	  //std::cerr << "assigning mesh2PDEElem[SD[iElem]->elemNum - 1 ] = elemCounter: ";
	  //std::cerr << SD[iElem]->elemNum - 1 << "=" << elemCounter << std::endl;
	  //std::cerr << "assigning pde2MeshElem[elemCounter-1] = SD[iElem]->elemNum: ";
	  //std::cerr << elemCounter-1 << " = " << SD[iElem]->elemNum << std::endl;
	  mesh2PDEElem[SD[iElem]->elemNum - 1 ] = elemCounter;
	  pde2MeshElem[elemCounter-1] = SD[iElem]->elemNum;
	  elemCounter++;
	}
      
    }
  
  numElems_ = mesh2PDEElem.size();
  
				  
}
     

void BasePDE::GetElemCoords(const Vector<Integer> connect, Matrix<Double> &coordMat, const Integer level)
{

  ENTER_FCN( "BasePDE::GetElemCoords" );

  ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
  if (deltCoords_.GetSizeRow() != 0 && geoUpdate_ == TRUE)
    {
      for (Integer i=0; i<coordMat.GetSizeRow(); i++)
	for (Integer j=0; j<coordMat.GetSizeCol(); j++) 
	  coordMat(i,j) += deltCoords_(i,mesh2PDENode_ [connect[j] - 1] - 1);
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
  std::vector<Elem*>       elemssd;
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

void BasePDE::GetDerivSolVecOfElement(Vector<Double>& sol, Vector<Integer>& connect_PDE)
{

  ENTER_FCN( "BasePDE::GetDerivSolVecOfElement" );

  // displacement of element nodes
  sol.Resize(dofspernode_ * connect_PDE.GetSize());
  
  const Vector<Double> & sol_der1 = getS1();

  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof + actNode*dofspernode_] = sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];
}





void BasePDE::GetDerivSolOfElement(Matrix<Double>& sol, Vector<Integer>& connect_PDE)
{

  ENTER_FCN( "BasePDE::GetDerivSolOfElement" );

  // displacement of element nodes
  sol.Resize(dofspernode_, connect_PDE.GetSize());

  const Vector<Double>& sol_der1 = getS1();  

  for(Integer actNode=0; actNode<connect_PDE.GetSize(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof][actNode] =  sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];

}







void BasePDE::CalcLineNormalVec(Vector<Double>& n, const Elem& interfaceElem, const Elem& neighbour)
{

ENTER_FCN( "BasePDE::CalcLineNormalVec" );

  const Integer nrVecElem2d = 2;

  // get coords of interface elemmt ========================================
  BaseFE * ptElem = interfaceElem.ptElem;
  Vector<Integer> connecth = interfaceElem.connect;
  
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




} // end of namespace

