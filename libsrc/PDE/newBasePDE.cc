#ifdef NEWBASEPDE

#include <fstream>
#include <iostream>
#include <string>

#include <DataInOut/conffile.hh> 
#include <Estimator/spaceerror.hh>
#include <DataInOut/WriteInfo.hh>

#include "newBasePDE.hh"


namespace CoupledField
{

BasePDE::BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile, WriteResults * aOutFile, 
		 TimeFunc *aptTimeFunc)
  :nonLin_(FALSE),
   incStopCrit_(1e-2), 
   residualStopCrit_(1e-3),
   firstTimeStepStatic_(TRUE),
   isaxi_(FALSE)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

  nonLin_ = FALSE;
  incStopCrit_ = 1e-2;
  residualStopCrit_ = 1e-3;

  InFile_     = aInFile;
  OutFile_    = aOutFile;
  ptTimeFunc_ = aptTimeFunc;
  ptgrid_     = aptgrid;
  ptBCs_      = aptBCs;

  actlevel_ = 0;
  couplingBCsCounter_ = 0;
  numDirichletBCs_ = 0;
  PDEisCoupled_ = FALSE;
  updateCouplingBCs_ = FALSE;
  updateBCs_ = 0;
  Dim_ = ptgrid_->GetDim();
  InitMatrices_ = FALSE;

  //standard parameter for solver
  eps_         = 1.0e-8;
  dampiter_    = 1.0;
  maxnumiter_  = 100;
    solvertype_  = RealDirect;
    //solvertype_  = RealCG;
  precondtype_ = ID;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;

  savederiv1_ = FALSE;
  savederiv2_ = FALSE;
  

   //get analysis type
  std::string analysis;
  conf->get("analysis", analysis);

  //allocate according algebraic system
  algsys_ = new StandardSystem();


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


  //for adaptivity

  if (conf->get_option("adaptspace"))
    conf->get("tolerance_space_error", tolSpaceErr_);
  else
    tolSpaceErr_ = .0;
}


void BasePDE::ReadSavings()
{
#ifdef TRACE
  (*trace) << "entering BasePDE::ReadSavings" << std::endl;
#endif

  //set saving of solution to yes, if user has not used the nodalsave-command
  savesol_ = TRUE;

  //check for node saving
  std::vector<std::string> savings;
  

  //reset saving of solution, if user has used the nodalsave-command
  if (conf->ifgetliststr("nodalsave", savings, pdename_))
    savesol_ = FALSE;

  for (Integer i=0; i<savings.size(); i++)
    {
      if (savings[i] == "dof")  savesol_ = TRUE;
      if (savings[i] == "deriv1") savederiv1_ = TRUE;
      if (savings[i] == "deriv2") savederiv2_ = TRUE;
    }
}


void BasePDE::WriteGeneralPDEdefines()
{
#ifdef TRACE
  (*trace) << "entering BasePDE::WriteGeneralPDEdefines" << std::endl;
#endif

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

      Info->WriteInHomBC(pdename_, bcs_id_[i], val_id_[i], fncnames_id_[i], dof);	
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
#ifdef TRACE
  (*trace) << "entering BasePDE::SolveStepStatic" << std::endl;
#endif

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
#ifdef TRACE
  (*trace) << "entering BasePDE::StepStaticLin" << std::endl;
#endif

  Integer job;
  Double * ptsol;
  lasttimecalc_ = aTime; // for correct output in unv-file
  

  if ( PDEisCoupled_ == FALSE || iterCoupledCounter_ == 0 && firstTimeStepStatic_)
    assemble_->AssembleMatrices(level);


  //this has to be done each time!
  assemble_->AssembleSrcRHS(level, aTime);

  if ( PDEisCoupled_ == FALSE || iterCoupledCounter_ == 0)
    {
      //account for bcs
      SetBCs(level, updateBCs_, aTime);
      updateBCs_ = 0;

      if (firstTimeStepStatic_)
	job = 1; // calc new preconditioner
      else 
	job = 3;  // just set new BCs    
    }
  else
    job = 3;

  algsys_->CalcPrecond(job);
  

  algsys_->Solve();

  ptsol = algsys_->GetSolutionVal();

  // save solution
  Integer k=0;
  
  for (Integer i=0; i<numPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptsol[k++];

  firstTimeStepStatic_ = FALSE;
}


void BasePDE::SolveStepTrans(const Integer kstep, const Double asteptime, 
			     const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::SolveStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;
  Recalc_ = FALSE;

  if (laststepcalc_==kstep && kstep!=0) 
    Recalc_=TRUE;
  else 
    laststepcalc_= kstep;

  if (nonLin_)
    StepTransNonLin(kstep,asteptime, level,reset);
  else
    StepTransLin(kstep,asteptime,level,reset);
}


void BasePDE::PreStepTrans(const Integer kstep, const Double asteptime,
			   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::PreStepTrans" << std::endl;
#endif

  lasttimecalc_= asteptime;

  // due to coupling-pdes, the RHS has to be initialized BEFORE 
  // the coupling forces are assembled to the RHS

  algsys_->InitRHS();
  assemble_->AssembleSrcRHS(level,lasttimecalc_);
  
}



void BasePDE::PostStepTrans(const Integer kstep, const Double asteptime, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::PostStepTrans" << std::endl;
#endif

  if (PDEisCoupled_)
    {
      //save solution
      Integer k=0;
      Array<Double> solhelp(1, numPDENodes_*dofspernode_);

      for (Integer i=0; i< numPDENodes_; i++)
	for (Integer dim=0; dim<dofspernode_; dim++)
	  {
	    solhelp[0][k] = sol_[dim][i];
	    k++;
	  }
    
      TS_alg_->Corrector(solhelp);  //perform corrector step
    }
  
  if (PDEisCoupled_)
    iterCoupledCounter_++;
}


// initialize PDEs before iteration (done for each time step)
void BasePDE::InitStepTransCoupled(Double asteptime)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::InitStepTransCoupled" << std::endl;
#endif

  lasttimecalc_= asteptime;
  iterCoupledCounter_ = 0;
}



void BasePDE::StepTransLin(const Integer kstep, const Double asteptime,
			   const Integer level, const Boolean reset)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::StepTransLin" << std::endl;
#endif

  Double * ptsol;
  Integer update,job;

  //current problem with Array class
  Array<Double> solhelp;
  solhelp.reshape(1, numPDENodes_*dofspernode_);
  Integer k=0;


  if ( PDEisCoupled_ == FALSE || iterCoupledCounter_ == 0)
    {    
      for (Integer i=0; i< numPDENodes_; i++)
	for (Integer j=0; j<dofspernode_; j++)
	  {
	    solhelp[0][k] = sol_[j][i];
	    k++;
	  }
      
      //perform predictor step
      TS_alg_->Predictor(solhelp);
    }
  

  if (laststepcalc_==0)
    {
      update = 0;
      job = 3;

      if ( PDEisCoupled_ == FALSE || iterCoupledCounter_ == 0)
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
  algsys_->CalcPrecond(job);
  algsys_->Solve();
  ptsol = algsys_->GetSolutionVal();

  //save solution
  k=0;
  for (Integer i=0; i< numPDENodes_; i++)
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i]  = ptsol[k];
	solhelp[0][k] = ptsol[k];
	k++;
      }
    
  if (!PDEisCoupled_)
    TS_alg_->Corrector(solhelp);  //perform corrector step
}


void  BasePDE::SetBCs(const Integer level, const Integer update, const Double time)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetBCs" << std::endl;
#endif

  Integer node, dof;
  Double val, val_tfunc;

  std::list<Integer> nodes;

  Integer i;
  Integer j;
  if (PDEisCoupled_)
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
	  algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1],val, dof, SYSTEM);
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
	  algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1],
				val, dof, SYSTEM);
	}
    }
  
#ifdef TRACE
  (*trace) << " leaving BasePDE::SetBCs " << std::endl;
#endif 
}




void BasePDE::ReadBCs(const std::string eq)
{
#ifdef TRACE
  (*trace) << " entering BasePDE::ReadBCs " << std::endl;
#endif


  conf->ifgetliststr("homogeneous_dirichlet",bcs_hd_,eq); 
  conf->ifgetliststr("inhomogeneous_dirichlet",bcs_id_,eq);

  Integer i;

  val_id_.resize(bcs_id_.size());
  fncnames_id_.resize(bcs_id_.size());

  for(i=0; i<bcs_id_.size(); i++)
    conf->get2(bcs_id_[i], val_id_[i], fncnames_id_[i], eq,"bc_conditions","inhomogeneous_dirichlet");
}



void BasePDE::ReadMaterialData()
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::ReadMaterialData" << std::endl;
#endif

  // read material-file name from config-file
  std::string matFileName;
  conf->get("material_file", matFileName);
  charMaterialFileName_ = c_string(matFileName);
  loadMaterial_ = new LoadMaterialData(charMaterialFileName_);

  //read material data for each subdomain
  materialData_  = new MaterialData[subdoms_.size()];

  std::string matName;
  for (Integer i=0; i<subdoms_.size(); i++)
    {
      // load material data into array "materialData_"
      conf->getstr(subdoms_[i], matName, "list_subdomains");
      loadMaterial_->GetMaterial(materialData_[i], matName, pdematerialclass_);
    }

  delete loadMaterial_;
  
#ifdef TRACE
  (*trace) << "leaving BasePDE::ReadMaterialData" << std::endl;
#endif

}


Integer BasePDE::GetNumRestraints(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::GetNumRestraints" << std::endl;
#endif
    
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






BasePDE::~BasePDE()
{
#ifdef TRACE
  (*trace) << " entering BasePDE::~BasePDE() " << std::endl;
#endif

  // ATTENTION: Dummy value for as_id!!!!!!!!!!!!!!!!!!!!!!!!!!
  DeleteAlgSys(0);
}




  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================

  void BasePDE::SetAlgSys(int sysid)
  {
#ifdef TRACE
    (*trace) << "entering  Analysis::SetAlgSys" << std::endl;
#endif

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
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetSolverParameters" << std::endl;
#endif

  //if values are defined in conf-file, take these
  conf->ifget("eps",eps_,pdename_); // relative accuracy in the precond. energy
  conf->ifget("dampiter",dampiter_,pdename_); // damping parameter for Jacobi, SSOR
  conf->ifget("maxnumit",maxnumiter_,pdename_); // maximal number of iterations
  conf->ifget("solvertype",solvertype_,pdename_); // solver
  conf->ifget("precondtype", precondtype_,pdename_); //preconditioner
  conf->ifget("numeqcoarse",numeqcoarse_,pdename_); // number of equation for coarsing
  conf->ifget("coarsealpha",coarsealpha_,pdename_); // coarsing parameter for AMG
  
  if (solvertype_==RealDirect && precondtype_!=ID)  precondtype_=ID;

  //communicate with algebraic system
  algsys_->CreateParameter();
  algsys_->SetAccuracy(eps_);
  algsys_->SetMaxNumIter(maxnumiter_);
  algsys_->SetPrecond(precondtype_);
  algsys_->SetSolver(solvertype_);
  algsys_->SetDampIter(dampiter_);
  algsys_->SetCoarseSystem(numeqcoarse_);
  algsys_->SetAlpha(coarsealpha_);

} 




void BasePDE::CreateMatrices_Solver()
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::CreateMatrices_Solver" << std::endl;
#endif
  // create and initialize matrices 
  assemble_->CreateMatrices();  
  assemble_->InitMatrices();  

  //create solver and preconditioner
  algsys_->CreateSolver();
  algsys_->CreatePrecond(assemble_->GetMatrixType());

  //now reset AlgebraicSystem 
  algsys_->InitRHS();
  algsys_->InitSol();
}





void BasePDE::StoreToSolArray(Double * ptSol)
{
  Integer k=0;

  for (Integer i=0; i<numPDENodes_; i++)   
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptSol[k++];
}



void BasePDE::StoreVecToSolArray(std::vector<Double>& sol)
{
  Integer k=0;

  for (Integer i=0; i<numPDENodes_; i++)   
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i] = sol[k];
	k++;
      }
}



  // ======================================================
  // COUPLING SECTION
  // ======================================================


void BasePDE::CalcInputCoupling()
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::CalcInputCoupling" << std::endl;
#endif

  std::vector<Integer> * nodes;
  std::vector<Elem*> * elements;
  Array<Double> * val;
  Integer PDEnode;
  
  // Reset counter for boundary conditions
  couplingBCsCounter_ = 0;
  
  for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++)
    {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);

      switch(ptCoupling_->GetInputType(i))
	{
	  
	case COORD:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(Coord)" << std::endl;
	  InitMatrices_ = TRUE;
	  ptCoupling_->GetInputNodes(i, nodes);
	  deltCoords_.reshape(Dim_, numPDENodes_);

	  // set ptr of deltCoords to assembly-object
	  assemble_->SetPtrDeltaCoordinates(&deltCoords_);

	  for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	    for (Integer j=0; j<nodes->size(); j++)
	      {
		//std::cerr << "processing dim = " << dim << ", j = " << j << std::endl;
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1]-1;
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);

		deltCoords_[dof][PDEnode] = (*val)[dof][j];
	      }
	      
	  break;

	case RHS:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(RHS)" << std::endl;
	  ptCoupling_->GetInputNodes(i, nodes);

	  for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	    for (Integer j=0; j<nodes->size(); j++)
	      {
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  {
		    std::cerr << "PDENODE: "  << PDEnode << "Node[" << (*nodes)[j] << "][" 
			      << dof+1 << "]= " << (*val)[dof][j] << std::endl; 
		    Error("Node not assigned to coupling domain: see mesh- and config-file"
			  , __FILE__,__LINE__);
		  }
		algsys_->SetNodeRHS((*val)[dof][j], PDEnode, dof+1);
	      }
	  
	  break;

	case ID_BC:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(ID_BC)" << std::endl;
 	  ptCoupling_->GetInputNodes(i, nodes);
	  //std::cerr << "found " << nodes->size() << " IDBCs input nodes" << std::endl;
	  for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
	    for (Integer j=0; j<nodes->size(); j++, couplingBCsCounter_++)
	      {
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);
		
		if (updateCouplingBCs_)
		  {
		    //		    std::cerr << "updating BC[" << dof << "][" << (*nodes)[j] << "] = " 
		    //		      << (*val)[dof][j] << std::endl;
		    algsys_->UpdateDirichlet(couplingBCsCounter_+1, (*val)[dof][j], SYSTEM);
		  }
		else
		  {	
		    //  std::cerr << "BC[" << dim << "][" << (*nodes)[j] << "] = " << (*val)[dim][j] << std::endl;
		    algsys_->SetDirichlet(couplingBCsCounter_+1, PDEnode, (*val)[dof][j], dof+1, SYSTEM);
		  }
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
#ifdef TRACE
  (*trace) << "entering BasePDE::Mesh2PDENode " << std::endl;
#endif

  PDENodes.Resize(MeshNodes.size());
  
  for (Integer i=0; i<MeshNodes.size(); i++) 
     PDENodes[i] = Mesh2PDENode[MeshNodes[i]-1];
}




void BasePDE::PDE2MeshNode(Vector<Integer> & MeshNodes, 
			   const Vector<Integer> & PDENodes,
			   const std::vector<Integer> & PDE2MeshNode)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::PDE2MeshNode " << std::endl;
#endif

  MeshNodes.Resize(PDENodes.size());

  for (Integer i=0; i<PDENodes.size(); i++)
    MeshNodes[i] = PDE2MeshNode[PDENodes[i]-1];

#ifdef DEBUG
  (*debug) << "--------------------" << std::endl;
  (*debug) << " PDE2MeshNode()" << std::endl;
  for (Integer i=0; i<PDENodes.size(); i++)
    (*debug) << "in: " << PDENodes[i] << " out: " << MeshNodes[i] << std::endl;
#endif
}




void BasePDE::AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			  std::vector<Integer> & PDE2MeshNode,
			  const std::vector<std::string> & subdoms)
{
#ifdef TRACE
  (*trace) << "entering BasePDE:AssignPDENodeNumbers:" << std::endl;
#endif

  // Initialize Mesh2PDENode and PDE2MeshNode
  Mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  std::vector<Elem*> SD;
  Integer NodeCounter = 1;

#ifdef ADAPTGRID
  std::cout << "NO MAPPING OF NODES!! " << std::endl << std::endl;
  
  PDE2MeshNode_.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  for (Integer i=0;i<ptgrid_->GetMaxnumnodes(actlevel_);i++)
    {
      Mesh2PDENode_[i] = i+1;
      PDE2MeshNode_[i] = i+1;
    }
  numPDENodes_ = PDE2MeshNode_.size();
  
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
	  for (Integer NumNodes=0; NumNodes<SD[j]->connect.size(); NumNodes++)
	    {
	      // Check if node was already assigned
	      if (Mesh2PDENode[SD[j]->connect[NumNodes] - 1] == -1)
		{
		  Mesh2PDENode[SD[j]->connect[NumNodes] - 1] = NodeCounter;
		  PDE2MeshNode.push_back(SD[j]->connect[NumNodes]);
		  NodeCounter++;
		}
	    }
	}
    }

  numPDENodes_ = PDE2MeshNode_.size();

#ifdef DEBUG
  (*debug) << "Mesh2PDENodes:" << std::endl << Mesh2PDENode << std::endl;
#endif
}



  

void BasePDE::GetElemCoords(const Vector<Integer> connect, Matrix<Double> &coordMat, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE:GetElemCoords:" << std::endl;
#endif

  ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
  if (deltCoords_.size() != 0 && GeoUpdate_ == TRUE)
    {
      for (Integer i=0; i<coordMat.size_row(); i++)
	for (Integer j=0; j<coordMat.size_col(); j++) 
	  coordMat(i,j) += deltCoords_[i][Mesh2PDENode_ [connect[j] - 1] - 1];
    }

}

void BasePDE::TransformNodeSolution(Array<Double> & MeshSol, 
				    Array<Double> & PDESol, 
				    const std::vector<Integer> & PDE2MeshNode)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformNodeSolution" << std::endl;
#endif

  Integer node, idx;
  MeshSol.reshape(PDESol.dim(), ptgrid_->GetMaxnumnodes(actlevel_));

  // loop over dimensions
  for (Integer dim=0; dim<PDESol.dim(); dim++)

      // Loop over all PDE nodes
      for (Integer i=0; i<PDE2MeshNode.size(); i++)

	  MeshSol[dim][PDE2MeshNode[i]-1] = PDESol[dim][i];

}


void BasePDE::TransformElemSolution(Array<Double> & MeshSol, 
				    Array<Double> & PDESol, 
				    const std::vector<Elem*> & Elems)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformElemSolution Elem*" << std::endl;
#endif

   MeshSol.reshape(PDESol.dim(), ptgrid_->GetMaxnumElem(actlevel_));
   std::cout << "dim= " << ptgrid_->GetMaxnumElem(actlevel_) << std::endl;
   
  // loop over all dimensions
  for (Integer dim=0; dim<PDESol.dim(); dim++)

    // loop over all elements
    for (Integer i=0; i<Elems.size(); i++) {
      MeshSol[dim][i] = PDESol[dim][i];
    }
  
}

void BasePDE::TransformElemSolution(Array<Double> & MeshSol, 
				    Array<Double> & PDESol, 
				    const std::vector<std::string> & SD)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformElemSolution string" << std::endl;
#endif

  MeshSol.reshape(PDESol.dim(), ptgrid_->GetMaxnumElem(actlevel_,SD));


  Integer elMesh=0;
  Integer elPDE=0;
  std::vector<std::string> AllSDs = ptgrid_->GetListSubDomains();

  // loop over all SubDomains of computational domain
  for (Integer isd=0; isd<AllSDs.size(); isd++)
  {
    Boolean SDbelongsToDomain = FALSE;
    for (Integer k=0; k<SD.size(); k++)
      if (SD[k] == AllSDs[isd]) 
	SDbelongsToDomain = TRUE;

    std::vector<Elem*> Elems;
    ptgrid_->GetElemSD(Elems, AllSDs[isd], actlevel_);    
    if (SDbelongsToDomain)
      {
	//computational subdomain belongs to PDE 
	// loop over all elements
	for (Integer i=0; i<Elems.size(); i++) 
	    {
	      // loop over dim
	      for (Integer dim=0; dim<PDESol.dim(); dim++)
		MeshSol[dim][elMesh] = PDESol[dim][elPDE]; 

	      elPDE++; elMesh++;
	    }
      }
    else
      elMesh += Elems.size();
    
  }
}


void BasePDE::NodeSolutionToCoupling(Array<Double>& CouplingSol,
				     const std::vector<Integer>& NodeNumbers)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::NodeSolutionToCoupling" << std::endl;
#endif
  
  CouplingSol.reshape(dofspernode_, NodeNumbers.size());
  
  for (Integer i=0; i<CouplingSol.dim(); i++)
    for (Integer j=0; j<CouplingSol.size(); j++)
      {
	//std::cerr << "processing dim: " << i <<", j:" << j << std::endl; 
	CouplingSol[i][j] = sol_[i][Mesh2PDENode_[NodeNumbers[j]-1 ] - 1];
      }
}





void BasePDE::ElemSolutionToCoupling(Array<Double>& CouplingSol,
				     const std::vector<Elem*>& NodeNumbers,
				     Vector<Double>& elemSol)
				     
{
#ifdef TRACE
  (*trace) << "entering BasePDE::ElemSolutionToCoupling" << std::endl;
#endif
  
  const Integer couplingDof = CouplingSol.dim();
  
  for (Integer actDof=0; actDof < couplingDof; actDof++)
    for (Integer actNode=0; actNode < NodeNumbers.size(); actNode++)
      {
	//std::cerr << "processing dim: " << i <<", j:" << j << std::endl; 
	CouplingSol[actDof][actNode] = elemSol[actDof + actNode*couplingDof];
      }
}






  // ======================================================
  // ADAPTIVITY SECTION 
  // ======================================================


#ifdef ADAPTGRID
void BasePDE::RefineMesh(const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::RefineMesh" << std::endl;
#endif
  
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
#ifdef TRACE
  (*trace) << "entering BasePDE::TestError" << std::endl;
#endif

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
#ifdef TRACE
  (*trace) << "entering BasePDE::ConstructorError" << std::endl;
#endif

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




void BasePDE::GetSolVecOfElement(Vector<Double>& sol, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::GetSolVecOfElement" << std::endl;
#endif

  // displacement of element nodes
  sol.Resize(dofspernode_ * connect_PDE.size());
  
  for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof + actNode*dofspernode_] = sol_[actDof][connect_PDE[actNode]-1];
}



void BasePDE::GetDerivSolVecOfElement(Vector<Double>& sol, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::GetDerivSolVecOfElement" << std::endl;
#endif

  // displacement of element nodes
  sol.Resize(dofspernode_ * connect_PDE.size());

  const Array<Double>& sol_der1Array = getS1();
  
  for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof + actNode*dofspernode_] = sol_der1Array[0][actDof + dofspernode_*(connect_PDE[actNode]-1)];
}





void BasePDE::GetDerivSolOfElement(Matrix<Double>& sol, Vector<Integer>& connect_PDE)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::GetDerivSolOfElement" << std::endl;
#endif

  // displacement of element nodes
  sol.Resize(dofspernode_, connect_PDE.size());

  const Array<Double>& sol_der1Array = getS1();  

  for(Integer actNode=0; actNode<connect_PDE.size(); actNode++)
    for(Integer actDof=0; actDof < dofspernode_; actDof++)
      sol[actDof][actNode] =  sol_der1Array[0][actDof + dofspernode_*(connect_PDE[actNode]-1)];
  //sol_der1Array[actDof][connect_PDE[actNode]-1];
}







void BasePDE::CalcLineNormalVec(Vector<Double>& n, const Elem& interfaceElem, const Elem& neighbour)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::CalcLineNormalVec" << std::endl;
#endif

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
  
  for(Integer actNode=0; actNode < connecth.size(); actNode++)
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
      (indexNode2-indexNode1)+connecth.size() == 1 )
    n *= -1;
  
  else
    // if not clockwise orientation of nodes (difference of node indizes is -1)
    if (! (indexNode2-indexNode1 == -1 || 
	   (indexNode2-indexNode1)-connecth.size()==-1) )
      Error("Nodes of interface don't lie beneath each other in neighbouring element!", __FILE__, __LINE__);  
}





// normal of line element: ATTENTION no defined sign!!
void BasePDE::CalcLineNormalVec(Vector<Double>& n, Matrix<Double>& ptCoord)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::CalcLineNormalVec" << std::endl;
#endif

  const Integer nrVecElem2d = 2;
  
  if (ptCoord.size_row()!=nrVecElem2d)
    Error("Calc element normal: no line element! ", __FILE__,__LINE__);

  n.Resize(nrVecElem2d);
  
  // normal of a vector: interchange x-coord and y-coord and take the new x-coord as negative

  n[0] = -(ptCoord[1][1] - ptCoord[1][0]);
  n[1] =  (ptCoord[0][1] - ptCoord[0][0]);
  n /= n.normL2();
}




} // end of namespace

#endif // #ifdef NEWBASEPDE
