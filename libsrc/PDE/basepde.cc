#ifndef NEWBASEPDE

#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
#include <DataInOut/conffile.hh>
 
#include <Estimator/spaceerror.hh>

namespace CoupledField
{

BasePDE::BasePDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile, WriteResults * aOutFile, 
		 TimeFunc *aptTimeFunc)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

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
  precondtype_ = ID;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;

  //axisymmetric case
  isaxi_ = FALSE;

  //get analysis type
  std::string analysis;
  conf->get("analysis", analysis);

  if (analysis=="static") 
    analysistype_ = STATIC;
  else if (analysis=="transient")
    analysistype_ = TRANSIENT;
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


void BasePDE::SetAlgSys(const Integer as_sysid)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetAlgSys" << std::endl;
#endif

  as_sysid_ = as_sysid;

  //allocate according algebraic system
  algsys_ = new StandardSystem();

  //set solver parameters  
  SetSolverParameters();

  //ask the PDE discrete form
  DiscreteParamsPDE();

  //set the graph type used for the system matrices
  SetupMatrixGraph(NumPDENodes_, GraphType_);

  //allocate the necessary matrices as well as solver and preconditioner
  CreateMatrices_Solver();
}

void BasePDE::DeleteAlgSys(const Integer as_sysid)
{
  if (algsys_)  delete algsys_;
}

void BasePDE::InitMatrices()
{
 //  //Initialize matrices in order to get BCs correct
  std::vector<Integer> matrixsystype(5,0);    
  if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
  if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
  if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
  if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix
  

  for (Integer i=0;i<5;i++)
    if (matrixsystype[i] !=0)
      algsys_->InitMatrix(i+1);
 
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
#ifdef TRACE
  (*trace) << "leaving BasePDE::ReadMaterialData" << std::endl;
#endif

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

void BasePDE::SetupMatrixGraph(Integer numeq, Integer graphtype)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetupMatrixGraph" << std::endl;
#endif

  //initialize matrix graph
  algsys_->CreateGraph(numeq,graphtype);

  // set the graph - connectivity matrix
  BaseFE * ptElem; 
  Integer nsub, iel, fe_type;
  Vector<Integer> connecth;

  for (nsub=0; nsub<subdoms_.size(); nsub++)
    {
      std::vector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[nsub],actlevel_);

      for (iel=0; iel < elemssd.size(); iel++)
	{  
	  ptElem=elemssd[iel]->ptElem;
	  
	  //Map Mesh Node numbers to PDE node numbers
	  Mesh2PDENode(connecth,elemssd[iel]->connect,Mesh2PDENode_);

	  fe_type=elemssd[iel]->ptElem->feType();
	  algsys_->SetElementPos(connecth.get(),connecth.size(),fe_type);
	}
    }

}


void BasePDE::CreateMatrices_Solver()
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::CreateMatrices_Solver" << std::endl;
#endif

  Integer matrixclass;
  Integer matrixsystype[5];    
  Integer graphtype; 
  Integer numdofpernode;
  Integer numconstraints;

  //ask the PDE discrete form
  // DiscreteParamsPDE();

  numdofpernode  = dofspernode_; 

  if (PDEisCoupled_)
    numDirichletBCs_  += GetNumRestraints(actlevel_);
  else
    numDirichletBCs_  = GetNumRestraints(actlevel_);

#ifdef DEBUG
  (*debug) << "Num Dirichlet nodes: " <<  numDirichletBCs_ << std::endl;
#endif

  numconstraints = 0;  // currently not handled
  matrixclass    = MatrixType_;
  graphtype      = GraphType_;
  
  if (SystemMatrix_    != 1)
    Error("One needs at least a system matrix!",__FILE__,__LINE__);

  if (SystemMatrix_     == 1) matrixsystype[0] = SYSTEM;      // memory for the system matrix
  if (StiffnessMatrix_  == 1) matrixsystype[1] = STIFFNESS;   // memory for the stiffness matrix
  if (DampingMatrix_    == 1) matrixsystype[2] = DAMPING;     // memory for the damping matrix
  if (ConvectionMatrix_ == 1) matrixsystype[3] = CONVECTION;  // memory for the convection matrix
  if (MassMatrix_       == 1) matrixsystype[4] = MASS;        // memory for the mass matrix

  //put to algebraic system
  algsys_->CreateMatrix(matrixsystype, matrixclass, graphtype, numdofpernode, numDirichletBCs_,
			numconstraints);

  //create solver and preconditioner
  algsys_->CreateSolver();
  algsys_->CreatePrecond(matrixclass);

  //now reset AlgebraicSystem 
  algsys_->InitRHS();
  algsys_->InitSol();

  InitMatrices();
  
}

void BasePDE::StoreToSolArray(Double * ptSol)
{
  Integer k=0;

  for (Integer i=0; i<NumPDENodes_; i++)   
    for (Integer dim=0; dim<dofspernode_; dim++)
      sol_[dim][i] = ptSol[k++];
}


void BasePDE::StoreVecToSolArray(std::vector<Double>& sol)
{
  Integer k=0;

  for (Integer i=0; i<NumPDENodes_; i++)   
    for (Integer dim=0; dim<dofspernode_; dim++)
      {
	sol_[dim][i] = sol[k];
	k++;
      }
}



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

      ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);

      switch(ptCoupling_->GetInputType(i))
	{
	  
	case COORD:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(Coord)" << std::endl;
	  InitMatrices_ = TRUE;
	  ptCoupling_->GetInputNodes(i, nodes);
	  deltCoords_.reshape(Dim_, NumPDENodes_);

	  for (Integer dim=0; dim<ptCoupling_->GetInputDim(i); dim++)
	    for (Integer j=0; j<nodes->size(); j++)
	      {
		//std::cerr << "processing dim = " << dim << ", j = " << j << std::endl;
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1]-1;
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);

		deltCoords_[dim][PDEnode] = (*val)[dim][j];
	      }
	      
	  break;

	case RHS:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(RHS)" << std::endl;
	  ptCoupling_->GetInputNodes(i, nodes);

	  for (Integer dim=0; dim<ptCoupling_->GetInputDim(i); dim++)
	    for (Integer j=0; j<nodes->size(); j++)
	      {
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);

		//	std::cerr << "PDENODE: "  << PDEnode << "Node[" << (*nodes)[j] << "][" << dim+1 << "]= " << (*val)[dim][j] << std::endl; 
		algsys_->SetNodeRHS((*val)[dim][j], PDEnode, dim+1);
	      }
	  
	  break;

	case ID_BC:
	  //std::cerr << "In " << pdename_ << "::CalcInputCoupling - Switch(ID_BC)" << std::endl;
 	  ptCoupling_->GetInputNodes(i, nodes);
	  //std::cerr << "found " << nodes->size() << " IDBCs input nodes" << std::endl;
	  for (Integer dim=0; dim<ptCoupling_->GetInputDim(i); dim++)
	    for (Integer j=0; j<nodes->size(); j++, couplingBCsCounter_++)
	      {
		PDEnode = Mesh2PDENode_[(*nodes)[j]-1];
		if (PDEnode==-1)
		  Error("Node not assigned to coupling domain: see mesh- and config-file",__FILE__,__LINE__);
		
		if (updateCouplingBCs_)
		  {
		    //std::cerr << "updating BC[" << dim << "][" << (*nodes)[j] << "] = " << (*val)[dim][j] << std::endl;
		    algsys_->UpdateDirichlet(couplingBCsCounter_+1, (*val)[dim][j], SYSTEM);
		  }
		else
		  {	
		    //std::cerr << "BC[" << dim << "][" << (*nodes)[j] << "] = " << (*val)[dim][j] << std::endl;
		    algsys_->SetDirichlet(couplingBCsCounter_+1, PDEnode, (*val)[dim][j], dim+1, SYSTEM);
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


void  BasePDE::SetBCs(const Integer level, const Integer update, const Double time)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetBCs" << std::endl;
#endif

//   Double BigConst;
//   if (SolverCFS_)
//     {
//       Integer matsize=sysmat_.getSize(),i;
//       Double max=sysmat_(0,0);
//       for (i=0;i<matsize;i++)
// 	if (sysmat_(i,i)>max) max=sysmat_(i,i);
//       BigConst=(10e+10)*max;
//     }
  
  Integer node;
  Double val, val_tfunc;

  val_tfunc = 1.0;
  if (ptTimeFunc_->GetmaxTimeFnc()!=0 && analysistype_!=HARMONIC)
      val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,level);

  std::list<Integer> nodes;

  Integer i;
  Integer j;
  if (PDEisCoupled_)
    j = couplingBCsCounter_;
  else
    j=0;

  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  val=0; 
#ifdef DEBUG
	      (*debug) << " node: " << node << " val: " << val << " number: " << j << std::endl;
#endif       
	  algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1],
				val, dofspernode_, SYSTEM);
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      nodes=ptBCs_->GetNodesLevel(bcs_id_[i]); 

      val=val_id_[i]*val_tfunc;
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
#ifdef DEBUG
	  (*debug) << " node: " << node << " val: " << val 
		   << " number: " << j << "PDEnode: " <<  Mesh2PDENode_[node-1] <<  std::endl;
#endif      

	  // Mesh node numbers are mapped to PDE node numbers
	  algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1],
				val, dofspernode_, SYSTEM);
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
  conf->ifgetliststr("loads",bcs_loads_,eq);

  Integer i;

  val_id_.resize(bcs_id_.size());

  for(i=0; i<bcs_id_.size(); i++)
    conf->get(bcs_id_[i],val_id_[i],eq,"bc_conditions","inhomogeneous_dirichlet");

  val_loads_.resize(bcs_loads_.size());
  for(i=0; i<bcs_loads_.size(); i++)
    conf->get(bcs_loads_[i],val_loads_[i],eq,"bc_conditions","loads");

}


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

#ifdef DEBUG
//   (*debug) << "--------------------" << std::endl;
//   (*debug) << " Mesh2PDENode()" << std::endl;
//   for (Integer i=0; i<MeshNodes.size(); i++)
//     (*debug) << "in: " << MeshNodes[i] << " out: " << PDENodes[i] << std::endl;
#endif
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
  NumPDENodes_ = PDE2MeshNode_.size();
  
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

  NumPDENodes_ = PDE2MeshNode_.size();
  //  (*cla) << "Mesh2PDENod   " << Mesh2PDENode << std::endl;
#ifdef DEBUG
  (*debug) << "Mesh2PDENodes:" << std::endl << Mesh2PDENode << std::endl;
#endif
}

  
void BasePDE::AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			  std::vector<Integer> & PDE2MeshNode,
			  const std::vector<Elem*> & Elements)
{
#ifdef TRACE
  (*trace) << "entering BasePDE:AssignPDENodeNumbers:" << std::endl;
#endif

  // Initialize Mesh2PDENode_ and PDE2MeshNode
  Mesh2PDENode.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  Integer NodeCounter = 1;
  
  // Iterate over all elements 
  for (Integer j=0; j<Elements.size(); j++)
    {
      // Iterate over all element nodes
      for (Integer NumNodes=0; NumNodes<Elements[j]->connect.size(); NumNodes++)
	{
	  // Check if node was already assigned
	  if (Mesh2PDENode[Elements[j]->connect[NumNodes] - 1] == -1)
	    {
	      Mesh2PDENode[Elements[j]->connect[NumNodes] - 1] = NodeCounter;
	      PDE2MeshNode.push_back(Elements[j]->connect[NumNodes]);
	      NodeCounter++;
	    }
	}
    }
}


void BasePDE::GetElemCoords(const Vector<Integer> connect, Matrix<Double> &coordMat, const Integer level)
{
#ifdef TRACE
  (*trace) << "entering BasePDE:GetElemCoords:" << std::endl;
#endif

  ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
  if (deltCoords_.size() != 0)
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

  std::vector<Elem*> Elems;

  MeshSol.reshape(PDESol.dim(), ptgrid_->GetMaxnumElem(actlevel_));


  // loop over all SubDomains
  for (Integer isd=0; isd<SD.size(); isd++)
  {
    ptgrid_->GetElemSD(Elems, SD[isd], actlevel_);
    
    // loop over all dims
    for (Integer dim=0; dim<PDESol.dim(); dim++)

	// loop over all elements
      for (Integer i=0; i<Elems.size(); i++) {

   	  MeshSol[dim][i] = PDESol[dim][i];
      }
	
  }

}


void BasePDE::NodeSolutionToCoupling(Array<Double>& CouplingSol,
				     const std::vector<Integer>& NodeNumbers)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformElemSolution" << std::endl;
#endif
  
  CouplingSol.reshape(dofspernode_, NodeNumbers.size());
  
  //std::cerr << "In " << pdename_ << "-NodeSolutiontoCoupling" << std::endl;
  //std::cerr << "CouplingSol size:" << CouplingSol.size() << ", dim: " << CouplingSol.dim() << std::endl;
  //std::cerr << "sol size:" << sol_.size() << ", dim: " << sol_.dim() << std::endl;
  //std::cerr << "Mesh2PDENode_.size = " << Mesh2PDENode_.size() << std::endl;
  //std::cerr << "NumPDENodes = " << NumPDENodes_ << std::endl;
  
  for (Integer i=0; i<CouplingSol.dim(); i++)
    for (Integer j=0; j<CouplingSol.size(); j++)
      {
	//std::cerr << "processing dim: " << i <<", j:" << j << std::endl; 
	CouplingSol[i][j] = sol_[i][Mesh2PDENode_[NodeNumbers[j]-1 ] - 1];
      }
  
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


void BasePDE::SetAlgSys_id(const Integer as_sysid)
{
  as_sysid_ = as_sysid;
}


void BasePDE::MassMultiDof(Matrix<Double>& massMultDof, const Matrix<Double>& massMatSingleDof,  const Integer nrDofs)
{
    
#ifdef TRACE
    (*trace) << "entering BasePDE::MassMatMultiDof" << std::endl;
#endif
    
    const Integer singleDofSize = massMatSingleDof.getSize();
    const Integer multDofSize = singleDofSize * nrDofs;
    
    Integer i, j, actDof;
    
    massMultDof.Resize(multDofSize);
    massMultDof.Init();
    
    for (i=0; i < singleDofSize; i++)
      for (j=0; j < singleDofSize; j++)
	for (actDof=0; actDof < nrDofs; actDof++)
	  massMultDof[i*nrDofs + actDof][j*nrDofs + actDof] = massMatSingleDof[i][j]; 
}

BasePDE::~BasePDE()
{
#ifdef TRACE
  (*trace) << " entering BasePDE::~BasePDE() " << std::endl;
#endif

  if (algsys_) delete algsys_;
}


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
  
  (*infofile) << " tolerance space error: " << tolSpaceErr_ << std::endl;
 
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
  
  (*infofile) << " space error: " << totalErr <<
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

} // end of namespace

#endif //#ifndef NEWBASEPDE
