#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
#include <DataInOut/conffile.hh>
 
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

  //standard parameter for solver
  eps_         = 1.0e-8;
  dampiter_    = 0.7;
  maxnumiter_  = 100;
  solvertype_  = RealDirect;
  precondtype_ = ID;
  numeqcoarse_ = 200;
  coarsealpha_ = 0.1;

  //get analysis type
  std::string analysis;
  conf->get("analysis", analysis);

  if (analysis=="static") 
    analysistype_ = STATIC;
  else if (analysis=="transient")
    analysistype_ = TRANSIENT;
  else
    Error("Analysis Type not supported",__FILE__,__LINE__);
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
	  Mesh2PDENode(connecth,elemssd[iel]->connect);
	  fe_type=elemssd[iel]->ptElem->feType();
	  algsys_->SetElementPos(connecth.get(),connecth.size(),fe_type);
#ifdef DEBUG
	  (*debug) << "Nodes to AlgSys, Element: " << iel+1 << std::endl;
	  (*debug) << connecth << std::endl;
#endif
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
  Integer numdirichlets;
  Integer numconstraints;

  //ask the PDE discrete form
  // DiscreteParamsPDE();

  numdofpernode  = dofspernode_; 
  numdirichlets  = GetNumRestraints(actlevel_);
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

   //  SpecifyMatrices(matrixclass, matrixsystype, graphtype, numdofpernode,numdirichlets, numconstraints);

  //put to algebraic system
  algsys_->CreateMatrix(matrixsystype, matrixclass, graphtype, numdofpernode,  numdirichlets,
			numconstraints);

  //create solver and preconditioner
  algsys_->CreateSolver();
  algsys_->CreatePrecond(matrixclass);

  //now reset AlgebraicSystem 
  algsys_->InitRHS();
  algsys_->InitSol();

  for (Integer i=0;i<5;i++)
    {
      if (matrixsystype[i] !=0)
 	algsys_->InitMatrix(i+1);
    }

}


void  BasePDE::SetBCs(const Integer level, const Integer update, const Double time)
{
#ifdef TRACE
  (*trace) << "entering  BasePDE::SetBCs" << std::endl;
#endif

  Integer node;
  Double val, val_tfunc;

  val_tfunc = 1.0;
  if (ptTimeFunc_->GetmaxTimeFnc()!=0)
      val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,level);

  Integer j=0;
  std::list<Integer> nodes;

  if (InfoPrint)
    (*infofile) << " ---------------- Dirichle boundary condition -------------" << std::endl;

  Integer i;
  for (i=0; i<bcs_hd_.size(); i++)
    {  
      nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
	  val=0; 
          if (update==1)
            {
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
              algsys_->UpdateDirichlet(j+1, val, SYSTEM);
            }
          else
            {
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
	      // Mesh node numbers are mapped to PDE node numbers
              algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, dofspernode_, SYSTEM);
            }
	}
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      nodes=ptBCs_->GetNodesLevel(bcs_id_[i]); 

      val=val_id_[i]*val_tfunc;
      for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++, j++)
	{
	  node=*p;
          if (update==1)
            {	     
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;

             algsys_->UpdateDirichlet(j+1, val, SYSTEM);
            }
          else
            {
	      if (InfoPrint)
		(*infofile) << " node: " << node << " val: " << val << std::endl;
	      // Mesh node numbers are mapped to PDE node numbers
              algsys_->SetDirichlet(j+1, Mesh2PDENode_[node-1], val, dofspernode_, SYSTEM);
            }
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



void BasePDE::Mesh2PDENode(Vector<Integer> & PDENodes, Vector<Integer> & MeshNodes)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::Mesh2PDENode " << std::endl;
#endif

  PDENodes.Resize(MeshNodes.size());
  
  for (Integer i=0; i<MeshNodes.size(); i++)
    PDENodes[i] = Mesh2PDENode_[MeshNodes[i]-1];

#ifdef DEBUG
  (*debug) << "--------------------" << std::endl;
  (*debug) << " Mesh2PDENode()" << std::endl;
  for (Integer i=0; i<MeshNodes.size(); i++)
    (*debug) << "in: " << MeshNodes[i] << " out: " << PDENodes[i] << std::endl;
#endif
}

void BasePDE::PDE2MeshNode(Vector<Integer> & MeshNodes, Vector<Integer> & PDENodes)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::PDE2MeshNode " << std::endl;
#endif

  MeshNodes.Resize(PDENodes.size());

  for (Integer i=0; i<PDENodes.size(); i++)
    MeshNodes[i] = PDE2MeshNode_[PDENodes[i]-1];

#ifdef DEBUG
  (*debug) << "--------------------" << std::endl;
  (*debug) << " PDE2MeshNode()" << std::endl;
  for (Integer i=0; i<PDENodes.size(); i++)
    (*debug) << "in: " << PDENodes[i] << " out: " << MeshNodes[i] << std::endl;
#endif
}


void BasePDE::AssignPDENodeNumbers()
{
#ifdef TRACE
  (*trace) << "entering BasePDE:AssignPDENodeNumbers:" << std::endl;
#endif

  // Initialize Mesh2PDENode_ and PDE2MeshNode
  Mesh2PDENode_.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
  std::vector<Elem*> SD;
  Integer NodeCounter = 1;

//  PDE2MeshNode_.resize(ptgrid_->GetMaxnumnodes(actlevel_),-1);
//   for (Integer i=0;i<ptgrid_->GetMaxnumnodes(actlevel_);i++)
//     {
//       Mesh2PDENode_[i] = i+1;
//       PDE2MeshNode_[i] = i+1;
//     }
  
  // Iterate over Subdomains
  for (Integer i=0; i<subdoms_.size(); i++)
    {
      ptgrid_->GetElemSD(SD,subdoms_[i],actlevel_);
      // Iterate over all elements in subdomain
      for (Integer j=0; j<SD.size(); j++)
	{
	  // Iterate over all element nodes
	  for (Integer NumNodes=0; NumNodes<SD[j]->connect.size(); NumNodes++)
	    {
	      // Check if node was already assigned
	      if (Mesh2PDENode_[SD[j]->connect[NumNodes] - 1] == -1)
		{
		  Mesh2PDENode_[SD[j]->connect[NumNodes] - 1] = NodeCounter;
		  PDE2MeshNode_.push_back(SD[j]->connect[NumNodes]);
		  NodeCounter++;
		}
	    }
	}
    }

  NumPDENodes_ = PDE2MeshNode_.size();

  //std::cout << "NumPDENodes = " << NumPDENodes_ << std::endl;
  //std::cout << "Mesh2PDENode_" << std::endl;
  //for (Integer i=0; i<Mesh2PDENode_.size(); i++)
  //  std::cout << i << ":" <<Mesh2PDENode_[i] << std::endl;

  //std::cout << "PDE2MeshNode" << std::endl;
  //for (Integer i=0; i<PDE2MeshNode_.size(); i++)
  //  std::cout << i << ": " << PDE2MeshNode_[i] << std::endl;
  
  

#ifdef TRACE
  (*trace) << " leaving BasePDE::AssignPDENodeNumbers" << std::endl;
#endif 
}

void BasePDE::TransformNodeSolution(Vector<Double> & MeshSol,  Vector<Double> & PDESol)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformNodeSolution" << std::endl;
#endif

  Integer k, node, idx;
  MeshSol.Resize(ptgrid_->GetMaxnumnodes(actlevel_)*dofspernode_);

  k=0;
  for (Integer i=0; i<PDE2MeshNode_.size(); i++)
    {
      node = PDE2MeshNode_[i];
      idx  = dofspernode_*(node-1);
      for (Integer j=0;j<dofspernode_; j++)
	{
	  MeshSol[idx+j] = PDESol[k];
	  k++;
	}
    }
}

void BasePDE::TransformElemSolution(Vector<Double> & MeshSol,  Vector<Double> & PDESol)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::TransformElemSolution" << std::endl;
#endif
  MeshSol.Resize(ptgrid_->GetMaxnumElem(actlevel_));
  //std::cerr << "MeshSol.size() = " << MeshSol.size() << std::endl;
  //std::cerr << "PDESol.size() = " << PDESol.size() << std::endl;
  for (Integer i=0; i<PDESol.size(); i++)
    MeshSol[i] = PDESol[i];
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
      res+=ptBCs_->GetNumNodesLevel(bcs_hd_[i],level);
    }

  for (i=0; i<bcs_id_.size(); i++)
    {
      res+=ptBCs_->GetNumNodesLevel(bcs_id_[i],level);
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

}


} // end of namespace
