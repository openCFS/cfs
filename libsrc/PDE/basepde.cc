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

  //standard paramewter for solver
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

  //set the graph type used for the system matrices
  Integer numnode = ptgrid_->GetMaxnumnodes(actlevel_);
  SetupMatrixGraph(numnode,NODEGRAPH);

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
      // load material data into array "matData"
      conf->getstr(subdoms_[i], matName, "list_subdomains");
      std::cout << "matname:" << matName << "  sub:" << subdoms_[i] << std::endl;	
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
	  connecth=elemssd[iel]->connect;
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
  Integer numdirichlets;
  Integer numconstraints;

  //ask the PDE discrete form
  DiscreteParamsPDE();

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
              algsys_->SetDirichlet(j+1, node, val, dofspernode_, SYSTEM);
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
              algsys_->SetDirichlet(j+1, node, val, dofspernode_, SYSTEM);
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

  conf->getliststr("homogeneous_dirichlet",bcs_hd_,eq); 
  conf->getliststr("inhomogeneous_dirichlet",bcs_id_,eq);

  Integer i;

  val_id_.resize(bcs_id_.size());

  for(i=0; i<bcs_id_.size(); i++)
    {
    conf->get(bcs_id_[i],val_id_[i],eq,"bc_conditions","inhomogeneous_dirichlet");
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


BasePDE::~BasePDE()
{
#ifdef TRACE
  (*trace) << " entering BasePDE::~BasePDE() " << std::endl;
#endif

}


} // end of namespace
