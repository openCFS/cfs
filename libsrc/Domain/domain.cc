#include <iostream>
#include <fstream>
#include <string>

#include "domain.hh"
#include "interface_piles.hh"
#include "elec2dPDE.hh"
#include "therm2dPDE.hh"
#include "acoustic2dPDE.hh"

namespace CoupledField
{

template <class Dim>
Domain<Dim> :: Domain(FileType * const aptFileType, WriteResults<Dim> * ptOut,  Material * materialdata, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif

 InFile_     = aptFileType; 
 OutFile_    = ptOut;
 ptmaterial_ = materialdata;
 ptTimeFunc_ = aptTimeFunc;

  // read type of output results from conf-file
  std::string libmesh;
  conf->get("mesh_library",libmesh);

 // initialize pointer to grid 
 //  if (libmesh =="gridlib") ptgrid=new InterfaceGridlib<Dim>(InFile);
//  else
  if (libmesh =="cfsgrid") ptgrid_=new GridInterfaceCFS<Dim>(InFile_);
   else
     Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);

  // allocate an object with an information about boundary condition
  ptBCs_=new BCs(InFile_);

 //read in the mesh information
 ptgrid_->Read();

 //read restraints information
 ptBCs_->ReadBCs();

 //set pointer to Algebraic system
 ptalgsys_ = new AlgSysPILES();
 if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");

 // it is important this order of these functions
 InitPDE();
 InitAlgSys();
}

template <class Dim>
Domain<Dim> :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << std::endl;
#endif

  if (!ptgrid_) delete ptgrid_;
  if (!ptBCs_) delete ptBCs_;
  if (!ptalgsys_) delete ptalgsys_;

  Integer i;
  if (!ptpde_) 
     {
      for (i=0; i<numpde_; i++)
        if (!ptpde_[i]) delete ptpde_[i];
      delete [] ptpde_;
     } 
}

template<class Dim>
void Domain<Dim> :: InitPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDE" << std::endl;
#endif

  // get numbers of PDEs in domain
  numpde_ = 1;

    //get analysis type
  Integer soltype;
  Integer statickey;
  InFile_->preReadTransAnal(soltype, statickey);

  Integer level = 0;
  Integer numnode = ptgrid_->GetMaxnumnodes(level);

  // read time step
  Integer numsteps;
  Double dt0;
  InFile_->ReadNumStepsAndTimeSteps(numsteps, dt0);

  //allocate all specific PDEs
  if (!ptalgsys_) Error("You try to allocate object BasePDE with null pointer to AlgSys");

//  ptpde_[0]=new Therm2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);

  ptpde_[0]=new Acoustic2dPDE(ptalgsys_,ptgrid_,ptmaterial_,ptTimeFunc_,InFile_,OutFile_);
  for (int i=0;i<numpde_;i++)
    {
      //      ptpde[i] = new Elec2dPDE(grid,ptmaterial,InFile,OutFile,statickey);      
 //     ptpde[i] = new Therm2dPDE(ptgrid,ptmaterial,ptTimeFunc,InFile,OutFile,statickey,numnode);
    }

//   ptpde_[0]->SetStepData();

}

template<class Dim>
void Domain<Dim> :: InitAlgSys()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitAlgSys" << std::endl;
#endif

/*
  //set pointer to Algebraic system
  ptalgsys_ = new AlgSysPILES();
  if (!ptalgsys_) Error("Can't allocate memory for algebraic system Piles");
*/

  //check, how much systems are needed and how much matrix graphs
  numsys_   = 1;
  numgraph_ = 1;
  ptalgsys_->InitAlgSys(numsys_, numgraph_);

    //set solver parameters
  Double eps, dampiter;
  Integer maxnumit;
  Integer solvertype;
  Integer precondtype;
  Integer insys;
  Integer numeqcoarse;

  for (insys=0;insys<numsys_;insys++)
    {
      ptpde_[insys]->SetAlgSys_id(insys);
      ptpde_[insys]->SpecifySolver(solvertype,precondtype,eps,dampiter,maxnumit,numeqcoarse);
      ptalgsys_->SetSolverParameter(insys,eps,dampiter,maxnumit,solvertype,precondtype,  numeqcoarse);
    }

  //init the algsys-graph
  Integer level=0;
  Integer numnode = ptgrid_->GetMaxnumnodes(level);

  //for each system: first diagonal blocks and then off-diagonalblocks
  for (insys=0;insys<numsys_;insys++)
   {
     ptalgsys_->InitAlgSysGraph(numnode,insys,insys);
   }

 // get the graph - connectivity matrix
  Integer nel, numelem, elemsize, pos[30];

  numelem = ptgrid_->GetMaxnumElem(level);

  for (insys=0;insys<numsys_;insys++)
    {
      for (nel=0; nel<numelem; nel++)
	{
	  elemsize = ptgrid_->GetNumNodesPerElem(nel,level);
          ptgrid_->GetConnection(pos, level, nel, elemsize); // global node numbers! pos[0],...,pos[elemsize-1]
	  ptalgsys_->SetAlgSysGraph(pos,elemsize,insys,insys);
	}
    }

  //now we can create all the necessary matrices
  Integer matrixtype;
  Integer matrixsystype[5];    
  Integer graphtype; 
  Integer numdofpernode;
  Integer numdirichlets;
  Integer numconstraints;

  for (insys=0;insys<numsys_;insys++)
    {
      ptpde_[insys]->SpecifyMatrices(matrixtype, matrixsystype, graphtype, numdofpernode,  numdirichlets, numconstraints);
      numdirichlets = ptBCs_->GetNumRestraints(level);

      ptalgsys_->CreateAlgSysMatrices(insys,insys,matrixsystype,matrixtype,graphtype, numdofpernode,  numdirichlets, numconstraints);
    }

  //now reset AlgebraicSystem 
  //matrix_id = 1: system matrix
  Integer matrix_id = 1;
  for (insys=0;insys<numsys_;insys++)
    {
      ptalgsys_->ResetAlgSys(insys,insys,matrix_id);
    }
}

template<class Dim>
void Domain<Dim> :: PrintGrid(Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintGrid" << std::endl;
#endif

 OutFile_->Init(ptgrid_);
 OutFile_->WriteGrid(level);
}

/*
template<class Dim>
void Domain<Dim> :: SetDBC(Integer apde, Integer level, Integer update)
{
#ifdef TRACE
  (*trace) << "entering Domain::SetDBC" << std::endl;
#endif

  ptpde_[apde]->SetBCs(ptalgsys, ptBCs, level, update,0);
}

template<class Dim>
void Domain<Dim> :: PrintSolution(Double * sol, Integer apde)
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintSolution" << std::endl;
#endif

  Integer level=0;
  Integer numnode = ptgrid->GetMaxnumnodes(level);

  Vector<Double> vecsol(numnode,sol);
  Integer step = 1;
  Double time  = 1.0;  
//  ptpde[apde]->PrintSolution(vecsol,step,time);
}
*/

template<class Dim>
void Domain<Dim> :: SetSubdomains()
{
#ifdef TRACE
  (*trace) << "entering Domain::SetSubdomains" << std::endl;
#endif

}

}
