#include <iostream>
#include <fstream>
#include <string>

#include "domain.hh"
#include "interface_piles.hh"
#include "elec2dPDE.hh"
#include "acoustic2dPDE.hh"

namespace CoupledField
{

template <class Dim>
Domain<Dim> :: Domain(FileType * const aptFileType, WriteResults<Dim> * ptOut,  Material * materialdata, TimeFunc * aptTimeFunc) 
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif
 InFile     = aptFileType; 
 OutFile    = ptOut;
 ptmaterial = materialdata;
 ptTimeFunc = aptTimeFunc;

  // read type of output results from conf-file
  std::string libmesh;
  conf->get("mesh_library",libmesh);

 // initialize pointer to grid 
 //  if (libmesh =="gridlib") ptgrid=new InterfaceGridlib<Dim>(InFile);
//  else
  if (libmesh =="cfsgrid") ptgrid=new GridInterfaceCFS<Dim>(InFile);
   else
     Error("Unknown type of mesh_library in conf-file",__FILE__,__LINE__);

  // allocate an object with an information about boundary condition
  ptBCs=new BCs(InFile);

 //read in the mesh information
 ptgrid->Read();

 //read restraints information
 ptBCs->ReadBCs();

 InitPDE();
 InitAlgSys();
}

template <class Dim>
Domain<Dim> :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << std::endl;
#endif

  if (!ptgrid) delete ptgrid;
  if (!ptBCs) delete ptBCs;
  if (!ptalgsys) delete ptalgsys;
}

template<class Dim>
void Domain<Dim> :: InitPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDE" << std::endl;
#endif

  // get numbers of PDEs in domain
  numpde = 1;

    //get analysis type
  Integer soltype;
  Integer statickey;
  InFile->preReadTransAnal(soltype, statickey);

  Integer level = 0;
  Integer numnode = ptgrid->GetMaxnumnodes(level);

  // read time step
  Integer numsteps;
  Double dt0;
  InFile->ReadNumStepsAndTimeSteps(numsteps, dt0);

  //allocate all specific PDEs
  ptpde[0]=new Acoustic2dPDE(ptgrid,ptmaterial,ptTimeFunc,InFile,OutFile);
  for (int i=0;i<numpde;i++)
    {
      //      ptpde[i] = new Elec2dPDE(grid,ptmaterial,InFile,OutFile,statickey);      
 //     ptpde[i] = new Therm2dPDE(ptgrid,ptmaterial,ptTimeFunc,InFile,OutFile,statickey,numnode);

    }

   
   ptpde[0]->SetStepData();

/*
  Integer * stepdata = new Integer[4];
  Double  * integrationdata = new Double[5];

  if (statickey==0)
    {
      InFile->ReadTransAnalTran(stepdata,integrationdata);
      for (int i=0;i<numpde;i++)
	{
	  ptpde[i]->SetStepData(stepdata,integrationdata);
	}
    }

  delete [] stepdata;
  delete [] integrationdata;
*/

}

template<class Dim>
void Domain<Dim> :: InitAlgSys()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitAlgSys" << std::endl;
#endif

  //set pointer to Algebraic system
  ptalgsys = new AlgSysPILES();


  //check, how much systems are needed and how much matrix graphs
  numsys   = 1;
  numgraph = 1;
  ptalgsys->InitAlgSys(numsys, numgraph);

    //set solver parameters
  Double eps, dampiter;
  Integer maxnumit;
  Integer solvertype;
  Integer precondtype;
  Integer insys;

  for (insys=0;insys<numsys;insys++)
    {
      ptpde[insys]->SetAlgSys_id(insys);
      ptpde[insys]->SpecifySolver(solvertype,precondtype,eps,dampiter,maxnumit);
      ptalgsys->SetSolverParameter(insys,eps,dampiter,maxnumit,solvertype,precondtype);
    }

  //init the algsys-graph
  Integer level=0;
  Integer numnode = ptgrid->GetMaxnumnodes(level);
  std::cout << "numnode" << numnode << std::endl;

  //for each system: first diagonal blocks and then off-diagonalblocks
  for (insys=0;insys<numsys;insys++)
   {
     ptalgsys->InitAlgSysGraph(numnode,insys,insys);
   }

 // get the graph - connectivity matrix
  Integer nel, numelem, elemsize, pos[30];

  numelem = ptgrid->GetMaxnumElem(level);

  for (insys=0;insys<numsys;insys++)
    {
      for (nel=0; nel<numelem; nel++)
	{
	  elemsize = ptgrid->GetNumNodesPerElem(nel,level);
          ptgrid->GetConnection(pos, level, nel, elemsize); // global node numbers! pos[0],...,pos[elemsize-1]
	  ptalgsys->SetAlgSysGraph(pos,elemsize,insys,insys);
	}
    }

  //now we can create all the necessary matrices
  Integer matrixtype;
  Integer matrixsystype[5];    
  Integer graphtype; 
  Integer numdofpernode;
  Integer numdirichlets;
  Integer numconstraints;

  for (insys=0;insys<numsys;insys++)
    {
      ptpde[insys]->SpecifyMatrices(matrixtype, matrixsystype, graphtype, numdofpernode, 
                                     numdirichlets, numconstraints);
      numdirichlets = ptBCs->GetNumDirichlet(level);

      ptalgsys->CreateAlgSysMatrices(insys,insys,matrixsystype,matrixtype,graphtype, numdofpernode, 
                               numdirichlets, numconstraints);
    }

  //now reset AlgebraicSystem 
  //matrix_id = 1: system matrix
  Integer matrix_id = 1;
  for (insys=0;insys<numsys;insys++)
    {
      ptalgsys->ResetAlgSys(insys,insys,matrix_id);
    }
}

template<class Dim>
void Domain<Dim> :: PrintGrid(Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintGrid" << std::endl;
#endif

 OutFile->Init(ptgrid);
 OutFile->WriteGrid(level);
}

template<class Dim>
void Domain<Dim> :: SetDBC(Integer apde, Integer level, Integer update)
{
#ifdef TRACE
  (*trace) << "entering Domain::SetDBC" << std::endl;
#endif

  ptpde[apde]->SetBCs(ptalgsys, ptBCs, level, update,0);
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

template<class Dim>
void Domain<Dim> :: SetSubdomains()
{
#ifdef TRACE
  (*trace) << "entering Domain::SetSubdomains" << std::endl;
#endif

}

template<class Dim>
void Domain<Dim> :: PrintDomain()
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintDomain" << std::endl;
#endif
}

}
