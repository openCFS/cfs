//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>

#include "domain.hh"
#include "elec2dPDE.hh"

namespace CoupledField
{

template <class Dim>
Domain<Dim> :: Domain(FileType * const aptFileType, OutResultUnverg * ptUnverg,  
                      Material * materialdata,  Grid<Dim> * ptgrid)
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif
 InFile     = aptFileType; 
 OutFile    = ptUnverg;
 ptmaterial = materialdata;

 grid = ptgrid;

 //read in the mesh information
 grid->Read();

}

template <class Dim>
Domain<Dim> :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << std::endl;
#endif

  ;
}

template<class Dim>
void Domain<Dim> :: InitPDE()
{
#ifdef TRACE
  (*trace) << "entering Domain::InitPDE" << std::endl;
#endif

  // get numbers of PDEs in domain
  numpde = 1;

  //allcocate all specific PDEs
  for (int i=0;i<numpde;i++)
    {
      ptpde[i] = new Elec2dPDE(grid,ptmaterial,InFile);
    }

}

template<class Dim>
void Domain<Dim> :: InitAlgSys(AbstractAlgebraicSys * AS)
{
#ifdef TRACE
  (*trace) << "entering Domain::InitAlgSys" << std::endl;
#endif

  Integer level=0;

  //set pointer to Algebraic system
  ptalgsys = AS;


  //check, how much systems are needed and how much matrix graphs
  numsys   = 1;
  numgraph = 1;
  AS->InitAlgSys(numsys, numgraph);


  //set solver parameters
  Double eps, dampiter;
  Integer maxnumit;
  Integer solvertype;
  Integer precondtype;
  Integer nsys;

  for (nsys=1;nsys<=numsys;nsys++)
    {
      ptpde[nsys-1]->SpecifySolver(solvertype,precondtype,eps,dampiter,maxnumit);
      AS->SetSolverParameter(nsys,eps,dampiter,maxnumit,solvertype,precondtype);
    }

  //init the algsys-graph
  Integer numnode = grid->GetMaxnumnodes(level);
  std::cout << "Number of Nodes: " << numnode << endl;

  AS->InitAlgSysGraph(numnode);

 // get the graph - connectivity matrix
  Integer nel, numelem, elemsize, pos[30];

  numelem = grid->GetMaxnumElem(level);

  for (nsys=1;nsys<=numsys;nsys++)
    {
      for (nel=1; nel<=numelem; nel++)
	{
	  elemsize = grid->GetNumNodesPerElem(nel-1,level);
          grid->GetConnection(pos, level, nel-1, elemsize); // global node numbers! pos[0],...,pos[elemsize-1]
	  AS->SetAlgSysGraph(pos,elemsize,nsys);
	}
    }

  //now we can create all the necessary matrices
  Integer matrixtype;
  Integer matrixsystype[6];    
  Integer graphtype; 
  Integer numdofpernode;
  Integer numdirichlets;
  Integer numconstraints;

  for (nsys=1;nsys<=numsys;nsys++)
    {
      ptpde[nsys-1]->SpecifyMatrices(matrixtype, matrixsystype, graphtype, numdofpernode, 
                                     numdirichlets, numconstraints);
      AS->CreateAlgSysMatrices(matrixtype, matrixsystype, graphtype, numdofpernode, 
                               numdirichlets, numconstraints);
    }
}


template<class Dim>
void Domain<Dim> :: PrintGrid(Integer level)
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintGrid" << std::endl;
#endif

 OutFile->Create(grid,level);

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
