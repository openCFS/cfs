#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

OutbackHierarchy :: OutbackHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs)
  : BaseHierarchy(param, asize, aoffset, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::OutbackHierarchy" << endl;
#endif

}
  
OutbackHierarchy :: ~OutbackHierarchy()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::~OutbackHierarchy" << endl;
#endif
  
}

void OutbackHierarchy :: Init()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::Init" << endl;
#endif
  
  Integer i;

  level = offset;

  delete node[level].vec1;
  delete node[level].vec2;
  delete node[level].vec3;
  delete node[level].vec4;
  node[level].vec1 = NULL;
  node[level].vec2 = NULL;
  node[level].vec3 = NULL;
  node[level].vec4 = NULL;
  
  delete node[level].syspro;
  node[level].syspro = NULL;
  node[level].sysres = NULL;
  node[level].auxpro = NULL;
  node[level].auxres = NULL;
  
  delete node[level].smooth;
  node[level].smooth = NULL;

  for (i=1; i<amglevel; i++)
    {
      level = offset-i;

      delete node[level].sysmat;
      node[level].sysmat = NULL;
      delete node[level].auxmat;
      node[level].auxmat = NULL;

      delete node[level].vec1;
      delete node[level].vec2;
      delete node[level].vec3;
      delete node[level].vec4;
      node[level].vec1 = NULL;
      node[level].vec2 = NULL;
      node[level].vec3 = NULL;
      node[level].vec4 = NULL;

      delete node[level].syspro;
      node[level].syspro = NULL;
      node[level].sysres = NULL;
      node[level].auxpro = NULL;
      node[level].auxres = NULL;

      delete node[level].smooth;
      node[level].smooth = NULL;
    }  

  amglevel = 0;
}

void OutbackHierarchy :: DeleteAuxMemory()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::DeleteAuxMemory" << endl;
#endif
  
  delete node[level].coarse;
  node[level].coarse = NULL;
}

void OutbackHierarchy :: MakeSmoother()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeSmoother" << endl;
#endif

  switch (param->GetSmoothingType())
    {
    case 1:
      node[level].smooth = new OutbackSmoother(*param,fsize_sys,numrhs);
      break;

    case 2:
      ;
      break;

    default:
      node[level].smooth = new OutbackSmoother(*param,fsize_sys,numrhs);
    }
}

void OutbackHierarchy :: MakeCoarsening()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeCoarsening" << endl;
#endif

  switch (param->GetCoarsening())
    {
    case 1:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
      break;

    case 2:
      //node[level].coarse = new RScalarCoarse();
      ;
      break;

    default:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
    }
}

void OutbackHierarchy :: MakeSysTransfer()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeSysTransfer" << endl;
#endif

  switch (param->GetTransfer())
    {
      case 1:
	node[level].syspro = node[level].auxpro;//new RScalarTransfer(fsize_sys, rsw);
	node[level].sysres = node[level].syspro;
      break;

    case 2:
      //node[level].syspro = new RScalarProlongation();
      //node[level].sysres = node[level].syspro;
      ;
      break;

    default:
      node[level].syspro = node[level].auxpro;//new RScalarTransfer(fsize_sys, rsw);
      node[level].sysres = node[level].syspro;
    }
}

void OutbackHierarchy :: MakeAuxTransfer()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeAuxTransfer" << endl;
#endif

  switch (param->GetTransfer())
    {
      case 1:
	node[level].auxpro = new RScalarTransfer(fsize_sys, rsw);
	node[level].auxres = node[level].auxpro;
      break;

    case 2:
      //node[level].auxpro = new RScalarProlongation();
      //node[level].auxres = node[level].auxpro;
      ;
      break;

    default:
      node[level].auxpro = new RScalarTransfer(fsize_sys, rsw);
      node[level].auxres = node[level].auxpro;
    }
}

void OutbackHierarchy :: MakeSysMatrix()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeSysMatrix" << endl;
#endif

  dir_sys = 0;

  node[level-1].sysmat = new OutbackMatrix(csize_sys, size, dir_sys, offset-level+1);
}

void OutbackHierarchy :: MakeAuxMatrix()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeAuxMatrix" << endl;
#endif

  dir_aux = 0;

  node[level-1].auxmat = new RScalarMatrix(csize_aux, cnne_aux, dir_aux); 
}

void OutbackHierarchy :: MakeVector()
{
#ifdef TRACE
  (*trace) << "entering OutbackHierarchy::MakeVector" << endl;
#endif

  node[level].vec1 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec2 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec3 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec4 = new RealVector(fsize_sys, numrhs, dof); 

}

}
