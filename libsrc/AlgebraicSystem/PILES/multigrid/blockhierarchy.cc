#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RBlockHierarchy :: RBlockHierarchy(BaseParameter & param, Integer asize, Integer aoffset, 
				   Integer anumrhs, Integer adof)
  : BaseHierarchy(param, asize, aoffset, anumrhs, adof)
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::RBlockHierarchy" << endl;
#endif

}
  
RBlockHierarchy :: ~RBlockHierarchy()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::~RBlockHierarchy" << endl;
#endif
  
}

void RBlockHierarchy :: Init()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::Init" << endl;
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
  delete node[level].auxpro;
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
      delete node[level].auxpro;
      node[level].auxpro = NULL;
      node[level].auxres = NULL;

      delete node[level].smooth;
      node[level].smooth = NULL;
    }  

  amglevel = 0;
}

void RBlockHierarchy :: DeleteAuxMemory()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::DeleteAuxMemory" << endl;
#endif

  //delete node[level].coarse;
  node[level].coarse = NULL;
}

void RBlockHierarchy :: MakeSmoother()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeSmoother" << endl;
#endif

  switch (param->GetSmoothingType())
    {
    case 1:
      if (param->GetSmoothingStepFor() == 1)
	{
	  node[level].smooth = new RBlockFastGSSmoother(*param,fsize_sys,numrhs,dof);
	}
      else
	{
	  node[level].smooth = new RBlockGSSmoother(*param,fsize_sys,numrhs,dof);
	}
      break;

    case 2:
      //node[level].smooth = new RBlockFastGSSmoother(*param,fsize_sys,numrhs,dof);
      ;
      break;

    default:
      node[level].smooth = new RBlockGSSmoother(*param,fsize_sys,numrhs,dof);
    }
}

void RBlockHierarchy :: MakeCoarsening()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeCoarsening" << endl;
#endif

  switch (param->GetCoarsening())
    {
    case 1:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
      break;

    case 2:
      //node[level].coarse = new RBlockCoarse();
      ;
      break;

    default:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
    }
}

void RBlockHierarchy :: MakeSysTransfer()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeSysTransfer" << endl;
#endif

  switch (param->GetTransfer())
    {
      case 1:
	node[level].syspro = new RBlockTransfer(fsize_sys, dof, rsw);
	node[level].sysres = node[level].syspro;
      break;

    case 2:
      //node[level].syspro = new RBlockProlongation();
      //node[level].sysres = node[level].syspro;
      ;
      break;

    default:
      node[level].syspro = new RBlockTransfer(fsize_sys, dof, rsw);
      node[level].sysres = node[level].syspro;
    }
}

void RBlockHierarchy :: MakeAuxTransfer()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeAuxTransfer" << endl;
#endif

  Boolean mem;

  if (param->GetAuxiliaryMatrix() == 0)
    {
      mem = FALSE;
    }
  else
    {
      mem = TRUE;
    }

  switch (param->GetTransfer())
    {
      case 1:
	node[level].auxpro = new RScalarTransfer(fsize_sys, rsw, mem);
	node[level].auxres = node[level].auxpro;
      break;

    case 2:
      //node[level].auxpro = new RBlockProlongation();
      //node[level].auxres = node[level].auxpro;
      ;
      break;

    default:
      node[level].auxpro = new RScalarTransfer(fsize_sys, rsw, mem);
      node[level].auxres = node[level].auxpro;
    }
}

void RBlockHierarchy :: MakeSysMatrix()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeSysMatrix" << endl;
#endif

  dir_sys = 0;
  
  node[level-1].sysmat = new RBlockMatrix(csize_sys, cnne_sys, dof, dir_sys); 
}

void RBlockHierarchy :: MakeAuxMatrix()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeAuxMatrix" << endl;
#endif

  dir_aux = 0;

  node[level-1].auxmat = new RScalarMatrix(csize_aux, cnne_aux, dir_aux); 

  if (param->GetAuxiliaryMatrix() == 0)
    {
      node[level-1].auxmat->Calculated();
      node[level-1].auxmat->SetGraph();
    }
}

void RBlockHierarchy :: MakeVector()
{
#ifdef TRACE
  (*trace) << "entering RBlockHierarchy::MakeVector" << endl;
#endif

  node[level].vec1 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec2 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec3 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec4 = new RealVector(fsize_sys, numrhs, dof); 
}

}
