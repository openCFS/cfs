#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RScalarHierarchy :: RScalarHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs)
  : BaseHierarchy(param, asize, aoffset, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::RScalarHierarchy" << endl;
#endif

}
  
RScalarHierarchy :: ~RScalarHierarchy()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::~RScalarHierarchy" << endl;
#endif
  
}

void RScalarHierarchy :: Init()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::Init" << endl;
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

void RScalarHierarchy :: DeleteAuxMemory()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::DeleteAuxMemory" << endl;
#endif
  
  delete node[level].coarse;
  node[level].coarse = NULL;

  if (node[level].sysmat != node[level].auxmat)
    {
      delete node[level].auxmat;
      node[level].auxmat = NULL;
    }
  else
    {
      node[level].auxmat = NULL;
    }
}

void RScalarHierarchy :: MakeSmoother()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeSmoother" << endl;
#endif

  switch (param->GetSmoothingType())
    {
    case 1:
      if (param->GetSmoothingStepFor() == 1)
	{
	  node[level].smooth = new RScalarFastGSSmoother(*param,fsize_sys,numrhs);
	}
      else
	{
	  node[level].smooth = new RScalarGSSmoother(*param,fsize_sys,numrhs);
	}
      break;

    case 2:
      //node[level].smooth = new RScalarFastGSSmoother(*param,fsize_sys,numrhs);
      ;
      break;

    default:
      node[level].smooth = new RScalarGSSmoother(*param,fsize_sys,numrhs);
    }
}

void RScalarHierarchy :: MakeCoarsening()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeCoarsening" << endl;
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

void RScalarHierarchy :: MakeSysTransfer()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeSysTransfer" << endl;
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

void RScalarHierarchy :: MakeAuxTransfer()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeAuxTransfer" << endl;
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

void RScalarHierarchy :: MakeSysMatrix()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeSysMatrix" << endl;
#endif

  dir_sys = 0;

  if (node[level].auxmat != node[level].sysmat)
    {
      node[level-1].sysmat = new RScalarMatrix(csize_sys, cnne_sys, dir_sys); 
    }
  else
    {
      node[level-1].sysmat = node[level-1].auxmat;
    }
}

void RScalarHierarchy :: MakeAuxMatrix()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeAuxMatrix" << endl;
#endif

  dir_aux = 0;

  node[level-1].auxmat = new RScalarMatrix(csize_aux, cnne_aux, dir_aux); 
}

void RScalarHierarchy :: MakeVector()
{
#ifdef TRACE
  (*trace) << "entering RScalarHierarchy::MakeVector" << endl;
#endif

  node[level].vec1 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec2 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec3 = new RealVector(fsize_sys, numrhs, dof); 
  node[level].vec4 = new RealVector(fsize_sys, numrhs, dof); 
}


///////////////////////////////////////// complex ///////////////////////////////////////

CScalarHierarchy :: CScalarHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs)
  : BaseHierarchy(param, asize, aoffset, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::CScalarHierarchy" << endl;
#endif

}
  
CScalarHierarchy :: ~CScalarHierarchy()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::~CScalarHierarchy" << endl;
#endif
  
}

void CScalarHierarchy :: Init()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::Init" << endl;
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

void CScalarHierarchy :: DeleteAuxMemory()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::DeleteAuxMemory" << endl;
#endif
  
  delete node[level].coarse;
  node[level].coarse = NULL;

  delete node[level].auxmat;
  node[level].auxmat = NULL;
}

void CScalarHierarchy :: MakeSmoother()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeSmoother" << endl;
#endif

  switch (param->GetSmoothingType())
    {
    case 1:
      if (param->GetSmoothingStepFor() == 1)
	{
	  node[level].smooth = new CScalarFastGSSmoother(*param,fsize_sys,numrhs);
	}
      else
	{
	  node[level].smooth = new CScalarGSSmoother(*param,fsize_sys,numrhs);
	}
      break;

    case 2:
      //node[level].smooth = new CScalarFastGSSmoother(*param,fsize_sys,numrhs);
      ;
      break;

    default:
      node[level].smooth = new CScalarGSSmoother(*param,fsize_sys,numrhs);
    }
}

void CScalarHierarchy :: MakeCoarsening()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeCoarsening" << endl;
#endif

  switch (param->GetCoarsening())
    {
    case 1:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
      break;

    case 2:
      //node[level].coarse = new CScalarCoarse();
      ;
      break;

    default:
      node[level].coarse = new StandardCoarse(node[level].auxmat, param->GetMaxNeighbour());
    }
}

void CScalarHierarchy :: MakeSysTransfer()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeSysTransfer" << endl;
#endif

  switch (param->GetTransfer())
    {
      case 1:
	node[level].syspro = node[level].auxpro;//new CScalarTransfer(fsize_sys, rsw);
	node[level].sysres = node[level].syspro;
      break;

    case 2:
      //node[level].syspro = new CScalarProlongation();
      //node[level].sysres = node[level].syspro;
      ;
      break;

    default:
      node[level].syspro = node[level].auxpro;//new CScalarTransfer(fsize_sys, rsw);
      node[level].sysres = node[level].syspro;
    }
}

void CScalarHierarchy :: MakeAuxTransfer()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeAuxTransfer" << endl;
#endif

  switch (param->GetTransfer())
    {
      case 1:
	node[level].auxpro = new RScalarTransfer(fsize_sys, rsw);
	node[level].auxres = node[level].auxpro;
      break;

    case 2:
      //node[level].auxpro = new CScalarProlongation();
      //node[level].auxres = node[level].auxpro;
      ;
      break;

    default:
      node[level].auxpro = new RScalarTransfer(fsize_sys, rsw);
      node[level].auxres = node[level].auxpro;
    }
}

void CScalarHierarchy :: MakeSysMatrix()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeSysMatrix" << endl;
#endif

  dir_sys = 0;

  node[level-1].sysmat = new CScalarMatrix(csize_sys, cnne_sys, dir_sys); 
}

void CScalarHierarchy :: MakeAuxMatrix()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeAuxMatrix" << endl;
#endif

  dir_aux = 0;

  node[level-1].auxmat = new RScalarMatrix(csize_aux, cnne_aux, dir_aux); 
}

void CScalarHierarchy :: MakeVector()
{
#ifdef TRACE
  (*trace) << "entering CScalarHierarchy::MakeVector" << endl;
#endif

  node[level].vec1 = new ComplexVector(fsize_sys, numrhs, dof); 
  node[level].vec2 = new ComplexVector(fsize_sys, numrhs, dof); 
  node[level].vec3 = new ComplexVector(fsize_sys, numrhs, dof); 
  node[level].vec4 = new ComplexVector(fsize_sys, numrhs, dof); 
}

}
