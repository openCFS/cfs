#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseMGPrecond :: BaseMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::RScalarMGPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  offset = 15;
  dof    = adof;

  Integer i;

  switch (param.GetCycle())
    {
    case 0:
      for (i=0; i<20; i++)
	{
	  cycle[i] = 0;
	}

      break;

    case 1:
      for (i=0; i<20; i++)
	{
	  cycle[i] = 1;
	}
      
      break;

    case 2:
      for (i=0; i<20; i++)
	{
	  cycle[i] = 2;
	}

      break;

    default:
      for (i=0; i<20; i++)
	{
	  cycle[i] = 1;
	}
    }
}

BaseMGPrecond :: ~BaseMGPrecond()
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::~RScalarMGPrecond" << endl;
#endif
}

void BaseMGPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::Step" << endl;
#endif
  
  Integer i;

  HierarchyNode * fnode, *cnode;

  fnode = mathier->GetNode(offset-level);
  cnode = mathier->GetNode(offset-level-1);

  if (level+1 == mathier->GetSize())
    {
      fnode->sysmat->Solve(rhs, sol);
      return;
    }
  else
    {
      fnode->smooth->StepForward(*fnode->sysmat, *fnode->auxmat, rhs, sol, *fnode->vec3, *fnode->vec4,level);
      cnode->vec2->Init();
      cnode->vec1->Init();
      fnode->sysres->MulthH(*fnode->vec3, *cnode->vec2);

      for (i=1; i<=cycle[level]; i++)
	{
	  Step(*cnode->vec2, *cnode->vec1, level+1);
	}

      fnode->syspro->MultHh(*cnode->vec1, sol);
      fnode->smooth->StepBackward(*fnode->sysmat, *fnode->auxmat, rhs, sol,level);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

RScalarMGPrecond :: RScalarMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BaseMGPrecond(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::RScalarMGPrecond" << endl;
#endif

  mathier = new RScalarHierarchy(param, size, offset, numrhs);
}
  
RScalarMGPrecond :: ~RScalarMGPrecond()
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::~RScalarMGPrecond" << endl;
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////

RBlockMGPrecond :: RBlockMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof)
  : BaseMGPrecond(param, asize, anumrhs, adof)
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::RScalarMGPrecond" << endl;
#endif

  mathier = new RBlockHierarchy(param, size, offset, numrhs, dof);
}
  
RBlockMGPrecond :: ~RBlockMGPrecond()
{
#ifdef TRACE
  (*trace) << "entering RScalarMGPrecond::~RScalarMGPrecond" << endl;
#endif
}

/////////////////////////////////////////////////////////////////////////////////////////

RFullMGPrecond :: RFullMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BaseMGPrecond(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RFullMGPrecond::RFullMGPrecond" << endl;
#endif

  mathier = new RFullHierarchy(param, size, offset, numrhs);
}
  
RFullMGPrecond :: ~RFullMGPrecond()
{
#ifdef TRACE
  (*trace) << "entering RFullMGPrecond::~RFullMGPrecond" << endl;
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////

CScalarMGPrecond :: CScalarMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BaseMGPrecond(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering CScalarMGPrecond::CScalarMGPrecond" << endl;
#endif

  mathier = new CScalarHierarchy(param, size, offset, numrhs);
}
  
CScalarMGPrecond :: ~CScalarMGPrecond()
{
#ifdef TRACE
  (*trace) << "entering CScalarMGPrecond::~CScalarMGPrecond" << endl;
#endif
}
}
