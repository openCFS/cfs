#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RIDScalarPrecond :: RIDScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering RIDScalarPrecond::RIDScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  omega  = param.GetDampIter();
}
  
RIDScalarPrecond :: ~RIDScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering RIDScalarPrecond::~RIDScalarPrecond" << endl;
#endif
}

void RIDScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RIDScalarPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i;

  for (i=1; i<=size; i++)
    {
      b.Elem(i) = omega*a.Get(i);
    }
}

void RIDScalarPrecond :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering RIDScalarPrecond::Calc" << endl;
#endif

}


CIDScalarPrecond :: CIDScalarPrecond(Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering CIDScalarPrecond::CIDScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
}
  
CIDScalarPrecond :: ~CIDScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering CIDScalarPrecond::~CIDScalarPrecond" << endl;
#endif
}

void CIDScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CIDScalarPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i,j;
  
  j = 1;

  for (i=1; i<=size; i++)
    {
      b.Elem(j)   = a.Get(j);
      b.Elem(j+1) = a.Get(j+1);

      j += 2;
    }
}

void CIDScalarPrecond :: Calc(BaseMatrix & sysmat)
{
#ifdef TRACE
  (*trace) << "entering CIDScalarPrecond::Calc" << endl;
#endif
}



RIDBlockPrecond :: RIDBlockPrecond(Integer asize, Integer adof)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering RIDBlockPrecond::RIDBlockPrecond" << endl;
#endif

  size = asize;
  dof  = adof;
}
  
RIDBlockPrecond :: ~RIDBlockPrecond()
{
#ifdef TRACE
  (*trace) << "entering RIDBlockPrecond::~RIDBlockPrecond" << endl;
#endif
}

void RIDBlockPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RIDBlockPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i;

  for (i=1; i<=size; i++)
    {
      b.Elem(i) = a.Get(i);
    }
}

void RIDBlockPrecond :: Calc(BaseMatrix & sysmat)
{
#ifdef TRACE
  (*trace) << "entering RIDBlockPrecond::Calc" << endl;
#endif
}


CIDBlockPrecond :: CIDBlockPrecond(Integer asize, Integer adof)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering CIDBlockPrecond::CIDBlockPrecond" << endl;
#endif

  size = asize;
  dof  = adof;
}
  
CIDBlockPrecond :: ~CIDBlockPrecond()
{
#ifdef TRACE
  (*trace) << "entering CIDBlockPrecond::~CIDBlockPrecond" << endl;
#endif
}

void CIDBlockPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CIDBlockPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i,j;
  
  j = 1;

  for (i=1; i<=size; i++)
    {
      b.Elem(j)   = a.Get(j);
      b.Elem(j+1) = a.Get(j+1);

      j += 2;
    }
}

void CIDBlockPrecond :: Calc(BaseMatrix & sysmat)
{
#ifdef TRACE
  (*trace) << "entering CIDBlockPrecond::Calc" << endl;
#endif
}

}
