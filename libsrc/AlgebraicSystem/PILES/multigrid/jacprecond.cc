#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RJACScalarPrecond :: RJACScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering RJACScalarPrecond::RJACScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  omega  = param.GetDampIter();
}
  
RJACScalarPrecond :: ~RJACScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering RJACScalarPrecond::~RJACScalarPrecond" << endl;
#endif
}

void RJACScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RJACScalarPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i;

  for (i=1; i<=size; i++)
    {
      b.Elem(i) = omega*a.Get(i)/mat->GetDiag(i);
    }
}

void RJACScalarPrecond :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering RJACScalarPrecond::Calc" << endl;
#endif

  RScalarMatrix & a = (RScalarMatrix &) sysmat;

  mat = &a;
}



CJACScalarPrecond :: CJACScalarPrecond(Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering CJACScalarPrecond::CJACScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;

  val = new Double[2*size*numrhs];
}
  
CJACScalarPrecond :: ~CJACScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering CJACScalarPrecond::~CJACScalarPrecond" << endl;
#endif

  delete [] val;
}

void CJACScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CJACScalarPrecond::Step" << endl;
#endif
  
  ComplexVector & a = (ComplexVector &) rhs;
  ComplexVector & b = (ComplexVector &) sol;

  Integer i,j;

  j = 1;

  for (i=0; i<size; i++)
    {
      b.Elem(j)   = 1.5*(val[j-1]*a.Get(j) - val[j]*a.Get(j+1));
      b.Elem(j+1) = 1.5*(val[j]*a.Get(j) + val[j-1]*a.Get(j+1));

      j += 2;
    }
}

void CJACScalarPrecond :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering CJACScalarPrecond::Calc" << endl;
#endif

  CScalarMatrix & a = (CScalarMatrix &) sysmat;

  mat = &a;

  Integer i,j;
  Double * d, h;

  j = 1;

  for (i=0; i<size; i++)
    {
      d = mat->GetDiag(i);
      h = d[0]*d[0]+d[1]*d[1];

      val[j-1] = d[0]/h;
      val[j]   = -d[1]/h;

      j += 2;
    }
}
}
