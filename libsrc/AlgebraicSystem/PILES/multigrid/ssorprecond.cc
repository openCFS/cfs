#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RSSORScalarPrecond :: RSSORScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering RSSORScalarPrecond::RSSORScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  omega  = param.GetDampIter();
}
  
RSSORScalarPrecond :: ~RSSORScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering RSSORScalarPrecond::~RSSORScalarPrecond" << endl;
#endif

}

void RSSORScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RSSORScalarPrecond::Step" << endl;
#endif
  
  RealVector & a = (RealVector &) rhs;
  RealVector & b = (RealVector &) sol;

  Integer i,j,rs,pos;
  Double sum;

  for(i=1; i<=size; i++)
    {	
      sum = 0;
      rs  = mat->GetRowSize(i);

      for (j=1; j<=rs; j++)
	{
	  pos  = mat->GetMatrixPos(i,j);
	  sum += mat->Get(i,j)*b.Get(pos);
	}
      
      sum = omega*(-sum+a.Get(i))/mat->GetDiag(i);
      
      b.Elem(i) += sum;
    }

  for(i=size; i>=1; i--)
    {	
      sum = 0;
      rs  = mat->GetRowSize(i);
      
      for (j=1; j<=rs; j++)
	{
	  pos  = mat->GetMatrixPos(i,j);
	  sum += mat->Get(i,j)*b.Get(pos);
	}
      
      sum = omega*(-sum+a.Get(i))/mat->GetDiag(i);
      b.Elem(i) += sum;
    }
}

void RSSORScalarPrecond :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering RSSORScalarPrecond::Calc" << endl;
#endif
  
  RScalarMatrix & c = (RScalarMatrix &) sysmat;
  mat = &c;
}




CSSORScalarPrecond :: CSSORScalarPrecond(Integer asize, Integer anumrhs)
  : BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering CSSORScalarPrecond::CSSORScalarPrecond" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  dof    = 1;
}
  
CSSORScalarPrecond :: ~CSSORScalarPrecond()
{
#ifdef TRACE
  (*trace) << "entering CSSORScalarPrecond::~CSSORScalarPrecond" << endl;
#endif

}

void CSSORScalarPrecond :: Step(BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CSSORScalarPrecond::Step" << endl;
#endif
  
  ComplexVector & a = (ComplexVector &) rhs;
  ComplexVector & b = (ComplexVector &) sol;

  Integer i,j,k,rs,pos;
  Double re,im, * h, d, p, q;

  k = 1;

  for(i=1; i<=size; i++)
    {	
      re = 0;
      im = 0;

      rs  = mat->GetRowSize(i);

      for (j=1; j<=rs; j++)
	{
	  pos  = 2*mat->GetMatrixPos(i,j)-1;
	  h    = mat->Get(i,j);

	  re += (h[0]*b.Get(pos) - h[1]*b.Get(pos+1));
	  im += (h[0]*b.Get(pos+1) + h[1]*b.Get(pos));
	}

      h = mat->GetDiag(i);
      d = h[0]*h[0]+h[1]*h[1];

      re += (-a.Get(k));
      im += (-a.Get(k+1));

      p = (re*h[0] + im*h[1])/d;
      q = (im*h[0] - re*h[1])/d;

      b.Elem(k)   += p;
      b.Elem(k+1) += q;

      k += 2;
    }

  k = 2*size-1;

  for(i=size; i>=1; i--)
    {	
      re = 0;
      im = 0;

      rs  = mat->GetRowSize(i);

      for (j=1; j<=rs; j++)
	{
	  pos  = 2*mat->GetMatrixPos(i,j)-1;
	  h    = mat->Get(i,j);

	  re += (h[0]*b.Get(pos) - h[1]*b.Get(pos+1));
	  im += (h[0]*b.Get(pos+1) + h[1]*b.Get(pos));
	}

      h = mat->GetDiag(i);
      d = h[0]*h[0]+h[1]*h[1];

      re += (-a.Get(k));
      im += (-a.Get(k+1));

      p = (re*h[0] + im*h[1])/d;
      q = (im*h[0] - re*h[1])/d;

      b.Elem(k)   += p;
      b.Elem(k+1) += q;

      k -= 2;
    }
}

void CSSORScalarPrecond :: Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
{
#ifdef TRACE
  (*trace) << "entering CSSORScalarPrecond::Calc" << endl;
#endif
  
  CScalarMatrix & c = (CScalarMatrix &) sysmat;
  mat = &c;
}
}
