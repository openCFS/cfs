#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RFullGSSmoother :: RFullGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : BaseSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RFullGSSmoother::RFullGSSmoother" << endl;
#endif

}
  
RFullGSSmoother :: ~RFullGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RFullGSSmoother::~RFullGSSmoother" << endl;
#endif
}

void RFullGSSmoother :: Calc(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RFullGSSmoother::Calc" << endl;
#endif
}

void RFullGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
				    BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RFullGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k;
  Double sum;

  RFullMatrix & mat = (RFullMatrix &) sys;
  RealVector & f    = (RealVector &) rhs;
  RealVector & u    = (RealVector &) sol;
  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=1; i<=numsmoothfor; i++)
    {
      for(j=1; j<=size; j++)
	{	
	  sum = 0;
      
	  for (k=1; k<=size; k++)
	    {
	      sum += mat.Get(j,k)*u.Get(k);
	    }
	  
	  sum = (f.Get(j)-sum)/mat.GetDiag(j);

	  u.Elem(j) += sum;
	}
    }

  sys.Mult(sol,res);
  def.Add(pone, rhs, mone, res);
}

void RFullGSSmoother :: StepBackward(BaseMatrix & sys, BaseMatrix & aux, 
				       BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RFullGSSmoother::StepBackward" << endl;
#endif

  Integer i,j,k;
  Double sum;

  RFullMatrix & mat = (RFullMatrix &) sys;
  RealVector & f      = (RealVector &) rhs;
  RealVector & u      = (RealVector &) sol;

  for (i=1; i<=numsmoothback; i++)
    {
      for(j=size; j>= 1; j--)
	{	
	  sum = 0;
	
	  for (k=1; k<=size; k++)
	    {
	      sum += mat.Get(j,k)*u.Get(k);
	    }
	  
	  sum = (f.Get(j)-sum)/mat.GetDiag(j);
	  u.Elem(j) += sum;
	}
    }
}

RFullFastGSSmoother :: RFullFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : RFullGSSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RFullFastGSSmoother::RFullFastGSSmoother" << endl;
#endif

}
  
RFullFastGSSmoother :: ~RFullFastGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RFullFastGSSmoother::~RFullFastGSSmoother" << endl;
#endif
}

void RFullFastGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
					  BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RFullFastGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k;
  Double sum;

  RFullMatrix & mat = (RFullMatrix &) sys;
  RealVector & f    = (RealVector &) rhs;
  RealVector & u    = (RealVector &) sol;
  RealVector & d    = (RealVector &) def;

  for (i=1; i<=size; i++)
    {
      d.Elem(i) = f.Get(i);
    }

  for (i=1; i<=size; i++)
    {
      sum = d.Get(i)/mat.GetDiag(i);
      u.Elem(i) = sum;

      for (j=1; j<=size; j++)
	{
	  d.Elem(j) -= (mat.Get(i,j)*sum);
	}
    }
}

}
