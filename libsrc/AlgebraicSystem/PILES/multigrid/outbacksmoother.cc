#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

OutbackSmoother :: OutbackSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : BaseSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering OutbackSmoother::OutbackSmoother" << endl;
#endif

}
  
OutbackSmoother :: ~OutbackSmoother()
{
#ifdef TRACE
  (*trace) << "entering OutbackSmoother::~OutbackSmoother" << endl;
#endif
}

void OutbackSmoother :: Calc(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering OutbackSmoother::Calc" << endl;
#endif
}

void OutbackSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
				    BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering OutbackSmoother::StepForward" << endl;
#endif

  Integer i,j,k;
  Double sum;
  Double damp = 0.7;

  OutbackMatrix & mat = (OutbackMatrix &) sys;
  RealVector & f    = (RealVector &) rhs;
  RealVector & u    = (RealVector &) sol;
  RealVector * des = new RealVector(size, 1, 1);
  RealVector & r = (RealVector &) *des;
  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=0; i<4; i++)
    {
      sys.Mult(u, r, 1);
      
      for (j=1; j<=size; j++)
	{
	  u.Elem(j) += (damp*(f.Get(j) - r.Get(j))/mat.GetDiag(j));
	}
    }

  sys.Mult(sol,res);
  def.Add(pone, rhs, mone, res);
}

void OutbackSmoother :: StepBackward(BaseMatrix & sys, BaseMatrix & aux, 
					BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering OutbackSmoother::StepBackward" << endl;
#endif

  Integer i,j,k;
  Double sum;
  Double damp = 0.7;

  OutbackMatrix & mat = (OutbackMatrix &) sys;
  RealVector & f      = (RealVector &) rhs;
  RealVector & u      = (RealVector &) sol;
  RealVector * res = new RealVector(size, 1, 1);
  RealVector & r = (RealVector &) *res;
  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=0; i<4; i++)
    {
      sys.Mult(u, r, 1);
      
      for (j=1; j<=size; j++)
	{
	  u.Elem(j) += (damp*(f.Get(j) - r.Get(j))/mat.GetDiag(j));
	}
    }
}
}
