#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RScalarGSSmoother :: RScalarGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : BaseSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RScalarGSSmoother::RScalarGSSmoother" << endl;
#endif

}
  
RScalarGSSmoother :: ~RScalarGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RScalarGSSmoother::~RScalarGSSmoother" << endl;
#endif
}

void RScalarGSSmoother :: Calc(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RScalarGSSmoother::Calc" << endl;
#endif
}

void RScalarGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
				      BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RScalarGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  RScalarMatrix & mat = (RScalarMatrix &) sys;
  RealVector & f      = (RealVector &) rhs;
  RealVector & u      = (RealVector &) sol;
  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=1; i<=numsmoothfor; i++)
    {
      for(j=1; j<=size; j++)
	{	
	  sum = 0;
	  rs  = mat.GetRowSize(j);
	  
	  for (k=1; k<=rs; k++)
	    {
	      pos  = mat.GetMatrixPos(j,k);
	      sum += mat.Get(j,k)*u.Get(pos);
	    }
	  
	  sum = (f.Get(j)-sum)/mat.GetDiag(j);

	  u.Elem(j) += sum;
	}
    }

  sys.Mult(sol,res);
  def.Add(pone, rhs, mone, res);
}

void RScalarGSSmoother :: StepBackward(BaseMatrix & sys, BaseMatrix & aux, 
				       BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RScalarGSSmoother::StepBackward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  RScalarMatrix & mat = (RScalarMatrix &) sys;
  RealVector & f      = (RealVector &) rhs;
  RealVector & u      = (RealVector &) sol;

  for (i=1; i<=numsmoothback; i++)
    {
      for(j=size; j>= 1; j--)
	{	
	  sum = 0;
	  rs  = mat.GetRowSize(j);
	  
	  for (k=1; k<=rs; k++)
	    {
	      pos  = mat.GetMatrixPos(j,k);
	      sum += mat.Get(j,k)*u.Get(pos);
	    }
	  
	  sum = (f.Get(j)-sum)/mat.GetDiag(j);
	  u.Elem(j) += sum;
	}
    }
}

RScalarFastGSSmoother :: RScalarFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : RScalarGSSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering RScalarFastGSSmoother::RScalarFastGSSmoother" << endl;
#endif

}
  
RScalarFastGSSmoother :: ~RScalarFastGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RScalarFastGSSmoother::~RScalarFastGSSmoother" << endl;
#endif
}

void RScalarFastGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
					  BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RScalarFastGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  RScalarMatrix & mat = (RScalarMatrix &) sys;
  RealVector & f      = (RealVector &) rhs;
  RealVector & u      = (RealVector &) sol;
  RealVector & d      = (RealVector &) def;

  for (i=1; i<=size; i++)
    {
      d.Elem(i) = f.Get(i);
    }

  for (i=1; i<=size; i++)
    {
      sum = d.Get(i)/mat.GetDiag(i);
      u.Elem(i) = sum;

      rs  = mat.GetRowSize(i);

      for (j=1; j<=rs; j++)
	{
	  pos = mat.GetMatrixPos(i,j);
	  d.Elem(pos) -= (mat.Get(i,j)*sum);
	}
    }
}

//////////////////////////////////// complex //////////////////////////////////////////////

CScalarGSSmoother :: CScalarGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : BaseSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering CScalarGSSmoother::CScalarGSSmoother" << endl;
#endif

}
  
CScalarGSSmoother :: ~CScalarGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering CScalarGSSmoother::~CScalarGSSmoother" << endl;
#endif
}

void CScalarGSSmoother :: Calc(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering CScalarGSSmoother::Calc" << endl;
#endif
}

void CScalarGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
				      BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CScalarGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  CScalarMatrix & mat = (CScalarMatrix &) sys;
  ComplexVector & f   = (ComplexVector &) rhs;
  ComplexVector & u   = (ComplexVector &) sol;

  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=1; i<=numsmoothfor; i++)
    {
//       for(j=1; j<=size; j++)
// 	{	
// 	  sum = 0;
// 	  rs  = mat.GetRowSize(j);
	  
// 	  for (k=1; k<=rs; k++)
// 	    {
// 	      pos  = mat.GetMatrixPos(j,k);
// 	      //sum += mat.Get(j,k)*u.Get(pos);
// 	    }
	  
// 	  sum = (f.Get(j)-sum)/mat.GetDiag(j);

// 	  u.Elem(j) += sum;
// 	}
    }

  sys.Mult(sol,res);
  def.Add(pone, rhs, mone, res);
}

void CScalarGSSmoother :: StepBackward(BaseMatrix & sys, BaseMatrix & aux, 
				       BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CScalarGSSmoother::StepBackward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  CScalarMatrix & mat = (CScalarMatrix &) sys;
  ComplexVector & f   = (ComplexVector &) rhs;
  ComplexVector & u   = (ComplexVector &) sol;

  for (i=1; i<=numsmoothback; i++)
    {
//       for(j=size; j>= 1; j--)
// 	{	
// 	  sum = 0;
// 	  rs  = mat.GetRowSize(j);
	  
// 	  for (k=1; k<=rs; k++)
// 	    {
// 	      pos  = mat.GetMatrixPos(j,k);
// 	      sum += mat.Get(j,k)*u.Get(pos);
// 	    }
	  
// 	  sum = (f.Get(j)-sum)/mat.GetDiag(j);
// 	  u.Elem(j) += sum;
//	}
    }
}

CScalarFastGSSmoother :: CScalarFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs) 
  : CScalarGSSmoother(param, asize, anumrhs)
{
#ifdef TRACE
  (*trace) << "entering CScalarFastGSSmoother::CScalarFastGSSmoother" << endl;
#endif

}
  
CScalarFastGSSmoother :: ~CScalarFastGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering CScalarFastGSSmoother::~CScalarFastGSSmoother" << endl;
#endif
}

void CScalarFastGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
					  BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering CScalarFastGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k,pos,rs;
  Double sum;

  CScalarMatrix & mat = (CScalarMatrix &) sys;
  ComplexVector & f   = (ComplexVector &) rhs;
  ComplexVector & u   = (ComplexVector &) sol;
  ComplexVector & d   = (ComplexVector &) def;

  for (i=1; i<=size; i++)
    {
      d.Elem(i) = f.Get(i);
    }

  for (i=1; i<=size; i++)
    {
//       sum = d.Get(i)/mat.GetDiag(i);
//       u.Elem(i) = sum;

//       rs  = mat.GetRowSize(i);

//       for (j=1; j<=rs; j++)
// 	{
// 	  pos = mat.GetMatrixPos(i,j);
// 	  d.Elem(pos) -= (mat.Get(i,j)*sum);
// 	}
    }
}

}
