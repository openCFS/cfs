#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RBlockGSSmoother :: RBlockGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof) 
  : BaseSmoother(param, asize, anumrhs, adof)
{
#ifdef TRACE
  (*trace) << "entering RBlockGSSmoother::RBlockGSSmoother" << endl;
#endif

  sum = new Double[dof];
  x   = new DenseMatrix(dof,dof);
}
  
RBlockGSSmoother :: ~RBlockGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RBlockGSSmoother::~RBlockGSSmoother" << endl;
#endif

  if (sum != NULL)
    {
      delete [] sum;
    }
}

void RBlockGSSmoother :: Calc(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RBlockGSSmoother::Calc" << endl;
#endif

  RBlockMatrix & mat = (RBlockMatrix &) amat;

  mat.SetDiagInv();
}

void RBlockGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
				     BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RBlockGSSmoother::StepForward" << endl;
#endif

  Integer i,j,k,l,m,n,p,pos,rs;
  Double v;

  RBlockMatrix & mat = (RBlockMatrix &) sys;
  RealVector & f     = (RealVector &) rhs;
  RealVector & u     = (RealVector &) sol;

  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  pone = +1;
  mone = -1;

  for (i=1; i<=numsmoothfor; i++)
    {
      for(j=1; j<=size; j++)
	{	
	  rs = mat.GetRowSize(j);
	  p  = (j-1)*dof+1;

	  for (k=0; k<dof; k++)
	    {
	      sum[k] = 0;
	    }
	  
	  for (k=1; k<=rs; k++)
	    {
	      x    = mat.Get(j,k);
	      pos  = mat.GetMatrixPos(j,k);
	      n    = (pos-1)*dof+1;

	      for (l=0; l<dof; l++)
		{
		  for (m=0; m<dof; m++)
		    {
		      sum[l] += x->Get(l+1,m+1)*u.Get(n+m);
		    }
		}
	    }

	  for (k=0; k<dof; k++)
	    {
	      sum[k] = f.Get(p+k) - sum[k];
	    }

	  x = mat.GetDiagInv(j);

	  for (k=0; k<dof; k++)
	    {
	      v = 0;

	      for (l=0; l<dof; l++)
		{
		  v += x->Get(k+1,l+1)*sum[l];
		}

	      u.Elem(p+k) += v;
	    }
	}
    }

  sys.Mult(sol,res);
  def.Add(pone, rhs, mone, res);
}

void RBlockGSSmoother :: StepBackward(BaseMatrix & sys, BaseMatrix & aux, 
				      BaseVector & rhs, BaseVector & sol, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RBlockGSSmoother::StepBackward" << endl;
#endif

  Integer i,j,k,l,m,n,p,pos,rs;
  Double v;

  RBlockMatrix & mat = (RBlockMatrix &) sys;
  RealVector & f     = (RealVector &) rhs;
  RealVector & u     = (RealVector &) sol;

  for (i=1; i<=numsmoothback; i++)
    {
      for(j=size; j>= 1; j--)
	{	
	  rs = mat.GetRowSize(j);
	  p  = (j-1)*dof+1;

	  for (k=0; k<dof; k++)
	    {
	      sum[k] = 0;
	    }
	  
	  for (k=1; k<=rs; k++)
	    {
	      x    = mat.Get(j,k);
	      pos  = mat.GetMatrixPos(j,k);
	      n    = (pos-1)*dof+1;

	      for (l=0; l<dof; l++)
		{
		  for (m=0; m<dof; m++)
		    {
		      sum[l] += x->Get(l+1,m+1)*u.Get(n+m);
		    }
		}
	    }

	  for (k=0; k<dof; k++)
	    {
	      sum[k] = f.Get(p+k) - sum[k];
	    }

	  x = mat.GetDiagInv(j);

	  for (k=0; k<dof; k++)
	    {
	      v = 0;

	      for (l=0; l<dof; l++)
		{
		  v += x->Get(k+1,l+1)*sum[l];
		}

	      u.Elem(p+k) += v;
	    }
	}
    }
}


RBlockFastGSSmoother :: RBlockFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof) 
  : RBlockGSSmoother(param, asize, anumrhs, adof)
{
#ifdef TRACE
  (*trace) << "entering RBlockFastGSSmoother::RBlockFastGSSmoother" << endl;
#endif

}
  
RBlockFastGSSmoother :: ~RBlockFastGSSmoother()
{
#ifdef TRACE
  (*trace) << "entering RBlockFastGSSmoother::~RBlockFastGSSmoother" << endl;
#endif

}

void RBlockFastGSSmoother :: StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
					 BaseVector & def, BaseVector & res, Integer level)
{
#ifdef TRACE
  (*trace) << "entering RBlockFastGSSmoother::StepForward" << endl;
#endif

  Integer j,k,l,m,p,pos,rs;

  RBlockMatrix & mat = (RBlockMatrix &) sys;
  RealVector & f     = (RealVector &) rhs;
  RealVector & u     = (RealVector &) sol;
  RealVector & d     = (RealVector &) def;

  for (j=1; j<=size*dof; j++)
    {
      d.Elem(j) = f.Get(j);
    }

//   for(j=1; j<=size; j++)
//     {	
//       x = mat.GetDiagInv(j);
//       p = (j-1)*dof+1;
      
//       for (k=0; k<dof; k++)
// 	{
// 	  sum[k] = 0;
	  
// 	  for (l=0; l<dof; l++)
// 	    {
// 	      sum[k] += x->Get(k+1,l+1)*d.Get(p+l);
// 	    }
	  
// 	  u.Elem(p+k) = sum[k];
// 	}

//       rs = mat.GetRowSize(j);
  
//       for (k=1; k<=rs; k++)
// 	{
// 	  x    = mat.Get(j,k);
// 	  pos  = mat.GetMatrixPos(j,k);
// 	  p    = (pos-1)*dof+1;
	  
// 	  for (l=0; l<dof; l++)
// 	    {
// 	      for (m=0; m<dof; m++)
// 		{
// 		  d.Elem(p+l) -= x->Get(m+1,l+1)*sum[m];
// 		}
// 	    }
// 	}    
//     }

  Integer n;
  Double v;

  for(j=1; j<=size; j++)
    {	
      rs = mat.GetRowSize(j);
      p  = (j-1)*dof+1;
      
      for (k=0; k<dof; k++)
	{
	  sum[k] = 0;
	}
      
      for (k=2; k<=rs; k++)
	{
	  x    = mat.Get(j,k);
	  pos  = mat.GetMatrixPos(j,k);
	  n    = (pos-1)*dof+1;
	  
	  for (l=0; l<dof; l++)
	    {
	      for (m=0; m<dof; m++)
		{
		  sum[l] += x->Get(l+1,m+1)*u.Get(n+m);
		}
	    }
	}
      
      for (k=0; k<dof; k++)
	{
	  sum[k] = f.Get(p+k) - sum[k];
	}
      
      x = mat.GetDiagInv(j);
      
      for (k=0; k<dof; k++)
	{
	  v = 0;
	  
	  for (l=0; l<dof; l++)
	    {
	      v += x->Get(k+1,l+1)*sum[l];
	    }
	  
	  u.Elem(p+k) = v;
	}

      for (k=1; k<=rs; k++)
	{
	  x    = mat.Get(j,k);
	  pos  = mat.GetMatrixPos(j,k);
	  n    = (pos-1)*dof+1;
	  
	  for (l=0; l<dof; l++)
	    {
	      for (m=0; m<dof; m++)
		{
		  d.Elem(n+l) -= x->Get(m+1,l+1)*u.Get(p+m);
		}
	    }
	}    
    }
}


}
