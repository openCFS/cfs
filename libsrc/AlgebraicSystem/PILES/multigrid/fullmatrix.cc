#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RFullMatrix :: RFullMatrix(Integer asize, Integer adir, Boolean mem)
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::RFullMatrix" << endl;
#endif

  size   = asize;
  dir    = adir;
  nne    = size*size;
  matfac = 1;
  
  pos     = NULL;
  start   = NULL;
  cm      = NULL;
  f       = NULL;
  perm_c  = NULL;
  perm_r  = NULL;
  fac     = NULL;

  if (mem)
    {
      val = new Double[nne];

      Integer i;

      for (i=0; i<nne; i++)
	{
	  val[i] = 0;
	}

      outback = FALSE;
    }
  else
    {
      val     = NULL;
      outback = TRUE;
    }

  if (dir != 0 && mem)
    {
      valdir = new Double[dir];
      numdir = new Integer[dir];
      compdir= new Integer[dir];
    }
  else
    {
      valdir  = NULL;
      numdir  = NULL;
      compdir = NULL;
    }

  calculated = FALSE;
  setgraph   = FALSE;
  buildindir = FALSE;
}
  
RFullMatrix :: ~RFullMatrix()
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::~RFullMatrix" << endl;
#endif

  if (val != NULL & !outback)
    {
      delete [] val;
      delete [] valdir;
      delete [] numdir;
    }
  
  if (fac != NULL)
    {
      delete [] fac;
      delete [] y;
    }
}

void RFullMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Mult" << endl;
#endif

  RealVector & u = (RealVector &) vec1;
  RealVector & v = (RealVector &) vec2;

  Integer i,j;
  Double sum;

  for (i=0; i<size; i++)
    {
      sum = 0;

      for (j=0; j<size; j++)
	{
	  sum += val[i*size+j]*u.Get(j+1);
	}

      v.Elem(i+1) = factor*sum;
    }
}

void RFullMatrix :: BuildInDirichlet()
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::BuildInDirichlet" << endl;
#endif
  
  if (buildindir)
    {
      return;
    }

  Integer i;
  Double maxdiag;
  
  maxdiag = val[0];

  for (i=1; i<size; i++)
    {
      if (maxdiag < val[i*size+i])
	{
	  maxdiag = val[i*size+i];
	}
    }

  maxdiag += 1e8;

  for (i=0; i<dir; i++)
    {
      val[(numdir[i]-1)*(size+1)] += maxdiag;
    }

  buildindir = TRUE;
}

void RFullMatrix :: BuildInDirichlet(BaseVector & rhs)
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::BuildInDirichlet" << endl;
#endif

  if (buildindir)
    {
      return;
    }

  RealVector & a = (RealVector &) rhs;

  Integer i;
  Double maxdiag;
  
  maxdiag = val[0];

  for (i=1; i<size; i++)
    {
      if (maxdiag < val[i*size+i])
	{
	  maxdiag = val[i*size+i];
	}
    }

  maxdiag += 1e8;

  for (i=0; i<dir; i++)
    {
      val[(numdir[i]-1)*(size+1)] += maxdiag;
      a.Elem(numdir[i]) += maxdiag*valdir[i];
    }

  buildindir = TRUE;
}

void RFullMatrix :: Factor()
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Factor" << endl;
#endif

  Integer i,j,k;
  Double h;
  
  fac = new Double[size*size];
  y   = new Double[size];

  /// factorization
  
  for (i=0; i<size; i++)
    {
      for (j=i; j<size; j++)
	{
	  h = 0;
	  
	  for (k=0; k<=i-1; k++)
	    {
	      h += fac[i*size+k]*fac[k*size+j];
	    }
	  
	  fac[i*size+j] = val[i*size+j] - h;
	}
      
      for (j=i+1; j<size; j++)
	{
	  h = 0;
	  
	  for (k=0; k<=i-1; k++)
	    {
	      h += fac[j*size+k]*fac[k*size+i];
	    }
	  
	  fac[j*size+i] = (val[i*size+j] - h)/fac[i*size+i];
	}
    }
}

void RFullMatrix :: Solve(BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Solve" << endl;
#endif

  Integer i,j;
  Double h;

  RealVector & f = (RealVector &) rhs;
  RealVector & u = (RealVector &) sol;

  y[0] = f.Get(1);
  
  for (i=0; i<size; i++)
    {
      h = 0;
      
      for (j=0; j<=i-1; j++)
	{
	  h += fac[i*size+j]*y[j];
	}
      
      y[i] = f.Get(i+1) - h;
    }
  
  u.Elem(size) = y[size-1]/fac[size*size-1];
  
  for (i=size-2; i>=0; i--)
    {
      h = 0;
      
      for (j=i+1; j<size; j++)
	{
	  h += fac[i*size+j]*u.Get(j+1);
	}
      
      u.Elem(i+1) = (y[i] - h)/fac[i*size+i];
    }
}

void RFullMatrix :: Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Galerkin" << endl;
#endif

  RScalarTransfer & res = (RScalarTransfer &) *ares;
  RScalarTransfer & pro = (RScalarTransfer &) *apro;
  RFullMatrix & mat     = (RFullMatrix &) amat;

  if (mat.IsCalculated())
    {
      return;
    }
  
  Integer i,j,k,p,rs,csize;

  csize = mat.GetSize();

  Double x;
  Double * z = new Double[csize];

  for (i=0; i<size; i++)
    {
      for (j=0; j<csize; j++)
	{
	  z[j] = 0;
	}

      for (j=0; j<size; j++)
	{
	  rs = pro.GetRowSize(j+1);

	  for (k=0; k<rs; k++)
	    {
	      p = pro.Pos(j+1, k+1);
	      z[p-1] += val[i*size+j]*pro.Get(j+1, k+1);
	    }
	}

      rs = res.GetRowSize(i+1);

      for (j=0; j<rs; j++)
	{
	  x = res.Get(i+1,j+1);
	  p = res.Pos(i+1,j+1);

	  for (k=0; k<csize; k++)
	    {
	      mat.Elem(p,k+1) += x*z[k];
	    }
	}
    }

  delete [] z;

  mat.Calculated();
}

void RFullMatrix :: SetAuxMatrix(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::SetAuxMatrix" << endl;
#endif
}

void RFullMatrix :: Copy(BaseMatrix * mat)
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Copy" << endl;
#endif

  RFullMatrix & a = (RFullMatrix &) mat;

  Integer i,j;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
	{
	  val[i*size+j] = a.Get(i+1,j+1);
	}
    }
}

void RFullMatrix :: Print() const
{
#ifdef TRACE
  (*trace) << "entering RFullMatrix::Print" << endl;
#endif

  Integer i,j;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
	{
	  cout << val[i*size+j] << " ";
	}
      cout << endl;
    }
}


CFullMatrix :: CFullMatrix()
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering CFullMatrix::CFullMatrix" << endl;
#endif

}
  
CFullMatrix :: ~CFullMatrix()
{
#ifdef TRACE
  (*trace) << "entering CFullMatrix::~CFullMatrix" << endl;
#endif
}

void CFullMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
}

void CFullMatrix :: Assemble(Double * v,Integer * p,Integer elemsize)
{
}

void CFullMatrix :: Print() const
{
}
}
