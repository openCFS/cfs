#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RScalarMatrix :: RScalarMatrix(Integer asize, Integer anne, Integer adir)
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::RScalarMatrix" << endl;
#endif

  size   = asize;
  nne    = anne;
  dof    = 1;
  length = size;
  dir    = adir;
  matfac = 1;

  val    = new Double[nne];
  pos    = new Integer[nne];
  start  = new Integer[size+1];

  if (dir != 0)
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

  cm     = NULL;
  f      = NULL;
  perm_c = NULL;
  perm_r = NULL;

  start[0] = 0;

  Integer i;

  for (i=0; i<nne; i++)
    {
      val[i] = 0;
    }

#ifdef MEMTRACE
  double dmb;
  double imb;

  dmb = (nne+dir)*8./1e6;
  imb = (nne+size+1+2*dir)*4./1e6;

  sumdmem += dmb;
  sumimem += imb;

  (*memtrace) << "+++ ALLOCATE MEMORY: double  ScalarMatrix     " << dmb << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer ScalarMatrix     " << imb << " MB" << endl;
#endif


  calculated = FALSE;
  setgraph   = FALSE;
  buildindir = FALSE;
  outback    = FALSE;
}
  
RScalarMatrix :: ~RScalarMatrix()
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::~RScalarMatrix" << endl;
#endif

  if (val != NULL)
    {
      delete [] val;
      delete [] pos;
      delete [] start;
    }

  if (valdir != NULL)
    {
      delete [] valdir;
      delete [] numdir;
      delete [] compdir;
    }

  if (cm != NULL)
    {
      delete [] cm;
    }

  if (perm_c != NULL)
    {
      delete [] perm_c;
      delete [] perm_r;
      delete [] f;
    }
}

void RScalarMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::Mult" << endl;
#endif

  Integer i,j,k,rs;
  Double sum;

  RealVector & u = (RealVector &) vec1;
  RealVector & v = (RealVector &) vec2;

  for (i=0; i<size; i++)
    { 
      sum = 0;
      
      rs = start[i+1]-start[i];
      k  = start[i];

      for (j=0; j<rs; j++)
	{
	  sum += val[k]*u.Get(pos[k]);
	  k++;
	}
	
	v.Elem(i+1) = factor*sum;
    }
}


double RScalarMatrix :: CalcEnergyNorm(Double * vec1) const
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::CalcEnergyNorm" << endl;
#endif

  Integer i,j,k,rs;
  Double sum, norm;

  norm = 0.0;
  for (i=0; i<size; i++)
    { 
      sum = 0;
      
      rs = start[i+1]-start[i];
      k  = start[i];

      for (j=0; j<rs; j++)
	{
	  sum += val[k]*vec1[pos[k]-1];
	  k++;
	}
	
      norm = norm + sum*vec1[i];
    }

  return norm;
}


void RScalarMatrix :: MultAdd(Double * vec1, BaseVector &vec2) const
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::MultAdd" << endl;
#endif

  Integer i,j,k,rs;
  Double sum;

  RealVector & v = (RealVector &) vec2;

  for (i=0; i<size; i++)
    { 
      sum = 0;
      
      rs = start[i+1]-start[i];
      k  = start[i];

      for (j=0; j<rs; j++)
	{
	  sum += val[k]*vec1[pos[k]-1];
	  k++;
	}
	
	v.Elem(i+1) += sum;
    }
}

void RScalarMatrix :: Assemble(Double * v,Integer * p,Integer elemsize)
{
  Integer i,j,k,p1,p2,p3;

  for (j=0; j<elemsize; j++)
    {
      if (p[j] != 0)
	{
	  if (p[j] > 0)
	    {
	      p1 = p[j];
	    }
	  else
	    {
	      p1 = cm[abs(p[j])];
	    }
	  
	  for (k=0; k<elemsize; k++)
	    {
	      if (p[k] != 0)
		{
		  if (p[k] > 0)
		    {
		      p2 = GetGraphPos(p1,p[k]);
		    }
		  else
		    {
		      p3 = cm[abs(p[k])];
		      p2 = GetGraphPos(p1,p3);
		    }

		  val[start[p1-1]+p2-1] += v[j*elemsize+k];
		}
	    }
	}
    }
}

void RScalarMatrix :: Print() const
{
  Integer i,j,rs;

  (*cla) << "start printing the matrix ###########################################################" << endl;
  (*cla) << size << endl;
  (*cla) << nne << endl;
  (*cla) << 1 << endl;

  for (i=0; i<size; i++)
    {
      rs = start[i+1]-start[i];

      for (j=0; j<rs; j++)
	{
	  (*cla) << i+1 << " " << pos[start[i]+j] << " " << val[start[i]+j] << endl;
	} 
    }
  
  (*cla) << "end printing the matrix #############################################################" << endl;
}

void RScalarMatrix :: BuildInDirichlet()
{
  if (buildindir)
    {
      return;
    }

  Integer i;
  Double maxdiag;
  
  maxdiag = val[start[0]];

  for (i=1; i<size; i++)
    {
      if (maxdiag < val[start[i]])
	{
	  maxdiag = val[start[i]];
	}
    }

  maxdiagentry = maxdiag;

  maxdiag *= 1e12;

  for (i=0; i<dir; i++)
    {
      val[start[numdir[i]-1]] += maxdiag;
    }

  buildindir = TRUE;
}

void RScalarMatrix :: BuildInDirichlet(BaseVector & rhs)
{
  if (buildindir)
    {
      return;
    }

  RealVector & a = (RealVector &) rhs;

  Integer i;
  Double maxdiag;
  
  maxdiag = val[start[0]];

  for (i=1; i<size; i++)
    {
      if (maxdiag < val[start[i]])
	{
	  maxdiag = val[start[i]];
	}
    }

  maxdiagentry = maxdiag;
  maxdiag *= 1e12;

  for (i=0; i<dir; i++)
    {
      val[start[numdir[i]-1]] += maxdiag;
      a.Elem(numdir[i])       += maxdiag*valdir[i];
    }

  buildindir = TRUE;
}

void RScalarMatrix :: UpdateDirichletRHS(BaseVector &rhs)
{
  RealVector & a = (RealVector &) rhs;

  Integer i;
  Double maxdiag;
  
  maxdiag = maxdiagentry*1e12;

  for (i=0; i<dir; i++)
    {
      a.Elem(numdir[i]) += maxdiag*valdir[i];
    }
}

void RScalarMatrix :: Factor()
{
  int i,j,numrhs;
  int permc_spec,info;

  numrhs = 1;

  f      = new double[size];
  perm_c = new int[size];
  perm_r = new int[size];

  for (i=0; i<size; i++)
    {
      f[i] = 1;
    }
  
  for (i=0; i<nne; i++)
    {
      pos[i]--;
    }

  dCreate_CompCol_Matrix(&A, size, size, nne, val, pos, start, NR, _D, GE);
  dCreate_Dense_Matrix(&B, size, numrhs, f, size, DN, _D, GE);

  permc_spec = 3;

  get_perm_c(permc_spec, &A, perm_c);
  dgssv(&A, perm_c, perm_r, &L, &U, &B, &info);

  for (i=0; i<nne; i++)
    {
      pos[i]++;
    }
}

void RScalarMatrix :: Solve(BaseVector & rhs, BaseVector & sol)
{
  RealVector & r = (RealVector &) rhs;
  RealVector & u = (RealVector &) sol;
  
  int i,numrhs,info;
  double * h, *g;
  char t[1];
  
  *t = 'T';
  
  numrhs = 1;
  h      = r.GetPointer();
  g      = u.GetPointer();

  for (i=0; i<size; i++)
    {
      f[i] = h[i];
    }
  
  dCreate_Dense_Matrix(&B, size, numrhs, f, size, DN, _D, GE);
  
  dgstrs (t, &L, &U, perm_r, perm_c, &B, &info);
  
  DNformat * vec1 = (DNformat *) B.Store;
  double * vec2 = (double *) vec1->nzval;
  
  for (i=0; i<size; i++)
    {
      u.Elem(i+1) = vec2[i];//B.Store;
    }
}

Double RScalarMatrix :: MinRowVal(Integer i) const
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::MinRowVal" << endl;
#endif

  Integer j,rs;
  Double min = 1e32;

  rs = start[i]-start[i-1];

  for (j=0; j<rs; j++)
    {
      if (val[start[i-1]+j] < min)
	{
	  min = val[start[i-1]+j];
	}
    }

  return min;
}

void RScalarMatrix :: Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::Galerkin" << endl;
#endif

  RScalarTransfer & res = (RScalarTransfer &) *ares;
  RScalarTransfer & pro = (RScalarTransfer &) *apro;
  RScalarMatrix & mat   = (RScalarMatrix &) amat;

  if (mat.IsCalculated())
    {
      return;
    }

  Integer i,j,k,l,m,n,p,q,rs,maxrs;
  Double x,y;

  maxrs = mat.MaxRowSize();

  LocalSet v(maxrs);

  for (i=0; i<size; i++)
    {
      q  = start[i];
      rs = start[i+1] - q;
      
      v.Init();

      for (j=0; j<rs; j++)
	{
	  n = pro.GetRowSize(pos[q]);

	  for (k=0; k<n; k++)
	    {
	      p = pro.Pos(pos[q], k+1);

	      v.Insert(p);
	    }

	  q++;
	}
      
      v.Sort();

      q  = start[i];
      rs = start[i+1] - q;

      for (j=0; j<rs; j++)
	{
	  n = pro.GetRowSize(pos[q]);
	  x = val[q];

	  for (k=0; k<n; k++)
	    {
	      p = pro.Pos(pos[q], k+1);
	      y = x*pro.Get(pos[q], k+1);

	      v.Insert(p, y);
	    }

	  q++;
	}

      rs = res.GetRowSize(i+1);
      m  = v.GetSize();

      for (j=0; j<rs; j++)
	{
	  x = res.Get(i+1,j+1);
	  p = res.Pos(i+1,j+1);
	  n = mat.GetRowSize(p);

	  for (k=1; k<=m; k++)
	    {
	      for (l=1; l<=n; l++)
		{
		  if (v.Pos(k) == mat.GetMatrixPos(p,l))
		    {
		      mat.Elem(p,l) += (x*v.Get(k));
		    }
		}
	    }
	}
    }

  mat.Calculated();
}

void RScalarMatrix :: Copy(BaseMatrix * mat)
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::Copy" << endl;
#endif
  
  if (mat == this)
    {
      return;
    }

  Integer i,j,k,rs;

  RScalarMatrix & a = (RScalarMatrix &) *mat;

  k = 0;

  for (i=0; i<size; i++)
    {
      rs = start[i+1] - start[i];

      for (j=0; j<rs; j++)
	{
	  val[k] = a.Get(i+1,j+1);
	  k++;
	}
    }

  for (i=0; i<dir; i++)
    {
      valdir[i] = a.GetValDirichlet(i+1);
      numdir[i] = a.GetIndDirichlet(i+1);
      compdir[i]= a.GetCmpDirichlet(i+1);
    }

  buildindir = FALSE;
}

void RScalarMatrix :: ConstructEffectiveMatrix(BaseMatrix ** amat, Double * matrix_fac)
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::ConstructEffectiveMatrix" << endl;
#endif

  Integer i,j,k,l,rs;

  for (l=0; l<4; l++)
    {
      if (amat[l] != NULL)
	{
	  RScalarMatrix & mat = (RScalarMatrix &) *amat[l];
	  
	  k = 0;
	  
	  for (i=0; i<size; i++)
	    {
	      rs = start[i+1] - start[i];
	      
	      for (j=0; j<rs; j++)
		{
		  val[k] += matrix_fac[l]*mat.Get(i+1,j+1);
		  k++;
		}
	    }
	}
    }
}

////////////////////////////////// complex ///////////////////////////////////////



CScalarMatrix :: CScalarMatrix(Integer asize, Integer anne, Integer adir)
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering CScalarMatrix::CScalarMatrix" << endl;
#endif
  
  size   = asize;
  nne    = anne;
  dof    = 1;
  length = size;
  dir    = adir;

  val   = new Double[2*nne];
  pos   = new Integer[nne];
  start = new Integer[size+1];

  valdir  = new Double[dir];
  numdir  = new Integer[dir];
  compdir = new Integer[dir];

  start[0] = 0;

  Integer i;

  for (i=0; i<2*nne; i++)
    {
      val[i] = 0;
    }

  calculated = FALSE;
  setgraph   = FALSE;
  buildindir = FALSE;
  outback    = FALSE;
}
  
CScalarMatrix :: ~CScalarMatrix()
{
#ifdef TRACE
  (*trace) << "entering CScalarMatrix::~CScalarMatrix" << endl;
#endif
 
  if (val != NULL && !outback)
    {
      delete [] val;
      delete [] pos;
      delete [] start;
      delete [] valdir;
      delete [] numdir;
      delete [] compdir;
    }
}

void CScalarMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
  Integer i,j,k,l,m,n,rs;
  Double re,im;

  ComplexVector & u = (ComplexVector &) vec1;
  ComplexVector & v = (ComplexVector &) vec2;

  l = 1;
  k = 0;

  for (i=0; i<size; i++)
    { 
      re = 0;
      im = 0;
      
      rs = start[i+1]-start[i];
      m  = start[i];

      for (j=0; j<rs; j++)
	{
	  n = 2*pos[m+j]-1;

	  re += (val[k]*u.Get(n) - val[k+1]*u.Get(n+1));
	  im += (val[k]*u.Get(n+1) + val[k+1]*u.Get(n));

	  k += 2;
	}
	
      v.Elem(l)   = factor*re;
      v.Elem(l+1) = factor*im;

      l += 2;
    }
}

void CScalarMatrix :: Assemble(Double * v,Integer * p,Integer elemsize)
{
  Integer i,j,k,p1,p2,p3,p4,p5;

  for (j=0; j<elemsize; j++)
    {
      if (p[j] != 0)
	{
	  if (p[j] > 0)
	    {
	      p1 = p[j];
	    }
	  else
	    {
	      p1 = cm[abs(p[j])];
	    }
	  
	  for (k=0; k<elemsize; k++)
	    {
	      if (p[k] != 0)
		{
		  if (p[k] > 0)
		    {
		      p2 = GetGraphPos(p1,p[k]);
		    }
		  else
		    {
		      p3 = cm[abs(p[k])];
		      p2 = GetGraphPos(p1,p3);
		    }

		  p4 = 2*start[p1-1]+2*(p2-1);
		  p5 = j*elemsize+k;

		  val[p4]   += v[p5];
		  val[p4+1] += v[p5+elemsize*elemsize];
		}
	    }
	}
    }
}

void CScalarMatrix :: Print() const
{
  Integer i,j,k,rs;

  (*test) << "printing the complex stiffness matrix" << endl;
  (*test) << "size " << size << endl;

  k = 0;

  for (i=0; i<size; i++)
    {
      rs = start[i+1]-start[i];

      for (j=0; j<rs; j++)
	{
	  (*test) << val[k] << " " << val[k+1] << ",";
	  k += 2;
	}
      
      (*test) << endl;

      for (j=0; j<rs; j++)
	{
	  (*test) << pos[start[i]+j] << " ";
	}
      
      (*test) << endl;
    }
}

void CScalarMatrix :: BuildInDirichlet()
{

}

void CScalarMatrix :: BuildInDirichlet(BaseVector & rhs)
{
  ComplexVector & a = (ComplexVector &) rhs;

  Integer i, j;

  for (i=0; i<dir; i++)
    {
      j = 2*start[numdir[i]-1];

      val[j]   += 1e+6;
      val[j+1] += 1e+6;
    }
}

void CScalarMatrix :: Factor()
{

}

void CScalarMatrix :: Solve(BaseVector & rhs, BaseVector & sol)
{

}

void CScalarMatrix :: Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat)
{

}

void CScalarMatrix :: SetAuxMatrix(BaseMatrix & amat)
{

}

void CScalarMatrix :: Copy(BaseMatrix * mat)
{

}
}
