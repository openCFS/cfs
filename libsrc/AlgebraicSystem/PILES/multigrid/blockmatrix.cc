#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RBlockMatrix :: RBlockMatrix(Integer asize, Integer anne, Integer adof, Integer adir)
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::RBlockMatrix" << endl;
#endif
  
  size   = asize;
  nne    = anne;
  dof    = adof;
  length = size;
  dir    = adir;
  matfac = 1;

  val    = new Double[nne*dof*dof];
  diaginv= new Double[size*dof*dof];
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

  for (i=0; i<nne*dof*dof; i++)
    {
      val[i] = 0;
    }

  sum = new Double[dof];
  dm  = new DenseMatrix(dof, dof);

  calculated = FALSE;
  setgraph   = FALSE;
  buildindir = FALSE;
  outback    = FALSE;

#ifdef MEMTRACE
  double dmb;
  double imb;

  dmb = (nne*dof*dof+dof+dir+size*dof*dof)*8./1e6;
  imb = (size+1+nne+2*dir)*4./1e6;

  sumdmem += dmb;
  sumimem += imb;

  (*memtrace) << "+++ ALLOCATE MEMORY: double  RBlockMatrix     " << dmb << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer RBlockMatrix     " << imb << " MB" << endl;
#endif
}
  
RBlockMatrix :: ~RBlockMatrix()
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::~RBlockMatrix" << endl;
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

//       (*cla) << "SuperLU" << endl;

//       Destroy_SuperNode_Matrix(&L);
//       Destroy_CompCol_Matrix(&U);
    }

  if (diaginv != NULL)
    {
      delete [] diaginv;
    }

  if (sum != NULL)
    {
      delete [] sum;
    }
}

void RBlockMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Mult" << endl;
#endif

  RealVector & u = (RealVector &) vec1;
  RealVector & v = (RealVector &) vec2;

  DenseMatrix x(dof,dof);

  Integer i,j,k,l,q,p,rs;

  for (i=0; i<size; i++)
    {
      rs  = start[i+1] - start[i];

      for (k=0; k<dof; k++)
	{
	  sum[k] = 0;
	}

      for (j=0; j<rs; j++)
	{
	  p = (pos[start[i]+j]-1)*dof+1;
	  x = *Get(i+1,j+1);

	  for (k=0; k<dof; k++)
	    {
	      for (l=0; l<dof; l++)
		{
		  sum[k] += x.Get(k+1,l+1)*u.Get(p+l);
		}
	    }
	}

      p = i*dof+1;

      for (k=0; k<dof; k++)
	{
	  v.Elem(p+k) = sum[k];
	}
    }
}

void RBlockMatrix :: MultAdd(Double * vec1, BaseVector &vec2) const
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::MultAdd" << endl;
#endif

  RealVector & v = (RealVector &) vec2;

  DenseMatrix x(dof,dof);

  Integer i,j,k,l,q,p,rs;

  for (i=0; i<size; i++)
    {
      rs  = start[i+1] - start[i];

      for (k=0; k<dof; k++)
	{
	  sum[k] = 0;
	}

      for (j=0; j<rs; j++)
	{
	  p = (pos[start[i]+j]-1)*dof+1;
	  x = *Get(i+1,j+1);

	  for (k=0; k<dof; k++)
	    {
	      for (l=0; l<dof; l++)
		{
		  sum[k] += x.Get(k+1,l+1)*vec1[p+l-1];
		}
	    }
	}

      p = i*dof+1;

      for (k=0; k<dof; k++)
	{
	  v.Elem(p+k) += sum[k];
	}
    }
}


void RBlockMatrix :: Assemble(Double * v,Integer * p, Integer elemsize)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Print" << endl;
#endif

  Integer i,j,k,l,m,n,q;
  Integer r,s,t;

  for (i=0; i<elemsize; i++)
    {
      s = start[p[i]-1]*dof*dof;
      n = i*dof*dof*elemsize;

      for (j=0; j<elemsize; j++)
	{
	  t = GetGraphPos(p[i], p[j]);
	  r = s+(t-1)*dof*dof;
	  q = j*dof;

	  for (k=0; k<dof; k++)
	    {
	      m = k*dof*elemsize;

	      for (l=0; l<dof; l++)
		{
		  val[r+k*dof+l] += v[n+q+m+l];
		}
	    }
	}
    }
}

void RBlockMatrix :: Print() const
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Print" << endl;
#endif

  Integer i,j,k,p,q,rs;

  (*cla) << "start printing the matrix ###########################################################" << endl;
  (*cla) << size << endl;
  (*cla) << nne << endl;
  (*cla) << dof << endl;

  for (i=0; i<size; i++)
    {
      p  = start[i];
      q  = p*dof*dof;
      rs = start[i+1]-p;

      for (j=0; j<rs; j++)
	{
	  (*cla) << i+1 << " " << pos[p+j] << " ";

	  for (k=0; k<dof*dof; k++)
	    {
	      (*cla) << val[q+j*dof*dof+k] << " ";
	    }
	  
	  (*cla) << endl;
	} 
    }
  
  (*cla) << "end printing the matrix #############################################################" << endl;
}

void RBlockMatrix :: SetAuxMatrix(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::SetAuxMatrix" << endl;
#endif

  amat.CopyGraph(start, pos);

  RScalarMatrix & mat = (RScalarMatrix &) amat;

  Double v;
  Integer i,j,k,q,rs;

  for (i=0; i<size; i++)
    {
      q  = start[i]*dof*dof;
      rs = start[i+1] - start[i];
      v  = 0;

      for (k=0; k<dof*dof; k++)
	{
	  v += val[q+k]*val[q+k];
	}
      
      v = sqrt(v);
      mat.Elem(i+1, 1) = v;

      for (j=1; j<rs; j++)
	{
	  v  = 0;
	  q += dof*dof;

	  for (k=0; k<dof*dof; k++)
	    {
	      v += val[q+k]*val[q+k];
	    }
	  
	  v = sqrt(v);
	  mat.Elem(i+1, j+1) = -v;
	}
    }
}

void RBlockMatrix :: Copy(BaseMatrix * mat)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Copy" << endl;
#endif
  
  if (mat == this)
    {
      return;
    }

  Integer i,j,k,l,m,rs;

  RBlockMatrix & a = (RBlockMatrix &) *mat;

  m = 0;

  for (i=1; i<=size; i++)
    {
      rs = start[i] - start[i-1];

      for (j=1; j<=rs; j++)
	{
	  for (k=1; k<=dof; k++)
	    {
	      for (l=1; l<=dof; l++)
		{
		  val[m] = a.Get(i,j)->Get(k,l);
		  m++;
		}
	    }
	}
    }

  m = 0;

  for (i=1; i<=size; i++)
    {
      for (j=1; j<=dof; j++)
	{
	  for (k=1; k<=dof; k++)
	    {
	      diaginv[m] = a.GetDiag(i)->Get(j,k);
	      m++;
	    }
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

void RBlockMatrix :: ConstructEffectiveMatrix(BaseMatrix ** amat, Double * matrix_fac)
{
#ifdef TRACE
  (*trace) << "entering RScalarMatrix::ConstructEffectiveMatrix" << endl;
#endif

  Integer i,j,k,l,m,n,rs;

  for (n=0; n<4; n++)
    {
      if (amat[n] != NULL)
	{
	  RBlockMatrix & mat = (RBlockMatrix &) *amat[n];
	  
	  m = 0;
	  
	  for (i=1; i<=size; i++)
	    {
	      rs = start[i] - start[i-1];
	      
	      for (j=1; j<=rs; j++)
		{
		  for (k=1; k<=dof; k++)
		    {
		      for (l=1; l<=dof; l++)
			{
			  val[m] += matrix_fac[n]*mat.Get(i,j)->Get(k,l);
			  m++;
			}
		    }
		}
	    }
	}
    }
}

void RBlockMatrix :: BuildInDirichlet()
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::BuildInDirichlet" << endl;
#endif

  if (buildindir)
    {
      return;
    }

  Integer i,p;
  Double maxdiag,absdiag;
  
  dm      = GetDiag(1);
  maxdiag = dm->MaxNorm();

  for (i=2; i<=size; i++)
    {
      dm      = GetDiag(i);
      absdiag = dm->MaxNorm();

      if (maxdiag < absdiag)
	{
	  maxdiag = absdiag;
	}
    }

  maxdiag *= 1e0;

  for (i=0; i<dir; i++)
    {
      p       = start[numdir[i]-1]*dof*dof+(compdir[i]-1)*(dof+1);
      val[p] += maxdiag;
    }

  buildindir = TRUE;
}

void RBlockMatrix :: BuildInDirichlet(BaseVector & rhs)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::BuildInDirichlet" << endl;
#endif

  if (buildindir || rhs.GetPointer() == NULL)
    {
      return;
    }

  RealVector & a = (RealVector &) rhs;

  Integer i,p,q;
  Double maxdiag,absdiag;
  
  dm      = GetDiag(1);
  maxdiag = dm->MaxNorm();

  for (i=2; i<=size; i++)
    {
      dm      = GetDiag(i);
      absdiag = dm->MaxNorm();

      if (maxdiag < absdiag)
	{
	  maxdiag = absdiag;
	}
    }

  maxdiag *= 1e0;

  for (i=0; i<dir; i++)
    {
      p          = start[numdir[i]-1]*dof*dof+(compdir[i]-1)*(dof+1);
      q          = (numdir[i]-1)*dof+compdir[i];
      val[p]    += maxdiag;
      a.Elem(q) += maxdiag*valdir[i];
    }

  buildindir = TRUE;
}

void RBlockMatrix :: Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Galerkin" << endl;
#endif
  RBlockTransfer & res = (RBlockTransfer &) *ares;
  RBlockTransfer & pro = (RBlockTransfer &) *apro;
  RBlockMatrix & mat   = (RBlockMatrix &) amat;

  if (mat.IsCalculated())
    {
      return;
    }

  Integer i,j,k,l,m,n,p,q,rs,maxrs,p1,p2;

  DenseMatrix * x = new DenseMatrix(dof,dof);
  DenseMatrix * y = new DenseMatrix(dof,dof,TRUE);
  DenseMatrix * z = new DenseMatrix(dof,dof);

  maxrs = mat.MaxRowSize();

  LocalSet v(maxrs,dof);

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
	  x = Get(i+1,j+1); // val[q];

	  for (k=0; k<n; k++)
	    {
	      p = pro.Pos(pos[q], k+1);
	      z = pro.Get(pos[q], k+1);
	      y->Mult(x,z);

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
		      z = v.GetBlock(k);
		      y->MultT(x,z); 

		      for (p1=1; p1<=dof; p1++)
			{
			  for (p2=1; p2<=dof; p2++)
			    {
			      mat.Elem(p,l,p1,p2) += y->Get(p1,p2);
			    }
			}
		    }
		}
	    }
	}
    }

  delete y;

  mat.Calculated();
}

void RBlockMatrix :: SetDiagInv()
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::SetDiagInv" << endl;
#endif

  DenseMatrix * x = new DenseMatrix(dof,dof);
  DenseMatrix * y = new DenseMatrix(dof,dof);
  DenseMatrix * z = new DenseMatrix(dof,dof,TRUE);
  Integer i,j,k,p;

  for (i=0; i<size; i++)
    {
      p = i*dof*dof;
      x = GetDiag(i+1);
      *y= &diaginv[p];
      y->Inverse(x);
    }

  delete z;
}

void RBlockMatrix :: Factor()
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Factor" << endl;
#endif

  int i,j,k,l,m,numrhs,dofsize,dofnne,ind,rs;
  int osb,osr,rsold;
  
  int permc_spec,info;

  numrhs = 1;
  dofsize= dof*size;
  dofnne = dof*dof*nne;

  f      = new double[dofsize];
  perm_c = new int[dofsize];
  perm_r = new int[dofsize];

  Double * dofval    = new Double[dofnne];
  Integer * dofstart = new Integer[dofsize+1];
  Integer * dofpos   = new Integer[dofnne];

  dofstart[0] = 0;
  k = 1;

  for (i=0; i<size; i++)
    {
      rs = start[i+1] - start[i];

      for (j=0; j<dof; j++)
	{
	  dofstart[k] = dofstart[k-1]+rs*dof;
	  k++;
	}
    }
  
  /// + construct the pos vector and the val vector

  int * os = new int[dof];
  DenseMatrix x(dof,dof);

  rsold = 0;
  osb   = 0;

  for (i=0; i<size; i++)
    {
      rs  = start[i+1] - start[i];
      osb += rsold*dof*dof;
      osr = dof*rs;

      for (j=0; j<dof; j++)
	{
	  os[j] = 0;
	}

      for (j=0; j<rs; j++)
	{
	  x = *Get(i+1,j+1);
	  m  = GetMatrixPos(i+1,j+1);

	  for (k=0; k<dof; k++)
	    {
	      for (l=0; l<dof; l++)
		{
		  ind         = osb+os[k]+k*osr+l;
		  dofpos[ind] = (m-1)*dof+l;
		  dofval[ind] = x.Get(k+1,l+1);
		}
	      
	      os[k] += dof;
	    }
	}

      rsold = rs;
    }

  for (i=0; i<dofsize; i++)
    {
      f[i] = 1;
    }

  dCreate_CompCol_Matrix(&A, dofsize, dofsize, dofnne, dofval, dofpos, dofstart, NR, _D, GE);
  dCreate_Dense_Matrix(&B, dofsize, numrhs, f, dofsize, DN, _D, GE);

  permc_spec = 3;

  get_perm_c(permc_spec, &A, perm_c);

  dgssv(&A, perm_c, perm_r, &L, &U, &B, &info);

  delete [] dofval;
  delete [] dofstart;
  delete [] dofpos;
  delete [] os;
}

void RBlockMatrix :: Solve(BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering RBlockMatrix::Solve" << endl;
#endif

  RealVector & r = (RealVector &) rhs;
  RealVector & u = (RealVector &) sol;

  int i,numrhs,info;
  int dofsize;
  double * h, *g;
  char t[1];


  *t = 'T';
  
  dofsize= size*dof;
  numrhs = 1;
  h      = r.GetPointer();
  g      = u.GetPointer();

  for (i=0; i<dofsize; i++)
    {
      f[i] = h[i];
    }

  dCreate_Dense_Matrix(&B, dofsize, numrhs, f, dofsize, DN, _D, GE);
  
  dgstrs (t, &L, &U, perm_r, perm_c, &B, &info);
  
  DNformat * vec1 = (DNformat *) B.Store;
  double * vec2 = (double *) vec1->nzval;
  
  for (i=0; i<dofsize; i++)
    {
      u.Elem(i+1) = vec2[i];//B.Store;
    }
}

/////////////////////////////////////////////////////////////////////////////////////

CBlockMatrix :: CBlockMatrix()
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering CBlockMatrix::CBlockMatrix" << endl;
#endif

}
  
CBlockMatrix :: ~CBlockMatrix()
{
#ifdef TRACE
  (*trace) << "entering CBlockMatrix::~CBlockMatrix" << endl;
#endif
}

void CBlockMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
}

void CBlockMatrix :: Assemble(Double * v,Integer * p, Integer elemsize)
{
}

void CBlockMatrix :: Print() const
{
}
}
