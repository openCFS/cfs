#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

extern "C"
{
  //extern void extern_get_sl_const_diag (long *, double **);
  //extern void extern_constr_hs_diag (long *, double **);
  //extern void extern_matrix_mult_sl(long , double *, double *);
  //extern void extern_matrix_mult_hs(long , double *, double *);
}

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

OutbackMatrix :: OutbackMatrix(Integer asize, Integer afsize, Integer adir, Integer alevel)
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::OutbackMatrix" << endl;
#endif

  Integer i;

  size   = asize;
  nne    = 0;
  fsize  = afsize;
  hsize  = fsize;
  dir    = adir;
  level  = alevel;
  matfac = 1;

  accu = NULL;
  rsw  = new Integer[fsize];
  r    = new RealVector(fsize, 1, 1);
  s    = new RealVector(fsize, 1, 1);

  (*cla) << "making the outbackmatrix" << endl;

  if (level == 0)
    {
      //extern_constr_hs_diag (&hsize, &diag);

//       diag = new Double[size];

//       for (i=0; i<size; i++)
// 	{
// 	  diag[i] = 1;
// 	}
    }
  else
    {
      diag = new Double[size];

      for (i=0; i<size; i++)
	{
	  diag[i] = 0;
	}
    }

  val     = NULL;
  pos     = NULL;
  start   = NULL;
  cm      = NULL;
  f       = NULL;
  perm_c  = NULL;
  perm_r  = NULL;
  valdir  = NULL;
  numdir  = NULL;
  compdir = NULL;
    
  outback    = TRUE;
  calculated = FALSE;
  setgraph   = TRUE;
  buildindir = TRUE;
}
  
OutbackMatrix :: ~OutbackMatrix()
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::~OutbackMatrix" << endl;
#endif

  delete [] diag;
  delete [] rsw;
  delete accu;
  delete r;
  delete s;
}

void OutbackMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Mult" << endl;
#endif

  Double * p1, * p2, * pr, * ps;

  p1 = vec1.GetPointer();
  p2 = vec2.GetPointer();
  pr = r->GetPointer();
  ps = s->GetPointer();

  if (level == 0)
    {
      //(*cla) << "mult on level " << level << endl;
      //extern_matrix_mult_sl(hsize, p1, p2);
      //extern_matrix_mult_hs(hsize, p1, p2);
    }
  else
    {
      //(*cla) << "mult on level " << level << endl;
      s->Init();
      accu->MultHh(vec1, *r);
      //extern_matrix_mult_sl(hsize, pr, ps);
      //extern_matrix_mult_hs(hsize, pr, ps);
      accu->MulthH(*s, vec2);
    }
}

void OutbackMatrix :: BuildInDirichlet()
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::BuildInDirichlet" << endl;
#endif
  
}

void OutbackMatrix :: BuildInDirichlet(BaseVector & rhs)
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::BuildInDirichlet" << endl;
#endif


}

void OutbackMatrix :: Factor()
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Factor" << endl;
#endif
}

void OutbackMatrix :: Solve(BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Solve" << endl;
#endif

  //(*cla) << "### coarse grid solver" << endl;

  Integer i,j,iter;
  Integer numsolve = 100;
  Double damp = 15;
  Double eps = 1e-2;
  Double epsmach = 1e-12;

  RealVector & f = (RealVector &) rhs;
  RealVector & u = (RealVector &) sol;
  RealVector * r = new RealVector(size, 1, 1);
  RealVector * d = new RealVector(size, 1, 1);
  RealVector * s = new RealVector(size, 1, 1);

  //(*cla) << "### preconditioned CG Solver" << endl;

  DenseVector alpha(1);
  DenseVector beta(1);
  DenseVector gamma(1);
  DenseVector wra(1);
  DenseVector wrn(1);
  DenseVector err(1);
  DenseVector rel(1);
  DenseVector scal(1);
  DenseVector pone(1);
  DenseVector mone(1);
  DenseVector zero(1);

  r->Init();
  s->Init();
  d->Init();

  Boolean loop;

  pone = +1.0;
  mone = -1.0;
  zero = 0;

  iter = 1;
  loop = TRUE;

  Mult(sol, *d, 1);

  r->Add(pone, rhs, mone, *d);
  s->Add(pone, *r, zero, *r);
  s->Inner(*r, wra);

  scal = wra;
  rel  = 1;
  err  = wra;
  err *= (eps*eps);

  if ((wra <= err) || (wra <= (epsmach*epsmach)))
    {
      loop = FALSE;
    }

  //(*cla) << 1 << " " << rel.Get(1) << " " << wra.Get(1) << endl;
  //(*conv) << rel << endl;

  while (loop && (iter <= size))
    {
      iter++;
      Mult(*s, *d, 1);
      d->Inner(*s, gamma);
      alpha = wra/gamma;

      sol.Add(alpha, *s);

      alpha *= -1;

      r->Add(alpha, *d);
      d->Add(pone, *r, zero, *r);

      d->Inner(*r, wrn);
      
      if ((wrn <= err) || (wrn <= (epsmach*epsmach)))
	{
	  loop =FALSE;
	}

      beta = wrn/wra;
      wra  = wrn;
      s->Add(beta, *s,pone, *d);

      if (wra <= -1e-32)
	{
	  (*cla) << "... EXIT in PCG-SOLVER " << wra.Get(1) << endl;
	  exit(1);
	}

      rel = wra/scal;
      rel.Sqrt();
	      
      //(*cla) << iter << " " << rel.Get(1) << " " << wra.Get(1) << endl;
    }

  delete r;
  delete d;
  delete s;

//   for (i=0; i<numsolve; i++)
//     {
//       r.Init();

//       Mult(u, r, 1);

//       for (j=1; j<=size; j++)
// 	{
// 	  u.Elem(j) += (damp*(f.Get(j) - r.Get(j)));///diag[j-1]);
// 	}
//     }

}

void OutbackMatrix :: Galerkin(BaseTransfer * ares, BaseTransfer * apro, BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Galerkin" << endl;
#endif

  Integer i,j,k,p,q,csize,rs,ss,maxrs;
  Double val,x,y;

  RScalarTransfer & res  = (RScalarTransfer &) *ares;
  RScalarTransfer & pro  = (RScalarTransfer &) *apro;
  RScalarTransfer & tra  = (RScalarTransfer &) *accu;
  OutbackMatrix & mat    = (OutbackMatrix &) amat;

  if (mat.IsCalculated())
    {
      return;
    }

  if (level == 0)
    {
      mat.SetTransfer(apro);

      /// calculate the diagonal matrix
      
      for (i=0; i<size; i++)
	{
	  rs = pro.GetRowSize(i+1);
	  
	  for (j=0; j<rs; j++)
	    {
	      val = diag[i]*pro.Get(i+1, j+1)*pro.Get(i+1, j+1);
	      p = pro.Pos(i+1, j+1);
	      
	      mat.Diag(p) += val;
	    }
	}
    
      return;
    }

  csize = mat.GetSize();

  maxrs = 1000;

  LocalSet v(maxrs);

  /// calc the rsw

  for (i=0; i<size; i++)
    {
      v.Init();

      rs = tra.GetRowSize(i+1);

      for (j=0; j<rs; j++)
	{
	  p  = tra.Pos(i+1, j+1);
	  ss = pro.GetRowSize(p);

	  for (k=0; k<ss; k++)
	    {
	      q = pro.Pos(p,k+1);

	      v.Insert(q);
	    }
	}

      rsw[i] = v.GetSize();
    }

  (*cla) << "before create transfer " << endl;

  mat.CreateTransfer(size, rsw);

  (*cla) << "after create transfer " << size << endl;

  /// perform the matrix multiplication

  for (i=0; i<size; i++)
    {
      v.Init();

      rs = tra.GetRowSize(i+1);

      for (j=0; j<rs; j++)
	{
	  p  = tra.Pos(i+1, j+1);
	  ss = pro.GetRowSize(p);

	  for (k=0; k<ss; k++)
	    {
	      q = pro.Pos(p,k+1);

	      v.Insert(q);
	    }
	}

      v.Sort();
      
      (*cla) << "i " << i << endl;

      rs = tra.GetRowSize(i+1);

      for (j=0; j<rs; j++)
	{
	  p  = tra.Pos(i+1, j+1);
	  x  = tra.Get(i+1, j+1);
	  ss = pro.GetRowSize(p);

	  for (k=0; k<ss; k++)
	    {
	      q = pro.Pos(p, k+1);
	      y = x*pro.Get(p, k+1);

	      v.Insert(q,y);
	    }
	}
      
      tra.SetRow(i+1, v);
    }

  /// calculate the diagonal matrix

//   for (i=0; i<size; i++)
//     {
//       rs = pro.GetRowSize(i+1);
      
//       for (j=0; j<rs; j++)
// 	{
// 	  val = diag[i]*pro.Get(i+1, j+1)*pro.Get(i+1, j+1);
// 	  p = pro.Pos(i+1, j+1);

// 	  mat.Diag(p) += val;
// 	}
//     }

//   for (i=0; i<csize; i++)
//     {
//       (*cla) << mat.GetDiag(i+1) << " ";
//     }
//   (*cla) << endl;

  mat.Calculated();
}

void OutbackMatrix :: SetAuxMatrix(BaseMatrix & amat)
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::SetAuxMatrix" << endl;
#endif
}

void OutbackMatrix :: Copy(BaseMatrix * mat)
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Copy" << endl;
#endif

}

void OutbackMatrix :: Print() const
{
#ifdef TRACE
  (*trace) << "entering OutbackMatrix::Print" << endl;
#endif

}
}
