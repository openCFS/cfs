#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

CGSolver :: CGSolver(BaseParameter & param, Integer size, Integer anumrhs, Integer dof)
  : BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering CGSolver::CGSolver" << endl;
#endif

  numrhs  = anumrhs;
  eps     = param.GetEps();
  epsmach = param.GetEpsMach();
  maxiter = param.GetMaxNumIter();

  r = new RealVector(size, numrhs, dof);
  s = new RealVector(size, numrhs, dof);
  d = new RealVector(size, numrhs, dof);

  aalpha = new DenseVector(numrhs);
  abeta  = new DenseVector(numrhs);
  agamma = new DenseVector(numrhs);
  awra   = new DenseVector(numrhs);
  awrn   = new DenseVector(numrhs);
  aerr   = new DenseVector(numrhs);
  arel   = new DenseVector(numrhs);
  ascal  = new DenseVector(numrhs);
  apone  = new DenseVector(numrhs);
  amone  = new DenseVector(numrhs);
}
  
CGSolver :: ~CGSolver()
{
#ifdef TRACE
  (*trace) << "entering CGSolver::~CGSolver" << endl;
#endif

  delete r;
  delete s;
  delete d;

  delete aalpha;
  delete abeta;
  delete agamma;
  delete awra;
  delete awrn;
  delete aerr;
  delete arel;
  delete ascal;
  delete apone;
  delete amone;
}

void CGSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		      BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering CGSolver::Calc" << endl;
#endif

  (*cla) << "### preconditioned CG Solver" << endl;

  DenseVector & alpha = (DenseVector &) *aalpha;
  DenseVector & beta  = (DenseVector &) *abeta;
  DenseVector & gamma = (DenseVector &) *agamma;
  DenseVector & wra   = (DenseVector &) *awra;
  DenseVector & wrn   = (DenseVector &) *awrn;
  DenseVector & err   = (DenseVector &) *aerr;
  DenseVector & rel   = (DenseVector &) *arel;
  DenseVector & scal  = (DenseVector &) *ascal;
  DenseVector & pone  = (DenseVector &) *apone;
  DenseVector & mone  = (DenseVector &) *amone;

  r->Init();
  s->Init();
  d->Init();

  Boolean loop;

  pone = +1.0;
  mone = -1.0;

  iter = 1;
  loop = TRUE;

  // generate an initial guess

  premat.Step(rhs, sol);

  sysmat.Mult(sol, *d);

  r->Add(pone, rhs, mone, *d);

  premat.Step(*r, *s);

  s->Inner(*r, wra);

  scal = wra;
  rel  = 1;
  err  = wra;
  err *= (eps*eps);

  if ((wra <= err) || (wra <= (epsmach*epsmach)))
    {
      loop = FALSE;
    }

  (*cla) << 1 << " " << rel.Get(1) << " " << wra.Get(1) << endl;
  //(*conv) << rel << endl;

  while (loop && (iter <= maxiter))
    {
      iter++;
      sysmat.Mult(*s, *d);
      d->Inner(*s, gamma);
      alpha.Div(wra, gamma);
      sol.Add(alpha, *s);

      alpha *= -1;

      r->Add(alpha, *d);
      d->Init();
 
      premat.Step(*r, *d);

      d->Inner(*r, wrn);
      
      if ((wrn <= err) || (wrn <= (epsmach*epsmach)))
	{
	  loop =FALSE;
	}

      beta.Div(wrn, wra);
      wra  = wrn;
      s->Add(beta, *s,pone, *d);

      if (wra <= -1e-32)
	{
	  (*cla) << "... EXIT in PCG-SOLVER " << wra.Get(1) << endl;
	  exit(1);
	}

      rel.Div(wra, scal);
      rel.Sqrt();
	      
      (*cla) << iter << " " << rel.Get(1) << " " << wra.Get(1) << endl;
      //(*conv) << rel << endl;
    }

  (*cla) << "... number of iterations = " << iter << endl;
}

void CGSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering CGSolver::CalcRate" << endl;
#endif
}

////////////////////////////////////////////////////////////////////////////////////

LanczosSolver :: LanczosSolver(BaseParameter & param, Integer size, Integer anumrhs, Integer dof)
  : BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering LanczosSolver::LanczosSolver" << endl;
#endif

  numrhs  = anumrhs;
  eps     = param.GetEps();
  epsmach = param.GetEpsMach();
  maxiter = param.GetMaxNumIter();

  s = new RealVector(size, numrhs, dof);
  d = new RealVector(size, numrhs, dof);

  agamma = new DenseVector(numrhs);
  aerr   = new DenseVector(numrhs);
  ascal  = new DenseVector(numrhs);
  apone  = new DenseVector(numrhs);
  azero  = new DenseVector(numrhs);
}
  
LanczosSolver :: ~LanczosSolver()
{
#ifdef TRACE
  (*trace) << "entering LanczosSolver::~LanczosSolver" << endl;
#endif

  delete s;
  delete d;

  delete agamma;
  delete aerr;
  delete ascal;
  delete apone;
  delete azero;

  delete [] diag;
  delete [] subd;
}

void LanczosSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		      BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering LanczosSolver::Calc" << endl;
#endif

  (*cla) << "### preconditioned Lanczos Solver" << endl;

  DenseVector & gamma = (DenseVector &) *agamma;
  DenseVector & err   = (DenseVector &) *aerr;
  DenseVector & scal  = (DenseVector &) *ascal;
  DenseVector & pone  = (DenseVector &) *apone;
  DenseVector & zero  = (DenseVector &) *azero;

  pone = 1;
  zero = 0;

  step = 200;
  prec = 1e-10;

  diag = new Double[step];
  subd = new Double[step];

  sol.Init();
  rhs.Rand();

  premat.Step(rhs, sol);

  sol.Inner(rhs, gamma);

  if (gamma.Get(1) < -1e-20)
    {
      (*cla) << "LanczosSolve:Calc" << endl;
      (*cla) << "preconditioner indefinite!" << endl;
      exit(1);
    }

  gamma.Sqrt();
  scal.Div(pone, gamma);

  sol.Scal(scal);
  rhs.Scal(scal);

  sysmat.Mult(sol, *d);

  iter    = 0;
  subd[1] = 0;
  err     = 1;

  while (!(err <= prec) && iter < step)
    {
      iter++;

      (*cla) << ".";

      sol.Inner(*d, scal);
      
      diag[iter] = scal.Get(1);

      scal *= -1;

      s->Add(pone, *d, scal, rhs);
      
      sol.Init();
      premat.Step(*s, sol);
      
      s->Inner(sol, gamma);

      if (gamma <= -1e-20)
	{
	  (*cla) << "LanczosSolve:Calc" << endl;
	  (*cla) << "preconditioner indefinite!" << endl;
	  (*cla) << "gamma " << gamma.Get(1) << endl;
	  exit(1);
	}

      gamma.Sqrt();
      scal.Div(pone, gamma);

      subd[iter+1] = gamma.Get(1);

      sol.Scal(scal);
      s->Scal(scal);

      sysmat.Mult(sol, *d);
    
      gamma *= -1;

      d->Add(gamma, rhs, pone, *d);
      rhs.Add(pone, *s, zero, rhs);

      gamma *= (1./diag[iter]);
      gamma.Abs();
      err   *= gamma.Get(1);
    }

  (*cla) << endl;

  ////////////////////////////////////////////////////////////////////

  if (iter < step)
    {
      CalcEigenValues();
    }
  else
    {
      (*cla) << "maxnumer iter too large!!!!" << endl;
    }
}

void LanczosSolver :: CalcEigenValues()
{
#ifdef TRACE
  (*trace) << "entering LanczosSolver::CalcEigenValues" << endl;
#endif

  Integer i,j,zbelow;
  Double xl, xu, x, q,lmin,lmax;

  for (j=1; j<=iter; j++)
    {
      xu = 0;
      
      for (i=1; i<=iter; i++)
	{
	  q = fabs(diag[i]) + fabs(subd[i]) + ((i < iter) ? fabs(subd[i+1]) : 0);
	  
	  if (q > xu)
	    {
	      xu = q;
	    }
	}
      
      xl = -xu;
      
      while (xu-xl > 1e-15 * fabs(xu) && xu-xl > 1e-100)
	{
	  x      = 0.5*(xu+xl);
	  zbelow = 0;
	  q      = 1;
	  
	  for (i=1; i<=iter; i++)
	    {
	      if (fabs(q) > 1e-100)
		{
		  q = diag[i]-x-subd[i]*subd[i]/q;
		}
	      else
		{
		  q = diag[i]-x-fabs(subd[i])*1e100;
		}
	      
	      if (q < 0)
		{
		  zbelow++;
		}
	    }
	  
	  if (zbelow < j)
	    {
	      xl = x;
	    }
	  else
	    {
	      xu = x;
	    }
	}

      if (j == 1)
	{
	  lmin = 0.5*(xl+xu);
	}

      if (j == iter)
	{
	  lmax = 0.5*(xl+xu);
	}
      
      (*cla) << "+++ eigenvalue number " << j << " " << 0.5*(xl+xu) << endl;
    }
  (*cla) << endl;
  (*cla) << "condition number " << lmax/lmin << endl;

}

void LanczosSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering LanczosSolver::CalcRate" << endl;
#endif


}

}
