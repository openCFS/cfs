#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

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
}
  
CGSolver :: ~CGSolver()
{
#ifdef TRACE
  (*trace) << "entering CGSolver::~CGSolver" << endl;
#endif

  delete r;
  delete s;
  delete d;
}

void CGSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		      BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering CGSolver::Calc" << endl;
#endif

  cout << "### preconditioned CG Solver" << endl;

  DenseVector alpha(numrhs);
  DenseVector beta(numrhs);
  DenseVector gamma(numrhs);
  DenseVector wra(numrhs);
  DenseVector wrn(numrhs);
  DenseVector err(numrhs);
  DenseVector rel(numrhs);
  DenseVector scal(numrhs);
  DenseVector pone(numrhs);
  DenseVector mone(numrhs);

  r->Init();
  s->Init();
  d->Init();

  Boolean loop;

  pone = +1.0;
  mone = -1.0;

  iter = 1;
  loop = TRUE;

  // generate an initial guess

  //premat.Step(rhs, sol);

  sysmat.Mult(sol, *d);

  //sysmat.Print();

  r->Add(pone, rhs, mone, *d);

  premat.Step(*r, *s);

  wra  = *s->Inner(*r);

  scal = wra;
  rel  = 1;
  err  = wra;
  err *= (eps*eps);

  if ((wra <= err) || (wra <= (epsmach*epsmach)))
    {
      loop = FALSE;
    }

  cout << 1 << " " << rel.Get(1) << " " << wra.Get(1) << endl;
  //(*conv) << rel << endl;

  while (loop && (iter <= maxiter))
    {
      iter++;
      sysmat.Mult(*s, *d);
      gamma = *d->Inner(*s);
      alpha = wra/gamma;

      sol.Add(alpha, *s);

      alpha *= -1;

      r->Add(alpha, *d);
      d->Init();
      
      premat.Step(*r, *d);

      wrn = *d->Inner(*r);
      
      if ((wrn <= err) || (wrn <= (epsmach*epsmach)))
	{
	  loop =FALSE;
	}

      beta = wrn/wra;
      wra  = wrn;
      s->Add(beta, *s,pone, *d);

      if (wra <= -1e-32)
	{
	  cout << "... EXIT in PCG-SOLVER " << wra.Get(1) << endl;
	  exit(1);
	}

      rel = wra/scal;
      rel.Sqrt();
	      
      cout << iter << " " << rel.Get(1) << " " << wra.Get(1) << endl;
      //(*conv) << rel << endl;
    }

  cout << "... number of iterations = " << iter << endl;
}

void CGSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering CGSolver::CalcRate" << endl;
#endif
}
}
