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

RRSolver :: RRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof)
  : BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering RRSolver::RRSolver" << endl;
#endif

  r = new RealVector(size, numrhs, dof);
  s = new RealVector(size, numrhs, dof);

  

  err   = new DenseVector(numrhs);
  rel   = new DenseVector(numrhs);
  scal  = new DenseVector(numrhs);
  wra   = new DenseVector(numrhs);
  pone  = new DenseVector(numrhs);
  mone  = new DenseVector(numrhs);
  omega = new DenseVector(numrhs);

  eps     = param.GetEps();
  epsmach = param.GetEpsMach();
  maxiter = param.GetMaxNumIter();
  //omega[0]= param.GetDampIter();
}
  
RRSolver :: ~RRSolver()
{
#ifdef TRACE
  (*trace) << "entering RRSolver::~RRSolver" << endl;
#endif

  delete err;
  delete rel;
  delete scal;
  delete wra;
  delete pone;
  delete mone;
}

void RRSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		      BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering RRSolver::Calc" << endl;
#endif

//   Boolean loop;

//   pone[0] = +1;
//   mone[0] = -1;

//   iter = 1;
//   loop = TRUE;

//   // generate an initial guess

//   premat.Step(rhs, sol);

//   sysmat.Mult(sol, *s);
//   r->Add(pone, rhs, mone, *s);
//   s->Init();

//   premat.Step(*r, *s);
  
//   wra = s->Inner(*r);

//   scal[0] = wra[0];
//   rel[0]  = 1;
//   err[0]  = eps*eps*scal[0];
  
//   if ((wra[0] <= err[0]) || (sqrt(wra[0]) <= epsmach))
//     {
//       loop = FALSE;
//     }

//   (*cla) << 1 << " " << rel[0] << endl;
//   (*conv) << rel[0] << endl;

//   while (loop && (iter <= maxiter))
//     {
//       iter++;

//       sol.Add(omega, *s);
//       sysmat.Mult(sol, *s);
//       r->Add(pone, rhs, mone, *s);
//       s->Init();
      
//       premat.Step(*r, *s);
      
//       wra = s->Inner(*r);

//       if (wra[0] < 0)
// 	{
// 	  (*cla) << "... EXIT in PCG-SOLVER" << endl;
// 	  exit(1);
// 	}

//       if ((wra[0] <= err[0]) || (sqrt(wra[0]) <= epsmach))
// 	{
// 	  loop =FALSE;
// 	}

//       rel[0] = sqrt(wra[0]/scal[0]);
	      
//       (*cla) << iter << " " << rel[0] << endl;
//       (*conv) << rel[0] << endl;
//     }

  (*cla) << "... number of iterations = " << iter << endl;
}

void RRSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering RRSolver::CalcRate" << endl;
#endif
}


CRSolver :: CRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof)
  : BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering CRSolver::CRSolver" << endl;
#endif

  r = new ComplexVector(size, numrhs, dof);
  s = new ComplexVector(size, numrhs, dof);

  err   = new DenseVector(numrhs);
  rel   = new DenseVector(numrhs);
  scal  = new DenseVector(numrhs);
  wra   = new DenseVector(numrhs);
  pone  = new DenseVector(numrhs);
  mone  = new DenseVector(numrhs);
  omega = new DenseVector(numrhs);

  eps     = param.GetEps();
  epsmach = param.GetEpsMach();
  maxiter = param.GetMaxNumIter();
  //omega[0]= param.GetDampIter();
}
  
CRSolver :: ~CRSolver()
{
#ifdef TRACE
  (*trace) << "entering CRSolver::~CRSolver" << endl;
#endif

  delete err;
  delete rel;
  delete scal;
  delete wra;
  delete pone;
  delete mone;
}

void CRSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		      BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering CRSolver::Calc" << endl;
#endif

//   Boolean loop;

//   pone[0] = +1;
//   mone[0] = -1;

//   iter = 1;
//   loop = TRUE;

//   // generate an initial guess

//   premat.Step(rhs, sol);

//   sysmat.Mult(sol, *s);
//   r->Add(pone, rhs, mone, *s);
//   s->Init();

//   premat.Step(*r, *s);
  
//   wra = s->Inner(*r);

//   scal[0] = wra[0];
//   rel[0]  = 1;
//   err[0]  = eps*eps*scal[0];
  
//   if ((wra[0] <= err[0]) || (sqrt(wra[0]) <= epsmach))
//     {
//       loop = FALSE;
//     }

//   (*cla) << 1 << " " << rel[0] << endl;
//   (*conv) << rel[0] << endl;

//   while (loop && (iter <= maxiter))
//     {
//       iter++;

//       sol.Add(omega, *s);
//       sysmat.Mult(sol, *s);
//       r->Add(pone, rhs, mone, *s);
//       s->Init();
      
//       premat.Step(*r, *s);
      
//       wra = s->Inner(*r);

//       if (wra[0] < 0)
// 	{
// 	  (*cla) << "... EXIT in PCG-SOLVER" << endl;
// 	  exit(1);
// 	}

//       if ((wra[0] <= err[0]) || (sqrt(wra[0]) <= epsmach))
// 	{
// 	  loop =FALSE;
// 	}

//       rel[0] = sqrt(wra[0]/scal[0]);
	      
//       (*cla) << iter << " " << rel[0] << endl;
//       (*conv) << rel[0] << endl;
//     }

  (*cla) << "... number of iterations = " << iter << endl;
}

void CRSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering CRSolver::CalcRate" << endl;
#endif
}
}
