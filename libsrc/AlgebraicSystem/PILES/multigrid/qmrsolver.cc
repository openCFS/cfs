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

QMRSolver :: QMRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof)
  : BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering QMRSolver::QMRSolver" << endl;
#endif

  r = new ComplexVector(size, numrhs, dof);
  s = new ComplexVector(size, numrhs, dof);
  d = new ComplexVector(size, numrhs, dof);
  g = new ComplexVector(size, numrhs, dof);
  h = new ComplexVector(size, numrhs, dof);

  beta    = new DenseVector(numrhs);
  gamma   = new DenseVector(numrhs);
  eta     = new DenseVector(numrhs);
  delta   = new DenseVector(numrhs);
  epsilon = new DenseVector(numrhs);
  omega   = new DenseVector(numrhs);
  zeta    = new DenseVector(numrhs);

  r0        = new DenseVector(numrhs);
  relnorm   = new DenseVector(numrhs);
  rho_new   = new DenseVector(numrhs);
  rho_old   = new DenseVector(numrhs);
  theta_new = new DenseVector(numrhs);
  theta_old = new DenseVector(numrhs);
  c_new     = new DenseVector(numrhs);
  c_old     = new DenseVector(numrhs);
  znorm     = new DenseVector(numrhs);

  eps     = param.GetEps();
  epsmach = param.GetEpsMach();
  maxiter = param.GetMaxNumIter();
}
  
QMRSolver :: ~QMRSolver()
{
#ifdef TRACE
  (*trace) << "entering QMRSolver::~QMRSolver" << endl;
#endif

  delete beta;
  delete gamma;
  delete eta;
  delete delta;
  delete epsilon;
  delete omega;
  delete zeta;

  delete r0;
  delete relnorm;
  delete rho_new;
  delete rho_old;
  delete theta_new;
  delete theta_old;
  delete c_new;
  delete c_old;
  delete znorm;
}

void QMRSolver :: Calc(BaseMatrix & sysmat, BasePrecond & premat, 
		       BaseVector & rhs, BaseVector & sol)
{
#ifdef TRACE
  (*trace) << "entering QMRSolver::Calc" << endl;
#endif

//   Double re,im;
//   Boolean loop;

//   iter         = 0;
//   c_old[0]     = 1;
//   epsilon[0]   = 1;
//   epsilon[1]   = 0;
//   theta_old[0] = 0;
//   eta[0]       = -1;
//   eta[1]       = 0;
//   loop         = TRUE;

//   (*test) << "rhs" << endl;
//   rhs.Print();
//   (*test) << "sol"  << endl;
//   sol.Print();

//   sysmat.Mult(sol, *g);

//   (*test) << "g" << endl;
//   g->Print();

//   omega[0] = 1;
//   omega[1] = 0;
//   zeta[0]  = -1;
//   zeta[1]  = 0;

//   r->Add(omega, rhs, zeta,* g);

//   (*test) << "r" << endl;
//   r->Print();
 
//   r0         = r->L2Norm();
//   rho_old[0] = r0[0];
//   relnorm[0] = 1; 

//   omega[0] = 1./rho_old[0];
//   omega[1] = 0;
//   zeta[0]  = 0;
//   zeta[1]  = 0;

//   g->Add(omega, *r, zeta, *r);

//   //s->Init();
//   //d->Init();

//   znorm[0] = sqrt(epsilon[0]*epsilon[0]+epsilon[1]*epsilon[1]);
 
//   while (relnorm[0] > eps && iter <= maxiter && znorm[0] > epsmach && loop == TRUE)
//     {
//       h->Init();
//       premat.Step(*g, *h);

//       (*test) << "g" << endl;
//       g->Print();
//       (*test) << "h" << endl;
//       h->Print();

//       delta = h->Inner(*g);

//       (*test) << "delta" << endl;
//       (*test) << delta[0] << " " << delta[1] << endl;

//       if (sqrt(delta[0]*delta[0]+delta[1]*delta[1]) > epsmach)
// 	{
// 	  (*test) << "in the if case" << endl; 

// 	  znorm[0] = epsilon[0]*epsilon[0]+epsilon[1]*epsilon[1];

// 	  gamma[0] = -rho_old[0]*(delta[0]*epsilon[0]+delta[1]*epsilon[1])/znorm[0];
// 	  gamma[1] = -rho_old[0]*(epsilon[0]*delta[1]-delta[0]*epsilon[1])/znorm[0];

// 	  zeta[0]  = 1;
// 	  zeta[1]  = 0;

// 	  (*test) << "s" << endl;
// 	  s->Print();

// 	  s->Add(zeta, *h, gamma, *s);

// 	  sysmat.Mult(*s, *h);

// 	  (*test) << "s" << endl;
// 	  s->Print();
// 	  (*test) << "h" << endl;
// 	  h->Print();

// 	  epsilon = s->Inner(*h);

// 	  (*test) << "epsilons" << endl;
// 	  (*test) << epsilon[0] << " " << epsilon[1] << endl;

// 	  znorm[0]= delta[0]*delta[0]+delta[1]*delta[1];
// 	  beta[0] = -(epsilon[0]*delta[0]+epsilon[1]*delta[1])/znorm[0];
// 	  beta[1] = -(delta[0]*epsilon[1]-epsilon[0]*delta[1])/znorm[0];

// 	  zeta[0]  = 1;
// 	  zeta[1]  = 0;

// 	  g->Add(zeta, *h, beta, *g);
	  
// 	  rho_new      = g->L2Norm();
// 	  theta_new[0] = rho_new[0]/(c_old[0]*sqrt(beta[0]*beta[0]+beta[1]*beta[1]));
// 	  c_new[0]     = 1./(sqrt(1+theta_new[0]*theta_new[0]));
// 	  znorm[0]     = beta[0]*beta[0]+beta[1]*beta[1];

// 	  (*test) << "beta, znorm, c_new, c_old, rho_old, eta" << endl;
// 	  (*test) << beta[0] << " " << beta[1] << endl;
// 	  (*test) << znorm[0] << endl;
// 	  (*test) << c_new[0] << endl;
// 	  (*test) << c_old[0] << endl;
// 	  (*test) << rho_old[0] << endl;
// 	  (*test) << eta[0] << " " << eta[1] << endl;

// 	  re     = (rho_old[0]*c_new[0]*c_new[0]/(c_old[0]*c_old[0]))*(eta[0]*beta[0]+eta[1]*beta[1])/znorm[0];
// 	  im     = (rho_old[0]*c_new[0]*c_new[0]/(c_old[0]*c_old[0]))*(beta[0]*eta[1]-eta[0]*beta[1])/znorm[0];
// 	  eta[0] = re;
// 	  eta[1] = im;

// 	  gamma[0]  = theta_old[0]*theta_old[0]*c_new[0]*c_new[0];
// 	  gamma[1]  = 0;

// 	  (*test) << "s" << endl;
// 	  s->Print();
// 	  (*test) << "d" << endl;
// 	  d->Print();

// 	  (*test) << "eta,gamma" << endl;
// 	  (*test) << eta[0] << " " << eta[1] << endl;
// 	  (*test) << gamma[0] << " " << gamma[1] << endl;

// 	  d->Add(eta, *s, gamma, *d);

// 	  zeta[0]  = 1;
// 	  zeta[1]  = 0;
	  
// 	  (*test) << "d" << endl;
// 	  d->Print();

// 	  sol.Add(zeta, *d);

// 	  if (rho_new[0] < epsmach)
// 	    {
// 	      loop = FALSE;
// 	    }
// 	  else
// 	    {
// 	      zeta[0]  = 1./rho_new[0];
// 	      zeta[1]  = 0;
// 	      omega[0] = 0;
// 	      omega[1] = 0;

// 	      g->Add(zeta, *g, omega, *g);

// 	      rho_old[0]   = rho_new[0];
// 	      theta_old[0] = theta_new[0];
// 	      c_old[0]     = c_new[0];	      
// 	    }
// 	}
//       else
// 	{
// 	  loop = FALSE;
// 	}
      
//       // compute the residul in the l2 norm
      
//       (*test) << "sol" << endl;
//       sol.Print();

//       sysmat.Mult(sol ,*h);

//       (*test) << "h" << endl;
//       h->Print();

//       zeta[0]  = 1;
//       zeta[1]  = 0;
//       omega[0] = -1;
//       omega[1] = 0;

//       r->Add(zeta, rhs, omega, *h);

      
//       (*test) << "residual" << endl;
//       r->Print();

//       relnorm  = r->L2Norm();
//       znorm[0] = sqrt(epsilon[0]*epsilon[0]+epsilon[1]*epsilon[1]);

//       iter++;

//       if (iter == 1)
// 	{
// 	  r0[0] = relnorm[0];
// 	}

//       relnorm[0] = relnorm[0]/r0[0];
      
//       (*cla) << iter << " " << relnorm[0] << endl;
//       (*conv) << iter << " " << relnorm[0] << endl; 
//     }

  (*cla) << "number of iterations " << iter << endl;
}

void QMRSolver :: CalcRate(BaseMatrix & sysmat, BasePrecond & premat)
{
#ifdef TRACE
  (*trace) << "entering QMRSolver::CalcRate" << endl;
#endif
}
}
