#ifndef FILE_QMRSOLVER_CLA
#define FILE_QMRSOLVER_CLA

namespace CoupledField
{

class QMRSolver : public BaseSolver
{
public:
  ///
  QMRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof);

  ///
  virtual ~QMRSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol);

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat);

private:
  ///
  BaseVector * r, * s, * d, *g, *h;

  ///
  DenseVector *beta, *eta, *gamma, *epsilon, *omega, *delta, *zeta;

  ///
  DenseVector *r0, *relnorm, *rho_new, *rho_old, *theta_new, *theta_old, *c_new, *c_old, *znorm;

  ///
  Integer maxiter;

  ///
  Double eps, epsmach;
};

}

#endif // FILE_QMRSOLVER_CLA
