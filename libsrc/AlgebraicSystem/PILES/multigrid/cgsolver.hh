#ifndef FILE_CGSOLVER
#define FILE_CGSOLVER

namespace CoupledField
{

class CGSolver : public BaseSolver
{
public:
  ///
  CGSolver(BaseParameter & param, Integer size, Integer anumrhs, Integer dof);

  ///
  virtual ~CGSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol);

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat);

private:
  ///
  BaseVector * r, * s, * d;

  ///
  DenseVector * aalpha, * abeta, * agamma, * awra, * awrn, * aerr, * arel, *ascal, *apone, *amone;
 
  ///
  Integer maxiter, numrhs;

  ///
  Double eps, epsmach;
};

class LanczosSolver : public BaseSolver
{
public:
  ///
  LanczosSolver(BaseParameter & param, Integer size, Integer anumrhs, Integer dof);

  ///
  virtual ~LanczosSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol);

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat);

private:
  ///
  void CalcEigenValues();

  ///
  BaseVector * s, * d;

  ///
  DenseVector *agamma, * aerr, *ascal, *apone, *azero;
 
  ///
  Integer maxiter, numrhs,step;

  ///
  Double eps, epsmach,prec;

  ///
  Double * diag, *subd;
};

}

#endif // FILE_CGSOLVER
