#ifndef FILE_RSOLVER
#define FILE_RSOLVER

namespace CoupledField
{

class RRSolver : public BaseSolver
{
public:
  ///
  RRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof);

  ///
  virtual ~RRSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol);

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat);

private:
  ///
  BaseVector * r, * s;

  ///
  DenseVector * err, * rel, * wra, * scal, * mone, * pone, * omega;

  ///
  Integer maxiter;

  ///
  Double eps, epsmach;
};


class CRSolver : public BaseSolver
{
public:
  ///
  CRSolver(BaseParameter & param, Integer size, Integer numrhs, Integer dof);

  ///
  virtual ~CRSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol);

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat);

private:
  ///
  BaseVector * r, * s;

  ///
  DenseVector * err, * rel, * wra, * scal, * mone, * pone, * omega;

  ///
  Integer maxiter;

  ///
  Double eps, epsmach;
};

}

#endif // FILE_RSOLVER
