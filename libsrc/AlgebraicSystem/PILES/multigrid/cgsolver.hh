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
  Integer maxiter, numrhs;

  ///
  Double eps, epsmach;
};

}

#endif // FILE_CGSOLVER
