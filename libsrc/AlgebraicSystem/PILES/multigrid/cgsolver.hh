#ifndef FILE_CGSOLVER_CLA
#define FILE_CGSOLVER_CLA

namespace CoupledField
{

//! Solver class for Conjugate Gradient Method
class CGSolver : public BaseSolver
{
public:
  //! Constructor
  /*! 
      \param param solver parameters
      \param size number of equations
      \param anumrhs number of right hand sides
      \param dof number of dofs per node (edge)
  */
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

#endif // FILE_CGSOLVER_CLA
