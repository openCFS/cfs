#ifndef FILE_BASESOLVER_CLA
#define FILE_BASESOLVER_CLA

namespace CoupledField
{

  //! Base class for algebraic system solver
class BaseSolver
{
public:
  //! Constructor
  BaseSolver();

  //! Destructor
  virtual ~BaseSolver();

  ///
  virtual void Calc(BaseMatrix & sysmat, BasePrecond & premat,
		    BaseVector & rhs, BaseVector & sol) = 0;

  ///
  virtual void CalcRate(BaseMatrix & sysmat, BasePrecond & premat) = 0;

  ///
  Integer GetNumIter() const {return iter;};

protected:
  ///
  Integer iter;
};

}

#endif // FILE_BASESOLVER_CLA
