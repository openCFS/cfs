#ifndef FILE_BASESOLVER
#define FILE_BASESOLVER

namespace CoupledField
{

class BaseSolver
{
public:
  ///
  BaseSolver();

  ///
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

#endif // FILE_BASESOLVER
