#ifndef FILE_BASEPRECOND_PILES
#define FILE_BASEPRECOND_PILES

namespace CoupledField
{

class BasePrecond
{
public:
  ///
  BasePrecond();

  ///
  virtual ~BasePrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level=0) = 0;

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat) = 0;

  ///
  void InitPrecond()
    {mathier->Init();};

protected:
  ///
  Integer size, dof, numrhs, offset;

  ///
  Double omega;

  ///
  BaseHierarchy * mathier;

};

}

#endif // FILE_BASEPRECOND_PILES
