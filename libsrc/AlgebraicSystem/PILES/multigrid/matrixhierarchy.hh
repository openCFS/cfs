#ifndef FILE_BASEPRECOND
#define FILE_BASEPRECOND

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
  virtual void Step(BaseVector & rhs, BaseVector & sol) = 0;

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat) = 0;

protected:
  Integer size, dof, numrhs;

  ///
  MatrixHierarchy * mathier;

};

}

#endif // FILE_BASEPRECOND
