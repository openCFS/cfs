#ifndef FILE_BASEPRECOND_CLA
#define FILE_BASEPRECOND_CLA

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

#endif // FILE_BASEPRECOND_CLA
