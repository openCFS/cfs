#ifndef FILE_JACPRECOND_PILES
#define FILE_JACPRECOND_PILES

namespace CoupledField
{

class RJACScalarPrecond : public BasePrecond
{
public:
  ///
  RJACScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RJACScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

private:
  ///
  RScalarMatrix * mat;

};

class CJACScalarPrecond : public BasePrecond
{
public:
  ///
  CJACScalarPrecond(Integer asize, Integer anumrhs);

  ///
  virtual ~CJACScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

private:
  ///
  CScalarMatrix * mat;

  ///
  Double * val;

};

}

#endif // FILE_JACPRECOND_PILES
