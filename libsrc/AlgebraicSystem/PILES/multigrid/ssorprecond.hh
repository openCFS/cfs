#ifndef FILE_SSORPRECOND_CLA
#define FILE_SSORPRECOND_CLA

namespace CoupledField
{

class RSSORScalarPrecond : public BasePrecond
{
public:
  ///
  RSSORScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RSSORScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level=0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

private:
  ///
  RScalarMatrix * mat;
};



class CSSORScalarPrecond : public BasePrecond
{
public:
  ///
  CSSORScalarPrecond(Integer asize, Integer anumrhs);

  ///
  virtual ~CSSORScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level=0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

private:
  ///
  CScalarMatrix * mat;

};

}

#endif // FILE_SSORPRECOND_CLA
