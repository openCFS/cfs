#ifndef FILE_IDPRECOND_PILES
#define FILE_IDPRECOND_PILES

namespace CoupledField
{

class RIDScalarPrecond : public BasePrecond
{
public:
  ///
  RIDScalarPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RIDScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

private:
  ///
  Integer size;

};

class CIDScalarPrecond : public BasePrecond
{
public:
  ///
  CIDScalarPrecond(Integer asize, Integer anumrhs);

  ///
  virtual ~CIDScalarPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat);

private:
  ///
  Integer size;

};



class RIDBlockPrecond : public BasePrecond
{
public:
  ///
  RIDBlockPrecond(Integer asize, Integer adof);

  ///
  virtual ~RIDBlockPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat);

private:
  ///
  Integer size, dof;

};

class CIDBlockPrecond : public BasePrecond
{
public:
  ///
  CIDBlockPrecond(Integer asize, Integer adof);

  ///
  virtual ~CIDBlockPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level = 0);

  ///
  virtual void Calc(BaseMatrix & sysmat);

private:
  ///
  Integer size, dof;

};

}

#endif // FILE_IDPRECOND_PILES
