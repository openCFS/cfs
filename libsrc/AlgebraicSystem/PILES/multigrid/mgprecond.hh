#ifndef FILE_MGPRECOND_PILES
#define FILE_MGPRECOND_PILES

namespace CoupledField
{

class BaseMGPrecond : public BasePrecond
{
public:
  ///
  BaseMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof=1);

  ///
  virtual ~BaseMGPrecond();

  ///
  virtual void Step(BaseVector & rhs, BaseVector & sol, Integer level=0);

  ///
  virtual void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat)
    {mathier->Calc(sysmat, auxmat);};

protected:
  ///
  Integer cycle[20];
  
};

//////////////////////////////////////////////////////////////////////

class RScalarMGPrecond : public BaseMGPrecond
{
public:
  ///
  RScalarMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RScalarMGPrecond();
};

////////////////////////////////////////////////////////////////////////

class RBlockMGPrecond : public BaseMGPrecond
{
public:
  ///
  RBlockMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof);

  ///
  virtual ~RBlockMGPrecond();
};

////////////////////////////////////////////////////////////////////////

class RFullMGPrecond : public BaseMGPrecond
{
public:
  ///
  RFullMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RFullMGPrecond();
};


///////////////////////////////////////////////////////////////////////

class CScalarMGPrecond : public BaseMGPrecond
{
public:
  ///
  CScalarMGPrecond(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~CScalarMGPrecond();
};

}

#endif // FILE_MGPRECOND_PILES
