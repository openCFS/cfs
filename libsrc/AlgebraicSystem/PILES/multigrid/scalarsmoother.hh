#ifndef FILE_SCALARSMOOTHER_CLA
#define FILE_SCALARSMOOTHER_CLA

namespace CoupledField
{

class RScalarGSSmoother : public BaseSmoother
{
public:
  ///
  RScalarGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RScalarGSSmoother();

  ///
  virtual void Calc(BaseMatrix & amat);

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level);

};

class RScalarFastGSSmoother : public RScalarGSSmoother
{
public:
  ///
  RScalarFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RScalarFastGSSmoother();

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);
};


///////////////////////////////////////// complex ///////////////////////////////////////////

class CScalarGSSmoother : public BaseSmoother
{
public:
  ///
  CScalarGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~CScalarGSSmoother();

  ///
  virtual void Calc(BaseMatrix & amat);

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level);

};

class CScalarFastGSSmoother : public CScalarGSSmoother
{
public:
  ///
  CScalarFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~CScalarFastGSSmoother();

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);
};

}

#endif // FILE_SCALARSMOOTHER_CLA
