#ifndef FILE_FULLSMOOTHER_PILES
#define FILE_FULLSMOOTHER_PILES

namespace CoupledField
{

class RFullGSSmoother : public BaseSmoother
{
public:
  ///
  RFullGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RFullGSSmoother();

  ///
  virtual void Calc(BaseMatrix & amat);

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level);

};

class RFullFastGSSmoother : public RFullGSSmoother
{
public:
  ///
  RFullFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~RFullFastGSSmoother();

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);
};

}

#endif // FILE_FULLSMOOTHER_PILES
