#ifndef FILE_BASESMOOTHER_PILES
#define FILE_BASESMOOTHER_PILES

namespace CoupledField
{

class BaseSmoother
{
public:
  ///
  BaseSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof=1);

  ///
  virtual ~BaseSmoother();

  ///
  virtual void Calc(BaseMatrix & amat) = 0;

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level) = 0;

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level) = 0;

protected:
  ///
  Integer size, dof, numrhs;

  ///
  Integer numsmoothfor, numsmoothback;
};

}

#endif // FILE_BASESMOOTHER_PILES
