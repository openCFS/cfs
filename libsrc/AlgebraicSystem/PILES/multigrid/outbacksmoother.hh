#ifndef FILE_OUTBACKSMOOTHER_CLA
#define FILE_OUTBACKSMOOTHER_CLA

namespace CoupledField
{

class OutbackSmoother : public BaseSmoother
{
public:
  ///
  OutbackSmoother(BaseParameter & param, Integer asize, Integer anumrhs);

  ///
  virtual ~OutbackSmoother();

  ///
  virtual void Calc(BaseMatrix & amat);

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level);

};

}

#endif // FILE_OUTBACKSMOOTHER_CLA
