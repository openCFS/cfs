#ifndef FILE_BLOCKSMOOTHER_CLA
#define FILE_BLOCKSMOOTHER_CLA

namespace CoupledField
{

class RBlockGSSmoother : public BaseSmoother
{
public:
  ///
  RBlockGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof);

  ///
  virtual ~RBlockGSSmoother();

  ///
  virtual void Calc(BaseMatrix & amat);

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

  ///
  virtual void StepBackward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol,
			    Integer level);

protected:
  ///
  Double * sum;

  ///
  DenseMatrix * x;
};


class RBlockFastGSSmoother : public RBlockGSSmoother
{
public:
  ///
  RBlockFastGSSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof);

  ///
  virtual ~RBlockFastGSSmoother();

  ///
  virtual void StepForward(BaseMatrix & sys, BaseMatrix & aux, BaseVector & rhs, BaseVector & sol, 
			   BaseVector & def, BaseVector & res, Integer level);

};

}

#endif // FILE_BLOCKSMOOTHER_CLA
