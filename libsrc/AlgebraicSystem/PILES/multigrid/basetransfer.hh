#ifndef FILE_BASETRANSFER_CLA
#define FILE_BASETRANSFER_CLA

namespace CoupledField
{

class BaseTransfer
{
public:
  ///
  BaseTransfer(Integer asize, Integer adof=1);

  ///
  virtual ~BaseTransfer();

  ///
  virtual void Calc(BaseTopology * topology) = 0;

  ///
  virtual void MultHh(BaseVector & uH, BaseVector & uh) = 0;

  ///
  virtual void MulthH(BaseVector & uh, BaseVector & uH) = 0;

  ///
  virtual void Print() const = 0;

  ///
  void Calculated() {calculated = TRUE;};

  ///
  Boolean IsCalculated() const {return calculated;};

protected:
  ///
  Double * val;

  ///
  Integer size, dof;

  ///
  Integer * start, * pos;

  ///
  Boolean calculated;
};

}

#endif // FILE_BASETRANSFER_CLA
