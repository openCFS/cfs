#ifndef FILE_SCALARTRANSFER_CLA
#define FILE_SCALARTRANSFER_CLA

namespace CoupledField
{

class RScalarTransfer : public BaseTransfer
{
public:
  ///
  RScalarTransfer(Integer asize, Integer * rsw, Boolean mem=TRUE);

  ///
  virtual ~RScalarTransfer();

  ///
  virtual void Calc(BaseTopology * topology);

  ///
  virtual void MultHh(BaseVector & uH, BaseVector & uh);

  ///
  virtual void MulthH(BaseVector & uh, BaseVector & uH);

  ///
  virtual void Print() const;

  ///
  Integer GetRowSize(Integer p) const {return (start[p]-start[p-1]);};

  ///
  Integer Pos(Integer p, Integer q) const {return pos[start[p-1]+q-1];};

  ///
  Double Get(Integer p, Integer q) const {return val[start[p-1]+q-1];};
};

}

#endif // FILE_SCALARTRANSFER_CLA
