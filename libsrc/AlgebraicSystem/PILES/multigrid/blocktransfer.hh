#ifndef FILE_BLOCKTRANSFER_PILES
#define FILE_BLOCKTRANSFER_PILES

namespace CoupledField
{

class RBlockTransfer : public BaseTransfer
{
public:
  ///
  RBlockTransfer(Integer asize, Integer adof, Integer * rsw);

  ///
  virtual ~RBlockTransfer();

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
  DenseMatrix * Get(Integer i, Integer j) const 
    {*dm = &val[(start[i-1]+j-1)*dof*dof]; return dm;};

private:
  ///
  DenseMatrix * dm;

  ///
  Double * sum;
};

}

#endif // FILE_BLOCKTRANSFER_PILES
