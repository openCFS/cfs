#ifndef FILE_OUTBACKHIERARCHY_PILES
#define FILE_OUTBACKHIERARCHY_PILES

namespace CoupledField
{
class OutbackHierarchy : public BaseHierarchy
{
public:
  ///
  OutbackHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs);

  ///
  virtual ~OutbackHierarchy();

  ///
  virtual void Init();

  ///
  virtual void DeleteAuxMemory();

  ///
  virtual void MakeSmoother();

  ///
  virtual void MakeCoarsening();

  ///
  virtual void MakeSysTransfer();

  ///
  virtual void MakeAuxTransfer();

  ///
  virtual void MakeSysMatrix();

  ///
  virtual void MakeAuxMatrix();

  ///
  virtual void MakeVector();
};

}

#endif // FILE_OUTBACKHIERARCHY_PILES
