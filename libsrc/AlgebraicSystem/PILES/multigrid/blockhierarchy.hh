#ifndef FILE_BLOCKHIERARCHY_PILES
#define FILE_BLOCKHIERARCHY_PILES

namespace CoupledField
{
class RBlockHierarchy : public BaseHierarchy
{
public:
  ///
  RBlockHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs, Integer adof);

  ///
  virtual ~RBlockHierarchy();

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

#endif // FILE_BLOCKHIERARCHY_PILES
