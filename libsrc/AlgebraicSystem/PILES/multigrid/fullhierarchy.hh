#ifndef FILE_FULLHIERARCHY_PILES
#define FILE_FULLHIERARCHY_PILES

namespace CoupledField
{
class RFullHierarchy : public BaseHierarchy
{
public:
  ///
  RFullHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs);

  ///
  virtual ~RFullHierarchy();

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

#endif // FILE_FULLHIERARCHY_PILES
