#ifndef FILE_FULLHIERARCHY_CLA
#define FILE_FULLHIERARCHY_CLA

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

#endif // FILE_FULLHIERARCHY_CLA
