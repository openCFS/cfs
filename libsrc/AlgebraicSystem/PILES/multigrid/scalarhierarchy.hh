#ifndef FILE_SCALARHIERARCHY_PILES
#define FILE_SCALARHIERARCHY_PILES

namespace CoupledField
{
class RScalarHierarchy : public BaseHierarchy
{
public:
  ///
  RScalarHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs);

  ///
  virtual ~RScalarHierarchy();

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

class CScalarHierarchy : public BaseHierarchy
{
public:
  ///
  CScalarHierarchy(BaseParameter & param, Integer asize, Integer aoffset, Integer anumrhs);

  ///
  virtual ~CScalarHierarchy();

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

#endif // FILE_SCALARHIERARCHY_PILES
