#ifndef FILE_BASEHIERARCHY_PILES
#define FILE_BASEHIERARCHY_PILES

namespace CoupledField
{
  
  ///
  struct HierarchyNode
  {
    BaseMatrix * sysmat;
    BaseMatrix * auxmat;
    BaseVector * vec1;
    BaseVector * vec2;
    BaseVector * vec3;
    BaseVector * vec4;
    BaseTransfer * syspro;
    BaseTransfer * sysres;
    BaseTransfer * auxpro;
    BaseTransfer * auxres;
    BaseCoarse * coarse;
    BaseSmoother * smooth;
  };

  /////////////////////////////////////////////////////////

class BaseHierarchy
{
public:
  ///
  BaseHierarchy(BaseParameter & aparam, Integer asize, Integer aoffset, Integer anumrhs, Integer adof=1);

  ///
  virtual ~BaseHierarchy();

  ///
  void Calc(BaseMatrix & sysmat, BaseMatrix & auxmat);

  ///
  virtual void Init() = 0;

  ///
  HierarchyNode * GetNode(Integer l) {return &node[l];};

  ///
  Integer GetSize() const {return (amglevel+gmglevel);};

  ///
  Integer GetAMGSize() const {return amglevel;};

  ///
  virtual void DeleteAuxMemory() = 0;

  ///
  virtual void MakeSmoother() = 0;

  ///
  virtual void MakeCoarsening() = 0;

  ///
  virtual void MakeSysTransfer() = 0;

  ///
  virtual void MakeAuxTransfer() = 0;

  ///
  virtual void MakeSysMatrix() = 0;

  ///
  virtual void MakeAuxMatrix() = 0;

  ///
  virtual void MakeVector() = 0;

protected:
  ///
  HierarchyNode node[20];

  ///
  BaseParameter * param;

  ///
  Integer amglevel, gmglevel, offset, coarsesys, level, numrhs, dof, * rsw;

  ///
  Integer fsize_sys, fsize_aux, csize_sys, csize_aux;

  ///
  Integer fnne_sys, fnne_aux, cnne_sys, cnne_aux, dir_aux, dir_sys;

  ///
  Double alpha[20], epsmat, cputime;

  ///
  CPUClock cpuclock;
};

}

#endif // FILE_BASEHIERARCHY_PILES
