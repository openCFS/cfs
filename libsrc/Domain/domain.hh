#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include "interface_gridcfs.hh"
#include "outUnverg.hh"
#include "basepde.hh"
#include "material.hh"
#include "matrix.hh"
#include "abstractAlgSys.hh"

namespace CoupledField
{
class FileType;


/// contain information about calculation domain and according to different meshes create different grids

template<class Dim>
class Domain
{
public:
  //!
  Domain(FileType * const aptFileType, OutResultUnverg<Dim> * ptUnverg, Material * materialdata,
         Grid<Dim> * ptgrid);

  //!
  virtual ~Domain();

  //!
  void InitPDE();

  //!
  void InitAlgSys(AbstractAlgebraicSys * algsys);

  //!
  void PrintGrid(Integer level);

  //!
  void SetSubdomains();

  //!
  void PrintDomain();

  //!
  BasePDE * ptpde[20];

  //!
  AbstractAlgebraicSys * ptalgsys;

protected:

private:
  //!
  Integer numsubdomain;

  //!
  Integer numpde;

  //!
  Integer numsys;

  //!
  Integer numgraph;

  //! 
  Integer ** syscoupling;

  //!
  //  BasePDE * ptpde[20];

  //!
  Grid<Dim> * grid;

  //!
  //  AbstractAlgebraicSys * ptalgsys;
  
  //!
  Material * ptmaterial;

  //!
  FileType *InFile;

  //!
  OutResultUnverg<Dim> *OutFile;

};

#ifdef __GNUC__
  //template class Domain<Point3D>;
template class Domain<Point2D>;
#endif

}

#endif // FILE_DOMAIN
