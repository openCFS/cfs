#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include "interface_gridcfs.hh"
#include "basepde.hh"

namespace CoupledField
{
//class FileType;


/// contain information about calculation domain and according to different meshes create different grids

template<class Dim>
class Domain
{
public:
  //!
  Domain(FileType * const aptFileType, WriteResults<Dim> * ptOut, Material * materialdata,  TimeFunc * aptTimeFunc);

  //!
  virtual ~Domain();

  //!
  void PrintGrid(Integer level);

  //!
  void SetSubdomains();

  //!
  void PrintDomain();

  //! set dirichlet boundary condition
  void SetDBC(Integer apde, Integer level, Integer update);

  //! print solution
  void PrintSolution(Double * sol, Integer apde);

  //! get pde
  BasePDE * GetPDE(const Integer ipde){ return ptpde[ipde];}

  //! get algebraic system
  AbstractAlgebraicSys * GetAlgSys(){ return ptalgsys;}

  //! get pointer to input-file
  FileType * GetInFile(){ return InFile;}

  //! get pointer to output-file
  WriteResults<Dim> * GetOutFile(){ return OutFile;}

  //! get pointer to boundary condition
  BCs * GetBCs(){ return ptBCs;}

protected:

private:
   //! initialize pde
   void InitPDE();

   //! initialization of alg.sys.
   void InitAlgSys();

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
  BasePDE * ptpde[20];

  //!
  Grid<Dim> * ptgrid;

  //!
  BCs * ptBCs;

  //!
  AbstractAlgebraicSys * ptalgsys;
  
  //!
  Material * ptmaterial;

  //!
  TimeFunc * ptTimeFunc;

  //!
  FileType *InFile;

  //!
  WriteResults<Dim> * OutFile;

};

#ifdef __GNUC__
  //template class Domain<Point3D>;
template class Domain<Point2D>;
#endif

}

#endif // FILE_DOMAIN
