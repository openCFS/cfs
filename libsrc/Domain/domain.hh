#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

//#include "interface_gridcfs.hh"
//#include "interface_netgen.hh"
#include "grid.hh"
#include "basepde.hh"

namespace CoupledField
{

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

  //! set dirichlet boundary condition
//  void SetDBC(Integer apde, Integer level, Integer update);

  //! print solution
//  void PrintSolution(Double * sol, Integer apde);

  //! get pde
  BasePDE * GetPDE(const Integer ipde){ return ptpde_[ipde];}

  //! get algebraic system
  AbstractAlgebraicSys * GetAlgSys(){ return ptalgsys_;}

  //! get pointer to input-file
  FileType * GetInFile(){ return InFile_;}

  //! get pointer to output-file
  WriteResults<Dim> * GetOutFile(){ return OutFile_;}

  //! get pointer to boundary condition
  BCs * GetBCs(){ return ptBCs_;}

protected:

private:
   //! initialize pde
   void InitPDE();

   //! initialization of alg.sys.
   void InitAlgSys();

  //!
  Integer numsubdomain_;

  //!
  Integer numpde_;

  //!
  Integer numsys_;

  //!
  Integer numgraph_;

  //! 
  Integer ** syscoupling_;

  //!
  BasePDE * ptpde_[20];

  //!
//  Grid<Point2D> * ptgrid_;
  Grid<Dim> * ptgrid_;
  

  //!
  BCs * ptBCs_;

  //!
  AbstractAlgebraicSys * ptalgsys_;
  
  //!
  Material * ptmaterial_;

  //!
  TimeFunc * ptTimeFunc_;

  //!
  FileType *InFile_;

  //!
  WriteResults<Dim> * OutFile_;

};

#ifdef __GNUC__
template class Domain<Point3D>;
//template class Domain<Point2D>;
#endif

}

#endif // FILE_DOMAIN
