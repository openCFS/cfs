#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include "grid.hh"
#include "bcs.hh"
#include "basepde.hh"

namespace CoupledField
{

/// contain information about calculation domain and according to different meshes create different grids

class Domain
{
public:
  //!
  Domain(FileType * const aptFileType, WriteResults * ptOut, Material * materialdata,  TimeFunc * aptTimeFunc);

  //!
  virtual ~Domain();

  //!
  void PrintGrid(Integer level);

  //!
  void SetSubdomains();

  //! get pde
  BasePDE * GetPDE(const Integer ipde){ return ptpde_[ipde];}

  //! get algebraic system
  AbstractAlgebraicSys * GetAlgSys(){ return ptalgsys_;}

  //! get pointer to input-file
  FileType * GetInFile(){ return InFile_;}

  //! get pointer to output-file
  WriteResults * GetOutFile(){ return OutFile_;}

  //! get pointer to boundary condition
  BCs * GetBCs(){ return ptBCs_;}

  //! update algebraic system and bcs after refinement the mesh
  void Update();

  //!
  Grid * GetGrid(){ return ptgrid_;}

 void TestGrid();
protected:

private:
   //! initialize pde
   void InitPDE();

   //! initialization of alg.sys.
   void InitAlgSys(const Integer level);

  //! update alg. sys. in case of new mesh
  void UpdateAlgSys();

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
  Grid * ptgrid_;
  

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
  WriteResults * OutFile_;

};

}

#endif // FILE_DOMAIN
