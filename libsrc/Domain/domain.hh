#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#define MAXNUMPDE 20

#include "grid.hh"
#include "bcs.hh"
#include <PDE/basepde.hh>

namespace CoupledField
{

//! contain information about calculation domain and according to different meshes create different grids

class Domain
{
public:
  //! Constructor
  /*!
    \param aptFileType (input) input file (mesh-data)
    \param ptOut (input) output file
    \param aptTimeFunc (input) time function data base
  */
  Domain(FileType * const aptFileType, WriteResults * ptOut, TimeFunc * aptTimeFunc);

  //!
  virtual ~Domain();

  //!
  /*!
    \param level index into grid hierarchy
  */
  void PrintGrid(const Integer level);

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
  /*!
    \param level index into hierarchy (multilevel methods)
  */
  void Update(const Integer level);

  //!
  Grid * GetGrid(){ return ptgrid_;}

protected:

private:
  
  Integer newlevel;

   //! initialize pde
   void InitPDEs();

   //! initialization of alg.sys.
  /*!
    \param level index into hierarchy (multilevel methods)
  */
   void InitAlgSys(const Integer level);

  //! update alg. sys. in case of new mesh
  /*!
    \param level index into hierarchy (multilevel methods)
  */
  void UpdateAlgSys(const Integer level);

  Integer numsubdomain_;  //!< number of subdomains
  Integer numsys_;        //!< number of systems (matrix dimension for algebraic system)
  Integer numpde_;        //!< number of PDEs
  Integer numgraph_;      //!< number of graphs needed (node-graphs, edge-graphs, etc.)
  Integer ** syscoupling_; //!< matrix, containing coupling information between the systems
  BasePDE * ptpde_[MAXNUMPDE]; //!< pointers to PDEs
  Grid * ptgrid_; //!< pointer to grid object
  BCs * ptBCs_;   //!< pointer to object storing boundary conditions
  AbstractAlgebraicSys * ptalgsys_; //!< pointer to algebraic system
  TimeFunc * ptTimeFunc_; //!< pointer to object handling time functions
  FileType *InFile_;      //!< pointer to object handling input file (mesh data)
  WriteResults * OutFile_; //!<  pointer to object handling output file 

};

}

#endif // FILE_DOMAIN
