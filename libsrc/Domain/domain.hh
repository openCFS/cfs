#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#define MAXNUMPDE 20
#define MAXNUMCOUPLEDPDE 20

#include <General/environment.hh>
#include <Domain/grid.hh>
#include <Domain/bcs.hh>

#include <PDE/basePDE.hh>


namespace CoupledField
{

class TimeFunc;
class WriteResults;
class BasePDE;
class BaseCoupledPDE;
class PDECoupling;
class  AbstractAlgebraicSys;

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

  //! get coupled pde
  BaseCoupledPDE * GetCoupledPDE() {return ptcoupledpde_;}

  //! get algebraic system
  // AbstractAlgebraicSys * GetAlgSys(){ return ptalgsys_;}

  //! get pointer to input-file
  FileType * GetInFile(){ return InFile_;}

  //! get pointer to output-file
  WriteResults * GetOutFile(){ return OutFile_;}

  //! get pointer to boundary condition
  BCs * GetBCs(){ return ptBCs_;}

  //! delete pointer to PDEs and create them new
  void ResetPDEs();


  //! update algebraic system and bcs after refinement the mesh
  /*!
    \param level index into hierarchy (multilevel methods)
  */
  void Update(const Integer level);

  //! update alg. sys. in case of new mesh
  /*!
    \param level index into hierarchy (multilevel methods)
  */
  void UpdateAlgSys(const Integer level);

  //!
  Grid * GetGrid(){ return ptgrid_;}

  //! get number of pdes
  Integer GetNumPDE() {return numpde_;}

protected:

private:
  
  

  //! initialize pde
  void InitPDEs();
  
  //! initialize coupled pde
  void InitCoupledPDE();

   //! initialization of alg.sys.
  /*!
    \param level index into hierarchy (multilevel methods)
  */
   void InitAlgSys(const Integer level);

  Integer numlevel_; //!< number of levels
  Integer numsubdomain_;  //!< number of subdomains
  Integer numsys_;        //!< number of systems (matrix dimension for algebraic system)
  Integer numpde_;        //!< number of PDEs
  Integer numcoupledpde_; //!< number of coupled PDEs
  Integer numgraph_;      //!< number of graphs needed (node-graphs, edge-graphs, etc.)
  Integer ** syscoupling_; //!< matrix, containing coupling information between the systems
  StdVector<BasePDE*> ptpde_;   //!< pointers to PDEs
  BaseCoupledPDE * ptcoupledpde_; //!< pointer to coupled PDEs
  StdVector<BasePDE*> orderedpdes_; //!<pointer to PDEs in right order for coupling
  StdVector<PDECoupling*> couplings_; //!<pointer to coupling objects
  Grid * ptgrid_; //!< pointer to grid object
  BCs * ptBCs_;   //!< pointer to object storing boundary conditions
  // AbstractAlgebraicSys * ptalgsys_; //!< pointer to algebraic system
  TimeFunc * ptTimeFunc_; //!< pointer to object handling time functions
  FileType *InFile_;      //!< pointer to object handling input file (mesh data)
  WriteResults * OutFile_; //!<  pointer to object handling output file 

};

}

#endif // FILE_DOMAIN
