#ifndef FILE_ELEC2DPDE_2001
#define FILE_ELEC2DPDE_2001

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for electrostatic 2D (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class Elec2dPDE: public BasePDE
{
public:
    
  //!  Constructor. here we read integration parameters
  /*!
    \param  aptalgsys pointer to class Algebraic system
    \param aGrid pointer to grid
    \param aMatFile pointer to class Material. material data.
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Elec2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * aGrid , Material * aMatFile, TimeFunc * aTimeFunc ,
	    FileType * aInFile, WriteResults * aOutFile );

  //!  Deconstructor
  virtual ~Elec2dPDE();

  //! specify type of solver for algebraic system. it is read from config-file
  /*!
    \param solvertype  Richardson, CG, etc.
    \param precondtype ID, MG, etc.
    \param eps relative accuracy in the precond. energy
    \param dampiter damping parameter for iterative solver
    \param maxnumit max number of iterations
     \param numeqcoarse number of equation for coarsing
    \param coarsealpha coarsing parameter for AMG
  */
  virtual void SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, 
			     Integer &maxnumit, Integer &numeqcoarse, Double &coarsealpha);

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param matrixtype out: 0..NOCLASS, 1..RSPARSE, 2..CSPARSE, 3..RBLOCK, 4.. CBLOCK = 0,
    5..RFULL, 6..CFULL, 7..MIXED
    \param matrixsystype out:define need we memory for different types of element-matrix or not
    \param graphtype out: type of graph
    \param numdofpernode out: number of dof per node
    \param numdirichlets out:number of nodes for dirichlets conditions
    \param numconstraints out:number of nodes for constraints conditions
  */
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, 
			       Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //! set information for algebraic system about PDE. set matrix factors
  virtual void SetMatrixFactors();

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level, BCs * ptBCs=NULL);

  //! set boundary condition
  /*!
    \param ptBCs pointer to boundary condition
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atime time of calculation
  */
  virtual void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);


  //! solve one step for static problem 
  /*! 
    \param ptBCs pointer to class with data about boundary condition
    \param level level of grid
  */
  virtual void SolveStepStatic(BCs * ptBCs ,const Integer level);


  virtual void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat)
  { 
    Error("Makes no sense for Electrostatics to perform transient step",__FILE__,__LINE__);
  }

  //! write results in file
  virtual void WriteResultsInFile();


  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const { return sol_;};

protected:

  //! number of unknowns per node
  Integer dofspernode_;

  //! pointer to Grid
  Grid * ptgrid_;

  //! stores electric potential
  Vector<Double> sol_;

  //!  calculation of coefficient.
  /*!
    \param coeff electric permittivity
  */
  //!
  void CalcCoeff(Vector<Double> & coeff);

  //! size of solution and etc.
  Integer size_;

  //! function for RHS
  Integer arg_rhs_;
  pfn1var ptRHSFnc_;
  std::vector<Double> directionFnc_;

  //! Indicator: is there RHS function
  Boolean SetRHSFnc;

};

} // end of namespace
#endif
