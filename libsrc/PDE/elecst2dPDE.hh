#ifndef FILE_ELECTST2DPDE_2002
#define FILE_ELECTST2DPDE_2002

#include "basepde.hh"
 
namespace CoupledField
{

  //! Class for electrostatic equation in 2D
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 2D. 
  */

class Elecst2dPDE: virtual public BasePDE
{
public:

  //!  Constructor. here we read integration parameters
  /*!
    \param  aptalgsys pointer to class Algebraic system
    \param aGrid pointer to Grid
    \param aMatFile pointer to class Material. material data.
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Elecst2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * aGrid, Material * aMatFile, TimeFunc * aTimeFunc, FileType * aInFile , WriteResults * aOutFile);

  //!  Deconstructor
  virtual ~Elecst2dPDE();

  //! specify type of solver for algebraic system. it is read from config-file
  /*!
    \param asolvertype  Richardson or CG
    \param aprecondtype ID or MG
    \param aeps relative accuracy in the precond. energy
    \param adampiter damping parameter for Jacobi, SSOR
    \param amaxnumit max number of iterations
     \param numeqcoarse number of equation for coarsing
    \param coarsealpha coarsing parameter for AMG
  */
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
		     Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse, Double &coarsealpha);

  //!  specify type of system matrix for AlgebraicSystem
  /*!
    \param matrixtype out: 0..NOCLASS, 1..RSPARSE, 2..CSPARSE, 3..RBLOCK, 4.. CBLOCK = 0,
 5..RFULL, 6..CFULL, 7..MIXED
    \param matrixsystype out:define need we memory for different types of element-matrix or not                     
    \param graphtype out: type of graph
    \param numdofpernode out: number of dof per node
    \param numdirichlets out:number of nodes for dirichlets conditions
    \param numconstraints out:number of nodes for constraints conditions
  */
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!  set information for algebraic system about PDE. set matrix factors.
  void SetMatrixFactors();

  //!  setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level=0, BCs * ptBCs=NULL);

    //! set boundary condition
  /*!
    \param ptBCs pointer to boundary condition
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atime time step of claculation
  */
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //! comput RHS 
  void ComputeRHS();

  //! calculation derivates of solution 
  void CalculationDerivativesSol();

  //! Calculation of energy norm
  Double CalcEnergyNorm();

  //! solve one step for static problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param level level of grid
  */
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //! write results in file
  void WriteResultsInFile();

  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const { return sol_;}

  //! return size of solution
  virtual Integer getSize() const { return size_;}
  
  //! test error of solution
  virtual Boolean TestError();

  //! refine mesh
  virtual void RefineMesh();

  //! write additional info (marked elements, relative error) to files with mesh
  virtual void PrintMeshesInfo(WriteResults * ptMeshes);

private:

  //!
  void CalcCoeff(Vector<Double> & coeff);

  //! calculation of gradient in element at integr-points
  void CalcGradSolElemIP(const Elem * ptElem, const Integer level, Vector<Double>*gradElIP);

  //!
  void CalcElectricField();

  //!
  Integer dofspernode_;

  //!
  Grid * ptgrid_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_,*solGrad_;  

  //! abs. value of electric field at center point of elem; for each elem
  Vector<Double> absValueElectricField_;

  //! size of solution and etc.
  Integer size_;

  //! indicator that for problem absolute value of electric field at center of each element is calculated
  Boolean AbsValueElectricField_;

  // ------------------- stuff for adaptivity ------------------------------

  //! initialization of pointer to SpaceErrorEstimator
  void ConstructorError();

   //! indicator that used in function SetupMatrices. TRUE, if we need mass matrix; otherwise - FALSE
  Boolean NeedMassMatrix_;

  //! indicator for WriteResults. TRUE, if we want to print error map.
  Boolean WriteErrorMap_;

  //! array, in which we store error map
  Vector<Double> errorMap_;
  Vector<Double> gradSPRElemL2norm_;
  Vector<Double> testMap_;
  Vector<Double> relativeErrorMap_;
  Vector<Double> markingElems_;
  Double normError_;
  Vector<Double> sqrDiffGradNodes_; //<! parameter storing difference of FEM-gradient and SPR-gradient at nodes
  Vector<Double> sqrGradFEMNodes_; //<! parameter storing FEM-gradient at nodes

  //! pointer to class of error estimators
  SpaceErrorEstimator<2> * ptError_;

  //! recovered solution by SPR
  Vector<Double> solSPR_;

  //! error tolerance
  Double errorTol_;

  // for testing: using own solver
  Boolean SolverCFS_; //<! parameter indicator: TRUE, if you want to use Solver CFS. reading from config-file
  Matrix<Double> sysmat_;
  Vector<Double> vecrhs_;
  
};

} // end of namespace
#endif
