#ifndef FILE_ELECTST3DPDE_2002
#define FILE_ELECTST3DPDE_2002

#include "basepde.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class Elecst3dPDE: virtual public BasePDE
{
public:

  //! Constructor. here we read integration parameters
  /*!
    \param  aptalgsys pointer to class Algebraic system
    \param aMatFile pointer to class Material. material data.
    \param aGrid pointer to class Grid
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Elecst3dPDE(AbstractAlgebraicSys * aptalgsys, Grid * aGrid, Material * aMatFile, TimeFunc * aTimeFunc,FileType * aInFile , WriteResults * aOutFile);

  //! Deconstructor
  virtual ~Elecst3dPDE();

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
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //! set information for algebraic system about PDE. set matrix factors.
  void SetMatrixFactors();

  //! setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level=0, BCs * ptBCs=NULL);

    //!  set boundary condition
  /*!
    \param ptBCs pointer to boundary condition
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atime time step of claculation
  */
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //! compute rhs
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
  virtual const Vector<Vector<Double> > & getS() const { return sol_;}

  //! return size of solution
  virtual Integer getSize() const { return size_;}

  //! test error of solution
  virtual Boolean TestError();

  //! refine mesh
  virtual void RefineMesh();

  //! write additional info (marked elements, relative error) to files with mesh
  virtual void PrintMeshesInfo(WriteResults * ptMeshes);

private:
  //! initialization of pointer to SpaceErrorEstimator
  void ConstructorError();

  //!  calculation of error for each cell of mesh.
  void CalcErrorMap();

   //! calculation of error for the element of mesh
  void CalcErrorForElem(const Elem* elem, const Vector<Double>* gradSPR, Double & error, Double & normGradSPR);

  //! calculation of electric field
  void CalcElectricField();

  //!
  Integer dofspernode_;

  //! calculation of coefficient in equation
  void CalcCoeff(Vector<Double> & coeff);

  //! pointer to the class Grid
  Grid * ptgrid_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Vector<Double> > sol_;  

  //! size of solution and etc.
  Integer size_;

  //! indicator for WriteResults. TRUE, if we want to print error map.
  Boolean WriteErrorMap_;

  //! indicator that we have to calculate electric field
  Boolean calcElecField_;
  //! arrays we store values of electric field at the center of element
  Vector<Double> elecFieldAtCenterElem_[3];
 
  //! array, in which we store error map
  Vector<Double> errorMap_;
  //! array, in which we store l2 norm of gard SPR
  Vector<Double> gradSPRElemL2norm_;
  //! array, in which we store marked elements
  Vector<Double> markingElems_;
  //! array, in which we store value of relative error
  Vector<Double> relativeErrorMap_;
  //! error tolerance from conf-file
  Double errorTol_;
  //! l2 norm of SPR gradient
  Double normError_;
  
  //! pointer to class of error estimators
  SpaceErrorEstimator<3> * ptError_;
  
};

} // end of namespace
#endif
