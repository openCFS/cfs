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

  //!
  Elecst2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * , Material * , TimeFunc * ,FileType * , WriteResults * );

  //!
  virtual ~Elecst2dPDE();

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse, Double &coarsealpha);

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer level=0);

    //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS();

  //! calculation derivates of solution 
  void CalculationDerivativesSol();

  //! Calculation of energy norm
  Double CalcEnergyNorm();

  //!
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //!
  void WriteResultsInFile();

  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual Integer getSize() const { return size_;}
  
  //! test error of solution
  virtual Boolean TestError();

  //! refine mesh
  virtual void RefineMesh();

  //! write additional info (marked elements, relative error) to files with mesh
  virtual void PrintMeshesInfo(WriteResults * ptMehes);

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
