#ifndef FILE_ACOUSTIC2DPDE_2001
#define FILE_ACOUSTIC2DPDE_2001

#include "basepde.hh"
#include "abstractAS.hh"
 
namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step. In this class we work with complex WorkWithSystemMatrix. When we create object of this class, we determine what kind of matrix should be system matrix, for example, sparse matrix, band matrix or so on. We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class Acoustic2dPDE: virtual public BasePDE
{
public:

  //!
  Acoustic2dPDE(AbstractAlgebraicSys * aptalgsys, Grid<Point2D> * , Material * , TimeFunc * ,FileType * , WriteResults<Point2D> * );

  //!
  virtual ~Acoustic2dPDE();

  //!
  void SetAlgSys_id(const Integer as_sysid);

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse);

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer type);

    //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS();

  //! calculation derivates of solution 
  void CalculationDerivativesSol();

  //! create pointer to class for time error estimation
  virtual TimeErrorEstimator * CreatePtTimeError();  

  //! Calculation of energy norm
  Double CalcEnergyNorm();

  //!
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //!
  void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
   void WriteResultsInFile();

  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual const Vector<Double> & getS1() const { return sol_der1_;}

  //!
  virtual const Vector<Double> & getS2() const { return sol_der2_;}

  //!
  virtual const Vector<Double> & getS2old() const { return sol_der2_old_;}

  //!
  virtual Integer getSize() const { return size_;}

  //!
  Double getBeta() const { return beta_;}
  Double getGamma() const { return gamma_;}

private:
  //!
  Integer dofspernode_;

  //!
  Grid<Point2D> * ptgrid_;

  //!
  Integer AS_sysid_;

  //! Calculation parameters for Newmark method
  virtual void CalcParameters(const Double dt);

  //!
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;

  //! Integration parameters
  Double alpha_,gamma_, beta_;

  //! coefficient in equation
  Double coeff_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_, sol_der1_, sol_der2_, sol_old_, sol_der2_old_;  

  //! Last time on which we have calculated solution
  Double lasttimecalc_;

  //! Number of last timestep on which we have calculated our solution
  Integer laststepcalc_;

  //! size of solution and etc.
  Integer size_;
};

} // end of namespace
#endif
