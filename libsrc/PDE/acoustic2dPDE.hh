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
  Acoustic2dPDE(Grid<Point2D> * , Material * , TimeFunc * ,FileType * , WriteResults<Point2D> * );

  //!
  virtual ~Acoustic2dPDE();

  //!
  void SetAlgSys_id(Integer as_sysid);

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
Double &adampiter,  Integer &amaxnumit);

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(AbstractAlgebraicSys * algsys, Integer type);

    //!
  void SetBCs(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS(AbstractAlgebraicSys *ptalgsys);

  //! calculation derivates of solution 
  void CalculationDerivativesSol();

  //!
  void SolveStepStatic(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs ,Integer level);

  //!
  void SolveStepTrans(AbstractAlgebraicSys *ptalgsys, BCs * ptBCs ,Integer kstep, Double steptime,
                      Integer level);

  //!
   void WriteResultsInFile();

  //!
  virtual Vector<Double> & getS() { return sol_;}

  //!
  virtual Vector<Double> & getS1() { return sol_der1_;}

  //!
  virtual Vector<Double> & getS2() { return sol_der2_;}

private:

  //!
  Integer dofspernode_;

  //!
  Grid<Point2D> * ptgrid_;

  //!
  Integer AS_sysid_;

  //!
  void CalcParamForNewmarkMethod(const Double dt);

  //!
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;

  //! Intergration parameters
  Double alpha_,gamma_, beta_;

  //! coefficient in equation
  Double coeff_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_, sol_der1_, sol_der2_, sol_old_;  

  //! Last time on which we have calculated solution
  Double lasttimecalc_;

  //! Number of last timestep on which we have calculated our solution
  Integer laststepcalc_;
};

inline Acoustic2dPDE::~Acoustic2dPDE()
{
 if (ptTimeFunc_) delete ptTimeFunc_;
}

} // end of namespace
#endif
