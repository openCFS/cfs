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

  //!
    //! create pointer to class for time error estimation
  virtual SpaceErrorEstimator * CreatePtSpaceError();

private:
  //!
  Integer dofspernode_;

  //!
  void CalcCoeff(Vector<Double> & coeff);

  //!
  Grid * ptgrid_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_;  

  //! size of solution and etc.
  Integer size_;

};

} // end of namespace
#endif
