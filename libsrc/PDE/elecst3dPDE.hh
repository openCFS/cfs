#ifndef FILE_ELECTST3DPDE_2001
#define FILE_ELECTST3DPDE_2001

#include "basepde.hh"
#include "abstractAS.hh"
 
namespace CoupledField
{

  //! Class for electrostatic equation in 3D
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class Elecst3dPDE: virtual public BasePDE
{
public:

  //!
  Elecst3dPDE(AbstractAlgebraicSys * aptalgsys, Grid * , Material * , TimeFunc * ,FileType * , WriteResults * );

  //!
  virtual ~Elecst3dPDE();

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

private:
  //!
  Integer dofspernode_;

  //!
  void CalcCoeff(Double & coeff, const Integer numsubdom);

  //!
  Grid * ptgrid_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_;  

  //! Last time on which we have calculated solution
  Double lasttimecalc_;

  //! Number of last timestep on which we have calculated our solution
  Integer laststepcalc_;

  //! size of solution and etc.
  Integer size_;

  //! type of dof
  Integer doftype_;
};

} // end of namespace
#endif
