#ifndef FILE_THERM2dPDE_2001
#define FILE_THERM2dPDE_2001

#include "basepde.hh"

namespace CoupledField
{

  /*! 
    This class is derived from class BasePDE. It is used for solving thermal equation. Generalized trapezoidal method. 
  */

class Therm2dPDE: virtual public BasePDE
{
public:

  //!
  Therm2dPDE(AbstractAlgebraicSys * ptalgsys, Grid * agrid, Material * aMatFile, TimeFunc * aptTimeFunc, FileType * aInFile, WriteResults * ); 

  //!
  virtual ~Therm2dPDE();

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse); 

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

  //! Calculation of parameters
  virtual void CalcParameters(const Double dt);

  //!
  void SolveStepStatic(BCs * ptBCs , const Integer level);

  //!
  void SolveStepTrans(BCs * ptBCs , const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual const Vector<Double> & getS1() const { return sol_der1_;}

  //!
  void WriteResultsInFile();
 
private:

  //!
  void Predictor();

  //!
  Integer dofspernode_;

  //!
  Integer doftype_;

  //!
  Grid * ptgrid_;

  //!
  Vector<Double> sol_, sol_pred_,sol_old_,sol_der1_,sol_der1_old_;

  //! parameters
  Double a0_, factor1_,factor2_,factor_pred_;  

  //! last time on which we have calculated solution
  Double lasttimecalc_;

  //! number of last timestep on which we have calculated our solution
  Double laststepcalc_;

  //! size of solution and etc.
  Integer size_;  

  //!
  Double gamma_parab_;  

  //! 
  /*!
   \param 1 .. effective stiffness formulation
   \param 2 .. effective mass formulation
  */

  Integer effsysmat_;

};

} // end of namespace
#endif
