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
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps, Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse, Double &coarsealpha); 

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer level);

  //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS(const Double atime, BCs * ptBCs); 

  //!
  void SolveStepStatic(BCs * ptBCs , const Integer level);

  //!
  void SolveStepTrans(BCs * ptBCs , const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
  void SaveSolAsPrevStep();
  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual const Vector<Double> & getS1() const { return sol_der1_;}

  //!
  void WriteResultsInFile();
 
private:
  //! calculation derivates of solution
  void CalcDerSol();

  //!
  void Predictor();

  //!
  Integer dofspernode_;

  void CalcParameters(const Double dt);

  //!
  Grid * ptgrid_;

  //!
  Vector<Double> sol_, sol_der1_,sol_pred_,sol_old_,sol_der1_old_;

  //! parameters
  Double a0_, factor1_,factor2_,factor_pred_;  

  //! last time on which we have calculated solution
  Double lasttimecalc_;

  //! number of last timestep on which we have calculated our solution
  Integer laststepcalc_;

  //! size of solution and etc.
  Integer size_;  

  //!
  Double gamma_parab_;  

  //! Parameter for type of formulation
  /*!
   \param 0 .. effective stiffness formulation
   \param 1 .. effective mass formulation
  */
  Integer formulation_;

};

} // end of namespace
#endif
