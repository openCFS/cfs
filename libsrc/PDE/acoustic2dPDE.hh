#ifndef FILE_ACOUSTIC2DPDE_2001
#define FILE_ACOUSTIC2DPDE_2001

#include "basepde.hh"
#include "outUnverg.hh"
 
namespace CoupledField
{

  //! Class for acoustic equation
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class Acoustic2dPDE: virtual public BasePDE
{
public:

  //  WriteResultsUnverg * ptOutput;

  //!
  Acoustic2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * , Material * , TimeFunc * ,FileType * , WriteResults * );

  //!
  virtual ~Acoustic2dPDE();

  //!
  void SpecifySolver(Integer &asolvertype, Integer &aprecondtype, Double &aeps,
Double &adampiter,  Integer &amaxnumit, Integer &numeqcoarse);

  //!
  void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, Integer &numdirichlets, Integer &numconstraints);

  //!
  void SetMatrixFactors();

  //!
  void SetupMatrices(const Integer level);

    //!
  void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //!
  void ComputeRHS(const Double atime, BCs * ptBCs=NULL);

  //! calculation derivates of solution 
  void CalculationDerivativesSol(const Boolean Recalc);
  void CalcDerSol();
  
  //! create pointer to class for time error estimation
  virtual TimeErrorEstimator * CreatePtTimeError();  

  //! create pointer to class for time error estimation
  virtual SpaceErrorEstimator * CreatePtSpaceError();

  //! Calculation of energy norm
  virtual Double CalcEnergyNorm();

  //!
  void SolveStepStatic(BCs * ptBCs ,const Integer level);

  //!
  void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, const Boolean updatesysmat);

  //!
   void SolveStepTransNewMesh(BCs * ptBCs, const Integer kstep, const Double asteptime, const Integer level);

  //!
  void RestoreSol();

  //!
  void SaveSolAsPrevStep();

  //!
   void WriteResultsInFile();

  //!
  virtual const Vector<Double> & getS() const { return sol_;}

  //!
  virtual const Vector<Double> & getS1() const { return sol_der1_;}

  //!
  virtual const Vector<Double> & getS1old() const { return sol_der1_old_;}

  //!
  virtual const Vector<Double> & getS2() const { return sol_der2_;}

  //!
  virtual const Vector<Double> & getS2old() const { return sol_der2_old_;}

  //!
  virtual Integer getSize() const { return size_;}

  //!
  Double getBeta() const { return beta_;}
  Double getGamma() const { return gamma_;}

  //! We use this function for time-error estimation
  void CalcThirdDerivateFromEquation(Vector<Double> & result);

private:
  //!
  Integer dofspernode_;

  //!
  Grid * ptgrid_;

  //! Calculation parameters for Newmark method
  virtual void CalcParameters(const Double dt);

  //!
  void CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff);

  //!
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;

  //! Integration parameters
  Double alpha_,gamma_, beta_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_, sol_der1_, sol_der2_, sol_old_, sol_der1_old_, sol_der2_old_;  

  Vector<Double> s_oldold_, s_d1_oldold_, s_d2_oldold_;

  //! Last time on which we have calculated solution
  Double lasttimecalc_;

  //! Number of last timestep on which we have calculated our solution
  Integer laststepcalc_;

  //! size of solution and etc.
  Integer size_;

};

// class PutElemMatInAlgSys
// {
// public:

//   PutElemMatInAlgSys(AbstractAlgebraicSys * aptalgsys, Grid * aptgrid, const Double acoeffm, const Double acoeffs, const Integer as_sysid, const Integer alevel)
//   { sysid_=as_sysid; ptalgsys_=aptalgsys; ptgrid_=aptgrid;
//     coeffm_=acoeffm; coeffs_=acoeffs; level_=alevel;
//     matrix_stiff_=2;  matrix_mass_=5; }

//   ~PutElemMatInAlgSys() {}

//   // method
//   void operator() (Elem t);

// private:

//      //!
//   Grid * ptgrid_;

//      //!
//   AbstractAlgebraicSys * ptalgsys_;

//      //!
//   Integer sysid_, level_;

//      //!
//   Double coeffm_, coeffs_;

//   //!
//   Integer matrix_stiff_, matrix_mass_;

// };

} // end of namespace
#endif
