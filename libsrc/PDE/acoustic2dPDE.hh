#ifndef FILE_ACOUSTIC2DPDE_2001
#define FILE_ACOUSTIC2DPDE_2001

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for acoustic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  We set rules for assembling global system matrix according to weak form of PDE, define right hand side and set boundary conditions. Then we cause one of methods of LinSystem for solving linear system. On the last step we calculate first and second derivatives of the solution.
  */

class Acoustic2dPDE: public BasePDE
{

public:

  //!  Constructor. here we read integration parameters
  /*!
    \param  aptalgsys pointer to class Algebraic system
    \param aGrid pointer to grid
    \param aMatFile pointer to class Material. material data.
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
 Acoustic2dPDE(AbstractAlgebraicSys * aptalgsys, Grid * aGrid , Material * aMatFile, TimeFunc * aTimeFunc ,FileType * aInFile, 
	   WriteResults * aOutFile );

  //!  Deconstructor
  virtual ~Acoustic2dPDE() {;};

  //! specify type of solver for algebraic system. it is read from config-file
  /*!
    \param solvertype  Richardson or CG
    \param precondtype ID or MG
    \param eps relative accuracy in the precond. energy
    \param dampiter damping parameter for Jacobi, SSOR
    \param maxnumit max number of iterations
     \param numeqcoarse number of equation for coarsing
    \param coarsealpha coarsing parameter for AMG
  */
  virtual void SpecifySolver(Integer &solvertype, Integer &precondtype, Double &eps, Double &dampiter, Integer &maxnumit, 
		     Integer &numeqcoarse, Double &coarsealpha);

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
  virtual void SpecifyMatrices(Integer &matrixtype, Integer *matrixsystype, Integer &graphtype, Integer &numdofpernode, 
			       Integer &numdirichlets, Integer &numconstraints);

  //! set information for algebraic system about PDE. set matrix factors
  virtual void SetMatrixFactors();

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level, BCs * ptBCs=NULL);

    //! set boundary condition
  /*!
    \param ptBCs pointer to boundary condition
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atime time of calculation
  */
  virtual void SetBCs(BCs * ptBCs, const Integer level, const Integer update, const Double atime);

  //! compute rhs
  /*!
    \param atime time of calculation
    \param ptBCs pointer to class with BCs. can be omitted.
  */
  virtual void ComputeRHS(const Double atime, BCs * ptBCs=NULL);

  //! calculation derivates of solution  old stuff. not used.
  //  virtual void CalculationDerivativesSol(const Boolean Recalc);

  
  //!  solve one step for static problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param level level of grid
  */
  virtual void SolveStepStatic(BCs * ptBCs ,const Integer level)
  { 
    Error("Makes no sense for Acoustics to perform static step",__FILE__,__LINE__);
  }

  //! solve one step for transient problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param kstep number of calculating step
    \param steptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */
  virtual void SolveStepTrans(BCs * ptBCs ,const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat);

  //! save received solution as solution on the previous step
  virtual void SaveSolAsPrevStep();

  //! write results in file
   virtual void WriteResultsInFile();

  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const { return sol_;}

  //!  return pointer to vector with first derivative of solution
  virtual const Vector<Double> & getS1() const { return sol_der1_;}

  //!  return pointer to vector with first derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS1old() const { return sol_der1_old_;}

  //! return pointer to vector with second derivative of solution
  virtual const Vector<Double> & getS2() const { return sol_der2_;}

  //! return pointer to vector with second derivative of solution, calculated on previous step
  virtual const Vector<Double> & getS2old() const { return sol_der2_old_;}

  //! return size of solution
  virtual Integer getSize() const { return size_;}

  //! return parameter beta from Newmark method
  Double getBeta() const { return beta_;}

  //! return parameter gamma from Newmark method
  Double getGamma() const { return gamma_;}

 //! 
  virtual void CalcDerSol();

protected:
  //! some preliminary actions for calculation of RHS. they are the same for all steps. 
  virtual void preComputeRHS();

  //!
  Integer dofspernode_;

  //!
  Grid * ptgrid_;

  //! Calculation parameters for Newmark method
  virtual void CalcParameters(const Double dt);

  //!  calculation of coefficient.
  /*!
    \param coeffmass coefficient before mass matrix
    \param coeffstiff coefficient before stiffness matrix
  */
  void CalcCoeff(Vector<Double> & coeffmass, Vector<Double> & coeffstiff, Vector<Double> & coeffdamp);

  //! coefficients from Newmark method
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

  //! list of surfaces, on which we have force
  std::vector<std::string> rhs_surfaces_;

  //! list of bnds( for absorbing BCs)
  std::vector<std::string> bnd_absBCs_;

  //! function for RHS
  Integer arg_rhs_;
  pfn1var ptRHSFnc_;
  std::vector<Double> directionFnc_;

  //! Indicator: is there RHS function
  Boolean SetRHSFnc;

  //! Indicator: 
  Boolean without_absBCs_;

};

} // end of namespace
#endif
