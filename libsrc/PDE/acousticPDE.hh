#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for acoustic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  
  */

class AcousticPDE: public BasePDE
{

public:

  //!  Constructor. here we read integration parameters
  /*!
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  AcousticPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, 
	      WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~AcousticPDE() {;};

  //! define discrete PDE
  virtual void DiscreteParamsPDE();

  //! set information for algebraic system about PDE. set matrix factors
  virtual void SetMatrixFactors();

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level)
  { 
    Error("Not implemented in base class of Acoustics",__FILE__,__LINE__);
  }

  //! compute rhs
  /*!
    \param atime time of calculation
  */
  virtual void ComputeRHS(const Double atime);
  
  virtual void SolveStepStatic(const Integer level)
  { 
    Error("Makes no sense for Acoustics to perform static step",__FILE__,__LINE__);
  }

  //! solve one step for transient problem 
  /*!
    \param kstep number of calculating step
    \param steptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */
  virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
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

  //! Calculation parameters for Newmark method
  virtual void CalcParameters(const Double dt);

  //! coefficients from Newmark method
  Double a0_,a1_,a2_,a3_,a4_,a5_,a6_,a7_;

  //! Integration parameters
  Double alpha_,gamma_, beta_;

  //! store solution, 1st derivative , 2nd derivative solution
  Vector<Double> sol_, sol_der1_, sol_der2_, sol_old_, sol_der1_old_, sol_der2_old_;  

  Double lasttimecalc_;  //!< Last time on which we have calculated solution
  Integer laststepcalc_; //!< Number of last timestep on which we have calculated our solution

  Integer size_; //!< total number of unknowns (equations)

  Boolean with_absBCs_; //!< Indicator for absorbing boundary conditions 
  std::vector<std::string> bnd_absBCs_;   //!< list of bnds( for absorbing BCs)


};

} // end of namespace
#endif
