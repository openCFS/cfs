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
   virtual void SetMatrixFactors()
  {Error("Not implemented in AcousticPDE; uses TimeStepping-Class",__FILE__,__LINE__);}

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level);

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

  //! prepare for correct time stepping
  /*!
    \param dt time step
  */
  virtual void InitTimeStepping(const Double dt);

  //! write results in file
   virtual void WriteResultsInFile();

  //!  return pointer to vector with first derivative of solution
  virtual const Array<Double>& getS1() const { return TS_alg_->GetDeriv1();}

  //! return pointer to vector with second derivative of solution
  virtual const Array<Double>& getS2() const { return TS_alg_->GetDeriv2();}

  //! return size of solution
  virtual Integer getSize() const 
  { return size_;}

protected:

  Double lasttimecalc_;  //!< Last time on which we have calculated solution
  Integer laststepcalc_; //!< Number of last timestep on which we have calculated our solution

  Integer size_; //!< total number of unknowns (equations)

  Boolean with_absBCs_; //!< Indicator for absorbing boundary conditions 
  std::vector<std::string> bnd_absBCs_;   //!< list of bnds( for absorbing BCs)

  Boolean with_fracdamping_; //!< attenuation according to power law
  Integer frac_memory_;      //!< number of old time steps to be saved

  //General dimension of problem
  Integer dim_;
};

} // end of namespace
#endif
