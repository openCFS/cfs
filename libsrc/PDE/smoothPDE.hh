#ifndef FILE_SMOOTHPDE
#define FILE_SMOOTHPDE

#include "SinglePDE.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
    It is used for solving mechanic equation on one time step.  
  */

class SmoothPDE: public SinglePDE
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
  SmoothPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~SmoothPDE() {;};

  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);

  //! define the SoltionStep-Driver
  virtual void DefineSolveStep();

  //! Obtain information on desired output quantities from parameter file
  virtual void ReadStoreResults();

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! initialize time stepping: nothing to do in smoother!
  virtual void InitTimeStepping(){;};

  //! set time step
  //! \params dt Current time step
  virtual void SetTimeStep(const Double dt){};
  

  //! calculate coupling terms
  virtual void CalcOutputCoupling();

  //! write results in file
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time  
  virtual void WriteResultsInFile(Integer stepOffset = 0,
				  Double timeOffset = 0.0);

  //! perform postprocessing on data
  void PostProcess(Integer level) {;};
  
  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(SolutionType output);
  
protected:

  Integer size_;        //!< total number of unknowns (equations)

private:

  Integer GetNrBCDof (const std::string & dofStartString);

  Integer GetBCDof(const std::string dofString);

  std::string method_;

  Boolean firstTurn_;

  Vector<Double> factor_;

};

} // end of namespace
#endif
