#ifndef FILE_SMOOTHPDE
#define FILE_SMOOTHPDE

#include "SinglePDE.hh"

#include "CoupledPDE/pdecoupling.hh"
 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  class SmoothPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    SmoothPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

    //!  Deconstructor
    virtual ~SmoothPDE() {;};

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators();

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    //! Obtain information on desired output quantities from parameter file
    virtual void ReadStoreResults();

    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);

    //! initialize time stepping: nothing to do in smoother!
    virtual void InitTimeStepping(){;};

    //! set time step
    //! \param dt Current time step
    virtual void SetTimeStep(const Double dt){};
  

    //! calculate coupling terms
    virtual void CalcOutputCoupling();

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const UInt kstep = 0,
                                    const Double asteptime = 0.0,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);

    //! perform postprocessing on data
    void PostProcess() {;};
  
    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output);
  
  protected:

    UInt size_;        //!< total number of unknowns (equations)

  private:

    UInt GetNrBCDof (const std::string & dofStartString);

    UInt GetBCDof(const std::string dofString);

    std::string method_;

    Boolean firstTurn_;

    Vector<Double> factor_;

  };
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SmoothPDE
  //! 
  //! \purpose 
  //! This class implements a simple, quasi-mechanic grid smoother.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace
#endif
