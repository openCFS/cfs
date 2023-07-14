#ifndef FILE_FULLWAVEMAXWELLEPDE
#define FILE_FULLWAVEMAXWELLEPDE

#include <map>
#include "SinglePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{


  /** Class for 3D simulation with fullwave Maxwells equations in E-form */
  class FullWaveMaxwellEPDE : public SinglePDE
  {
  public:

    //! Constructor
    FullWaveMaxwellEPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~FullWaveMaxwellEPDE();


  protected:
    
    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs(){};

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! Define all RHS linearforms for load / excitation
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    // =======================================================================
    //  Initialization
    // =======================================================================

    //! Define available primary result types
    void DefinePrimaryResults();

    //! Define available postprocessing results
    void DefinePostProcResults();
    
    void FinalizePostProcResults();

    //! Init the time stepping
    void InitTimeStepping();

    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string& formulation,
                    PtrParamNode infoNode );

    //! Coefficient function, containing the overall permittivity
    shared_ptr<CoefFunctionMulti> eps_;

    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> reluc_;

    //! Coefficient function, containing the overall conductivity
    shared_ptr<CoefFunctionMulti> conduc_;

    //! Map containing the remanence (B excitation on RHS)
    //! needed for calculating H field
    std::map<RegionIdType,PtrCoefFct> bRHSRegions_;

    //! Set containing all regions with regularized conductivity

    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;

    /** String containing the information, which Darwin formulation we are using */
    std::string formulationType_;

  private:

    //! Define integrators for general coils/inductors
    void DefineGeneralCoilIntegrator();

    //! Define integrators for classical cylindrical coils
    void DefineCylindricalCoilIntegrator();

    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;
    
    //! flag for transient PML
    bool isTimeDomPML_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class FullWaveMaxwellEPDE
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve Encapsulate the defintiion of magnets within an own struct
  //! 

#endif

} // end of namespace
#endif

