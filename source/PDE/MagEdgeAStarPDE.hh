#ifndef FILE_MAGNETICEDGEASTARPDE
#define FILE_MAGNETICEDGEASTARPDE

#include <map>
#include "MagBasePDE.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class CurlCurlEdgeInt;
  class nLinCurlCurlEdgeInt;

  //! Class for 3D magnetics using edge elements
  class MagEdgeAStarPDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgeAStarPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeAStarPDE();

//    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    //stores the field intensity for energy-based hystersis/nonlinear models
    shared_ptr<CoefFunctionMulti> nonlinear_field_intensity_coef_;

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_POTENTIAL; }

  protected:
    
//    //! Initialize NonLinearities
//    virtual void InitNonLin();

//    //! read special boundary conditions (coils, magnets)
//    void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    void DefineStandardIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators(){};

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

//    //! define the SoltionStep-Driver
//    void DefineSolveStep();
    
    // =======================================================================
    //  Initialization
    // =======================================================================

    //! Define available primary result types
    void DefinePrimaryResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    //! Init the time stepping
    void InitTimeStepping();

    // =======================================================================
    //   COILS
    // =======================================================================
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );

    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation );

    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity );


    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    
//    //! Coefficient function, containing the overall reluctivity
//    shared_ptr<CoefFunctionMulti> reluc_;
//
//    //! Coefficient function, containing the overall conductivity
//    shared_ptr<CoefFunctionMulti> conduc_;
    
    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;

  private:
    //! Use gradient fields in shape functions (Edge elements of second kind)
    bool useGradFields_;
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeAStarPDE
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

