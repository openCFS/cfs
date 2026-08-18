#ifndef FILE_MAGNETICEDGEPDE
#define FILE_MAGNETICEDGEPDE

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
  class MagEdgePDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgePDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgePDE();

//    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    //stores the field intensity for energy-based hystersis/nonlinear models
    std::map<RegionIdType, shared_ptr<CoefFunctionMulti> >  nlFieldIntensitym_;
    //shared_ptr<CoefFunctionMulti> nonlinear_field_intensity_coef_;

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
    void DefineNcIntegrators();

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
    
//    //! Query parameter object for information on coils
//    void ReadCoils();
    
//    //! Init the time stepping
//    void InitTimeStepping();
//
    // =======================================================================
    //   COILS
    // =======================================================================
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );
    
    // =======================================================================
    //   HYSTERESIS
    // =======================================================================
    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation );

    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity );
//
//    //@{ \name Attributes related to coils
//    //! Map CoilID to coil definition
//    std::map<Coil::IdType, shared_ptr<Coil> > coils_;
//
//    //! Map regionIds to coil definitions
//    typedef std::map<RegionIdType, shared_ptr<Coil> > CoilRegionMap;
//    CoilRegionMap coilRegions_;
//
//    //! Coefficients holding the current density for each coil
//    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;
//
//    //! Tells if there are coils excited by voltage
//    bool hasVoltCoils_;
//
//    //! Storage for CoefFunctions of external current density as source
//    std::map<shared_ptr<Coil::Part>, PtrCoefFct> coilPartsExtJ_;
//    //@}


    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    
//    //! Coefficient function, containing the overall reluctivity
//    shared_ptr<CoefFunctionMulti> reluc_;
//
//    //! Coefficient function, containing the overall conductivity
//    shared_ptr<CoefFunctionMulti> conduc_;

    //! Map containing the remanence (B excitation on RHS)
    //! needed for calculating H field
    std::map<RegionIdType,PtrCoefFct> bRHSRegions_;

    //! Set containing all regions with regularized conductivity
    
    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;
    
    //! Coefficient function for the multiharmonic material adaptions
    shared_ptr<CoefFunction> multiHarmCoef_;

  private:
    //! Use gradient fields in shape functions (Edge elements of second kind)
    bool useGradFields_;

    //! For the calculation of the source field Hs (magnetic scalar potential PDE),
    //! we need to consider every region as vacuum
    bool onlyVacuum_;

    //! store velocity bilinear forms
    std::map<RegionIdType, BaseBDBInt*> velocityInts_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgePDE
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

