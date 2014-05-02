#ifndef FILE_MAGNETICEDGEPDE
#define FILE_MAGNETICEDGEPDE

#include <map>
#include "SinglePDE.hh" 
#include "Driver/SolveSteps/SolveStepMagEdge.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class Coil;
  class CurlCurlEdgeInt;
  class nLinCurlCurlEdgeInt;
  class CoefFunctionMulti;

  //! Class for 3D magnetics using edge elements
  class MagEdgePDE : public SinglePDE
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
    
  protected:
    
    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators() {
      EXCEPTION("ncInterfaces are not implemented for MagEdgePDE");
    }

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();
    
    // =======================================================================
    //  Initialization
    // =======================================================================

    //! Define available primary result types
    void DefinePrimaryResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    //! Query parameter object for information on coils
    void ReadCoils();
    
    //! Init the time stepping
    void InitTimeStepping();

    // =======================================================================
    //   COILS
    // =======================================================================

    //@{ \name Attributes related to coils
    //! Map CoilID to coil definition
    std::map<Coil::IdType, shared_ptr<Coil> > coils_;

    //! Map regionIds to coil definitions 
    typedef std::map<RegionIdType, shared_ptr<Coil> > CoilRegionMap;
    CoilRegionMap coilRegions_;
    
    //! Coefficients holding the current density for each coil
    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;
    //@}


    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    
    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> reluc_;
    
    //! Coefficient function, containing the overall conductivity
    shared_ptr<CoefFunctionMulti> conduc_;
    
    //! Set containing all regions with regularized conductivity
    
    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;
    
  private:
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

