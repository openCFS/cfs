#ifndef FILE_MAGNETICEDGEHPDE
#define FILE_MAGNETICEDGEHPDE

#include <map>
#include "MagBasePDE.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class CurlCurlEdgeInt;
  class nLinCurlCurlEdgeInt;

  //! Class for 3D magnetics H-formulation using edge elements
  class MagEdgeHPDE : public MagBasePDE 
  {
  public:

    //! Constructor
    MagEdgeHPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeHPDE();

//    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_FIELD_INTENSITY; }

  protected:
    

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
    void DefineStandardIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators(){};

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){}; // not implemented

    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ){}; // not implemented
    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ){}; // not implemented

    //! Define available primary result types
    void DefinePrimaryResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    // =======================================================================
    //   COILS
    // =======================================================================
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );


    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    

    //! Set containing all regions with regularized conductivity
    
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

  //! \class MagEdgeHPDE
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