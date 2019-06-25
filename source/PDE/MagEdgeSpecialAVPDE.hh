/*!
 *       \file     MagEdgeSpecialAVPDE.hh
 *       \brief    A-V formulation of the eddy current problem, whereas the
 *                 magnetic vector potential A is discretized with edge finite elements
 *                 and the electric scalar potential V with nodal finite elements.
 *                 This is the specialized version of the classical (full fledged) A-V
 *                 formulation, where we do not solve for A and V in the coupled way.
 *                 Instead we utilize that the grad of H1 functions is in Hcurl functionspace
 *                 and therefore the V-problem decouples from the A-V formulation and we can
 *                 use the solution of the electrokinetic problem (Laplace-type) for the electric
 *                 source potential and use it in the A-V formulation.
 *
 *
 *       \date     Apr. 22, 2019
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_MAGNETICEDGESPECIALAVPDE
#define FILE_MAGNETICEDGESPECIALAVPDE

#include <map>
#include "SinglePDE.hh" 
#include "MagEdgePDE.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{

  // Forward declaration of classes
  class CurlCurlEdgeInt;
  class nLinCurlCurlEdgeInt;

  class MagEdgeSpecialAVPDE : public MagEdgePDE
  {
  public:

    //! Constructor
    MagEdgeSpecialAVPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeSpecialAVPDE();


  protected:
    
    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    void DefineCoilIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

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

    bool useModifiedAVVoltageFormulation_;
    bool useModifiedAVCurrentFormulation_;

    //! Storage for CoefFunctions of external current density as source
//    std::map<shared_ptr<Coil::Part>, PtrCoefFct> coilPartsExtJ_;


    //! Storage for integral grad(V) \cdot grad(V) as source
    std::map<shared_ptr<Coil::Part>, Double> gradVsource_;


    //@}


    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    
  private:
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeSpecialAVPDE
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

