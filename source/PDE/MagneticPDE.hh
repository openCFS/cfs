// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_MAGNETIC_PDE_HH
#define FILE_CFS_MAGNETIC_PDE_HH

#include <map>
#include <set>
#include "MagBasePDE.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{
  //! Class for magnetic field calculation (eddy current case)
  
  //! This class is utilized to solve the magnetic problem (eddy current case)
  //! in 2D and 3D using nodal finite elements. The primary unknown is the 
  //! magnetic vector potential \[\vec{A}\]. 
  //! In case of 3D time dependent / harmonic simulations an additional scalar
  //! potential \[\Phi\] is utilized to ensure that the solution is embedded
  //! in the space HCurl.
  class MagneticPDE: public MagBasePDE {
    
  public:
    
    //! Constructor
    MagneticPDE( Grid *aGrid, PtrParamNode paramNode,
      PtrParamNode infoNode,
      shared_ptr<SimState> simState, Domain* domain );
    
    //!  Destructor
    virtual ~MagneticPDE();
    
    //! pass pointer to mechanicalPDE for later use in nonlinear material evaluation
    void SetMagnetoStrictCoupling(SinglePDE *mechanicPDE);
    
//    //! Get method for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);
    
    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_POTENTIAL; }

  protected:
    
//    //! Initialize NonLinearities
//    void InitNonLin();
    
//    //! read special boundary conditions (coils, magnets)
//    void ReadSpecialBCs();
    
    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
    
    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators();
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
//    //! define the SoltionStep-Driver
//    void DefineSolveStep();
    
    //! Read special results definition
    void ReadSpecialResults();
    
    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
    PtrParamNode infoNode );
    //! Flag for mixed formulation
    
//    //! In case of a transient / harmonic 3D simulation, we need
//    //! an additional scalar potential to ensure HCurl compatibility.
//    bool isMixed_;
    
    // =======================================================================
    //   COILS
    // =======================================================================
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );
    
    // =======================================================================
    //   HYSTERESIS
    // =======================================================================
    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation );

    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity );

    // =======================================================================
    //   NETWORK
    // =======================================================================
    void DefineLemIntegrators( void );

    //! Map containing the remanence (B excitation on RHS)
    //! needed for calculating H field
    std::map<RegionIdType,PtrCoefFct> bRHSRegions_;

//    //! Query parameter object for information on coils
//    void ReadCoils();
    
//    //! Initialize time stepping method
//    void InitTimeStepping();
    
    //! read in softening types
    void ReadSoftening();
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
    
    //! flag for magn_strict coupling
    bool isMagnetoStrictCoupled_;
        
    void CheckForConductivity();
    bool anyRegionHasConductivity_;

    SinglePDE *mechanicPDE_;
    
    void FinalizeAfterTimeStep();

    // FROM NACS
    //! In this set we store all coefficient functions which are related to
    //! the penalty term in case of 3D mixed nodal magnetics, which computes
    //! the grad-grad-term. We need it for the calculation of the eddy current
    //! density.
    std::set<shared_ptr<CoefFunctionFormBased> > mixedFormFunctor_;
    
  private:
    //! store velocity bilinear forms
    std::map<RegionIdType, BaseBDBInt*> velocityInts_;
    
  };
  
} // end of namespace
#endif

