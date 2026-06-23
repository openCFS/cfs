// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================

#ifndef FILE_MAGNETICEDGEREGULSFGPDE
#define FILE_MAGNETICEDGEREGULSFGPDE

#include <map>
#include "MagBasePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{
  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeRegulSFGPDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgeRegulSFGPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor
    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeRegulSFGPDE();

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_FIELD_INTENSITY; }

  protected:
    // =======================================================================
    //  DEFINE ALL BILINEAR- AND LINEARFORMS 
    // =======================================================================
    void DefineIntegrators();
    void DefineStandardIntegrators();
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );
    void DefineCoilIntegrators();

    // =======================================================================
    //  NOT IMPLEMENTED BUT NECESSARY METHODS
    // =======================================================================
    void DefineNcIntegrators(){}; // not implemented
    void DefineSurfaceIntegrators( ){}; // not implemented
    LinearForm* GetRHSMagnetizationInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ){return NULL;}; // not implemented
    BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ){return NULL;}; // not implemented
    void InitMagnetization(){}; // not implemented
    virtual void InitNonLin(){}; // not implemented
    
    // =======================================================================
    // INITIALIZATION
    // =======================================================================
    //! Define available primary result types
    void DefinePrimaryResults();
    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
    //! Define available postprocessing results
    void DefinePostProcResults();
    //! Init the time stepping
    void InitTimeStepping();
    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );

  private:
    Double mu0_;
    Double rho0_;
  };
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeMixedAVPDE
  //! 
  //! \purpose 
  //! This class is used to solve for a two-dimensional magnetic field 
  //! a solution that satisfies Ampere's law (curlh = j).
  //! Since the curlcurl equation alone gives a singular stiffness matrix,
  //! gauging is introduced via Lagrange multipliers to ensure a
  //! Coulomb-like gauge condition, i.e., div h_{sc} = 0.

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

