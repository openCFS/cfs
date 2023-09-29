// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MagEdgeMixedSFGPDE.hh
 *       \brief    This class is used to solve for a magnetic field a solution that satisfies
 *                 Ampere's law (curlh = j), to obtain a proper source field hs. 
 *                 Since, the curlcurl equation alone gives a singular stiffness matrix a sort of
 *                 gauging is introduced that sets the divergence of h to zero (divh = 0)
 * 
 *                 This gives the following system of equations
 *                 \int curlh curlh' + \int gradu h' = \int j curlh' for all h' in H(curl)  (I)
 *                 \int h gradu'                     = 0             for all u' in H1       (II)
 * 
 *       \date     September 2023
 *       \author   ldomenig
 */
//================================================================================================

#ifndef FILE_MAGNETICEDGEMIXEDSFGPDE
#define FILE_MAGNETICEDGEMIXEDSFGPDE

#include <map>
#include "MagBasePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{
  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeMixedSFGPDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgeMixedSFGPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor
    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeMixedSFGPDE();

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
  };
#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeMixedAVPDE
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

