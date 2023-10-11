// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MagEdgeMixedGaugedAPDE.hh
 *       \brief    In 3D the bilinear form (curlu,curlv) delivers singular FEM stiffness matrices
 *                 which is due to the not explicitly gauged formulation of the Maxwell equations
 *                 This PDE introduces the Coulomb gauging (diva = 0) by Lagrange multipliers,
 *                 that leads to the following system of equations

 *                 \int curla curla' + \int gradu a' = \int j a'     for all a' in H(curl)  (I)
 *                 \int a gradu'                     = 0             for all u' in H1       (II)
 * 
 *                 There exist a faster version that uses a regularization which produces similar
 *                 results as this formulation and is faster, because it only has one unknown.
 *                 For this reason this formulation should only be used as a reference solution 
 *                 to varify faster ones and not for real life problems.
 * 
 *       \date     October 2023
 *       \author   ldomenig
 */
//================================================================================================

#ifndef FILE_MAGNETICEDGEMIXEDGAUGEDAPDE
#define FILE_MAGNETICEDGEMIXEDGAUGEDAPDE

#include <map>
#include "MagBasePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{
  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeMixedGaugedAPDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgeMixedGaugedAPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor
    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeMixedGaugedAPDE();

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MAG_POTENTIAL; }

  protected:
    // =======================================================================
    //  DEFINE ALL BILINEAR- AND LINEARFORMS 
    // =======================================================================
    void DefineIntegrators();
    void DefineStandardIntegrators();
    LinearForm* GetCurrentDensityInt( Double factor, PtrCoefFct coef, std::string coilId = "" );
    //void DefineCoilIntegrators();

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

    //! Coefficient function, containing the overall permeability
    shared_ptr<CoefFunctionMulti> reluc_;

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

