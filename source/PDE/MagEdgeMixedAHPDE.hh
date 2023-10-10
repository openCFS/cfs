// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MagEdgeMixedAHPDE.hh
 *       \brief    A-V formulation of the eddy current problem, whereas the
 *                 magnetic vector potential A is discretized with edge finite elements
 *                 and the electric scalar potential V with nodal finite elements.
 *                 This is the classical (full fledged) A-V formulation, where we
 *                 solve for A and V in a coupled way (mixed problem). There is actually
 *                 an alternative to solving the full coupled problem by utilizing the
 *                 properties of edge and nodal finite element function spaces, which is
 *                 realized in MagEdgeSpecialAVPDE (child-class of MagEdgePDE).
 *
 *       \date     Dec. 27, 2018
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_MAGNETICEDGEMIXEDAHPDE
#define FILE_MAGNETICEDGEMIXEDAHPDE

#include <map>
#include "MagBasePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{


  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeMixedAHPDE : public MagBasePDE
  {
  public:

    //! Constructor
    MagEdgeMixedAHPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeMixedAHPDE();

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

    
    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> perm_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MagEdgeMixedAHPDE
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

