// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MagEdgeMixedAVPDE.hh
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

#ifndef FILE_MAGNETICEDGEMIXEDAVPDE
#define FILE_MAGNETICEDGEMIXEDAVPDE

#include <map>
#include "MagBasePDE.hh"
#include "SinglePDE.hh" 
#include "MagEdgePDE.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{


  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeMixedAVPDE : public MagEdgePDE
  {
  public:

    //! Constructor
    MagEdgeMixedAVPDE( Grid * aptgrid, PtrParamNode paramNode,
                PtrParamNode infoNode,
                shared_ptr<SimState> simState, Domain* domain );

    //! Default Destructor

    //! The default destructor is responsible for freeing the Coil objects
    //! the ReadCoils() method brought into being.
    ~MagEdgeMixedAVPDE();

    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
    //! specifying coil results.
    //shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

  protected:
    
    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    void DefineAVIntegrators();
    
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
    //void ReadCoils();
    
    //! Init the time stepping
    void InitTimeStepping();

    //! \copydoc SinglePDE::CreateFeSpace
    std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string& formulation, 
                    PtrParamNode infoNode );
    
    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> reluc_;
    
    //! Coefficient function, containing the overall conductivity
    shared_ptr<CoefFunctionMulti> conduc_;

    //! Map containing the remanence (B excitation on RHS)
    //! needed for calculating H field
    std::map<RegionIdType,PtrCoefFct> bRHSRegions_;

    //! Set containing all regions with regularized conductivity
    
    //! This set contains all regions, which have no physical conductivity,
    //! but only a non-zero one due to regularization. This must be considered
    //! when calculating eddy currents, which would be "unphysical" in these
    //! regions.
    std::set<RegionIdType> regularizedRegions_;
    

  private:

    //! Define integrators for general coils/inductors
    //void DefineGeneralCoilIntegrator();

    //! Define integrators for classical cylindrical coils
    //void DefineCylindricalCoilIntegrator();

    //! For the calculation of the source field Hs (magnetic scalar potential PDE),
    //! we need to consider every region as vacuum
    //bool onlyVacuum_;

    //! store velocity bilinear forms
    std::map<RegionIdType, BaseBDBInt*> velocityInts_;
    std::map<RegionIdType, BaseBDBInt*> velocityInts_2;

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

