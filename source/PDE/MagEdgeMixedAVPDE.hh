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
 *
 *       \date     Dec. 27, 2018
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_MAGNETICEDGEMIXEDAVPDE
#define FILE_MAGNETICEDGEMIXEDAVPDE

#include <map>
#include "SinglePDE.hh" 
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Utils/Coil.hh"

namespace CoupledField
{


  //! Class for 3D magnetics using edge elements and scalar nodal elements
  class MagEdgeMixedAVPDE : public SinglePDE
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
    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

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

    //! Tells if there are coils excited by voltage
    bool hasVoltCoils_;

    //! String defining if we have a classical cylindric coil or a generic inductor
    std::string coilFormulation_;

    //! boolean if we consider voltage coil in mixed form
    bool mixedCoil_;

    //! Storage for CoefFunctions of external current density as source
    std::map<shared_ptr<Coil::Part>, PtrCoefFct> coilPartsExtJ_;
    //@}


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
    void DefineGeneralCoilIntegrator();

    //! Define integrators for classical cylindrical coils
    void DefineCylindricalCoilIntegrator();

    //! Temporary fix, needed for postprocessing the gradient of
    //! the electric scalar potential V
    BaseBDBInt* gradInt_;
    shared_ptr<ResultInfo> elecIntens_;
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

