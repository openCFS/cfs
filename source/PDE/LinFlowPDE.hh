// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LinFlowPDE.hh
 *       \brief    Linearized flow equations (incompressible and compressible version);
 *                 A background flow can also be considered
 *       \date     Oct 21, 2018
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#ifndef FILE_LINFLUIDMECHPDE
#define FILE_LINFLUIDMECHPDE

#include "SinglePDE.hh" 

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class LinearFormContext;
  class CoefFunctionMulti;
  
  //! Class for linearized flow equations (incompressible and compressible version).
  class LinFlowPDE : public SinglePDE {

  public:

    //! Constructor
    /*!
      \param grid pointer to grid
      \param paramNode pointer to the corresponding parameter node
    */
    LinFlowPDE( Grid* grid, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~LinFlowPDE(){};

    //! set coupling to Heat PDE and whteher to use a symmetric form
    void SetHeatPDECouplingFlags(bool useSymmetricForm);

    //! Return the integrator sign
    double GetBalanceOfMomentumSign() const;

  protected:

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this PDE
    void DefineIntegrators( );

    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! Defines integrators for Nitsche coupling of an unknown on one specific
    //! interface.
    template<UInt DIM, UInt D_DOF, UInt P_DOF>
    void DefineNitscheCoupling( NcInterfaceInfo &iface,
                                shared_ptr<CoefFunctionMulti> additionalCoef = NULL);

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators();

    //! Define all RHS linearforms for load / excitation / tractions
    void DefineRhsLoadIntegrators(PtrParamNode input = PtrParamNode());

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping();

    //! Read special boundary conditions (nothing to do)
    void ReadSpecialBCs() {;};

    //! Calculate field variables at arbitrary points
    void CalcField( SolutionType solType, StdVector<const Elem*>& elems,
                    StdVector<LocPoint>& points, SingleVector& values );

//    //! Finalize before a time step is done
//    void FinilizeBeforTimeStep();

    //! SubType (plane, axi, 3D)
    std::string subType_;

    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! \copydoc SinglePDE::FinalizePostProcResults
    void FinalizePostProcResults();
     
    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );


  private:

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);
    
    BaseBDBInt *
    GetStiffIntegrator( BaseMaterial* actSDMat,
                        RegionIdType regionId,
                        bool isComplex );

    //! Add surface integrators stemming from partially integrating the stress tensor
    //! We usually do not need these terms when prescribing the velocity directly, but we need them for e.g. LM-based constraints or the doNothingBC
    void AddSurfaceIntegratorFromPartialIntegration(UInt i, StdVector<std::string> volumeRegions,
                                                    StdVector<shared_ptr<EntityList> > ent, std::string integrator);

    //! Coefficient function for the flow field
    
    //! This coefficient function describes the flow field. As this 
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanVelocityCoef_;
    shared_ptr<CoefFunction> meanFlowCoefScattered_;
    shared_ptr<BaseFeFunction> meanFlowFeFct_;
    shared_ptr<FeSpace> meanFlowFeSpace_;

    //! Slip condition
    bool isSlip_;
    shared_ptr<CoefFunctionMulti> meanFreePathCoef_;

    //! Use lagrange multiplier to enforce constraints
    bool useLagrangeMultVec_;

    //! Use lagrange multiplier 1 to enforce scalar constraints
    bool useLagrangeMultScal_;

    //! Surface regions on which pressure surface integral has to be defined
    StdVector<RegionIdType> presSurfaces_;

    //! compressible formulation
    bool isCompressible_;

    //! order of FE basis function for pressure
    std::string presPolyId_;

    //! order of FE basis function for velocity
    std::string velPolyId_;

    //! order of FE basis function for lagrange multiplier
    std::string lagrangeMultPolyId_;
    
    //! order of integration for pressure
    std::string presIntegId_;

    //! order of integration for velocity
    std::string velIntegId_;

    //! order of integration lagrange multiplier
    std::string lagrangeMultIntegId_;

    //! considers a varying background flow
    bool enableC2_;
    bool enableC3_;

    //! actiavtes convective term in the conservation of mass caused by the grid velocity (ALE-stuff)
    bool enableGridVelC1_;

    //! actiavtes convective term in the conservation of momentum caused by the grid velocity (ALE-stuff)
    bool enableGridVelC2_;

    //! amplitude factor for first convective term in case of varying background flow
    Double factorC1_;

    //! true, if coupled to Heat PDE
    bool isHeatPDECoupled_;
    //! Whether to use a symmetric formulation or not in coupling (to HeatPDE)
    bool isCouplingFormulationSymmetric_;

    //! Sign of the integrators (used to controll the sign of the balance of momentum terms)
    double balanceOfMomentumSign_;

  private:

    //! Coefficient function for the moving mesh

    //! This coefficient function describes the moving mesh for the ALE formulation. As this
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> gridVelCoef_;
    shared_ptr<CoefFunction> gridVelCoef_Scattered_;
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class LinFlowPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! a linearized perturbed formulation of Navier-Stoke's equations in 3D.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace CoupledField

#endif // FILE_FLUIDMECHPERTURBEDPDE
