#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "SinglePDE.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

namespace CoupledField
{
  
  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class linElecInt;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  
  //! Class for electrostatic PDE
  class ElecPDE : public SinglePDE {
    
  public:
    
    //! Helper struct for defining a general impedance
    struct Impedance {
      
      //! Constructor
      Impedance() : 
      resistance(0.0),
      inductance(0.0), capacitance(0.0) {};
      
      //@{
      //! Node numbers the impedance is connected with
      shared_ptr<NodeList> node1, node2;
      //@}
      
      //@{
      //! Defining quantities of the impedance
      Double resistance, inductance, capacitance;
      //@}
    };
    
    //! Constructor
    /*!
     \param grid pointer to grid
     \param paramNode pointer to the corresponding parameter node
     */
    ElecPDE( Grid* grid, PtrParamNode paramNode,
    PtrParamNode infoNode,
    shared_ptr<SimState> simState, Domain* domain );
    
    //! Destructor
    virtual ~ElecPDE(){};
    
  protected:
    
    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();
    
    //! Initialize NonLinearities
    void InitNonLin();
    
    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );
    
    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();
    
    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators();
    
    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
    //! Define the SolveStep-Driver
    void DefineSolveStep();
    
    //! Init the time stepping:
    void InitTimeStepping();
    
    //! Nothing to do
    void SetTimeStep(const Double dt) {;};
    
    //! Read special boundary conditions
    void ReadSpecialBCs();
    
    void FinalizeAfterTimeStep();
    
    void InitHystCoefs();
    
    
    // ======================================================
    // COUPLING SECTION
    // ======================================================
    
  public:
    
    //! Turn the piezo coupling on
    
    //! Triggers the correct assembly of the electrostatic block in a 
    //! piezo-coupled simulation, because the coupled electrostatic block
    //! is negative compared to the normal one
    void SetPiezoCoupling();
    
  protected:
    
    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return ELEC_POTENTIAL; }
    
    //! SubType of electrostatic section
    std::string subType_;
    
    //! Return linear stiffness integrator for a given region
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, SubTensorType tensorType, RegionIdType regionId);
    
    //! Return linear stiffness integrator for a given region with the material tensor scaled by 'scalingFactor'
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, SubTensorType tensorType, RegionIdType regionId, PtrCoefFct scalingFactor);

    //! Return linear stiffness integrator for a given region with the material tensor scaled by 'scalingFactor'
    BaseBDBInt* GetStiffIntegratorInfMap(BaseMaterial* actSDMat, SubTensorType tensorType, RegionIdType regionId, PtrCoefFct scalingFactor);
    //! Return flux integrator used for Nitsche coupling
    template<typename DATA_TYPE>
    BiLinearForm* GetFluxIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec, Double factor,
    BiLinearForm::CouplingDirection cplDir, bool fluxOpA);
    
    //! Return penalty integrator used for Nitsche coupling
    template<typename DATA_TYPE>
    BiLinearForm* GetPenaltyIntegrator(PtrCoefFct scalCoefFunc, Double factor, BiLinearForm::CouplingDirection cplDir);
    
    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    //! Calculates the polarization vector
    void CalcPolarizationField( shared_ptr<BaseResult> vals );
    
    //! Calculate electric charges
    template <class TYPE>
    void CalcCharges( shared_ptr<BaseResult> vals );
    
    //! Contains the (Volume) subdomains next to the surface
    //! elements where the charges are computed
    std::map<shared_ptr<EntityList>,RegionIdType> chargeNeighborRegion_;
    
    //! Electric impedances
    StdVector<Impedance> impedances_;
    
    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
    PtrParamNode infoNode );
    
  private:
    
    bool regionApproxSet_;
    std::map<RegionIdType,shared_ptr<ElemList> > SDLists_;
    
    //! This map stores the polarization obtained by hysteresis coefFunctions
    shared_ptr<CoefFunctionMulti> hysteresisPolarization_;
    
    
    //! Read definitions for electric impedances
    void ReadImpedances();
    
    //! Define integrators for electric impedances
    void DefineImpedanceIntegrators();
    
    //! flag for piezo-coupling
    bool isPiezoCoupled_;
    
    //! Stores the dielectric permittivity for each region
    std::map<RegionIdType, PtrCoefFct > regionPermittivity_;
    
    //! Tensor type
    //    SubTensorType tensorType_;


  public:
    // =======================================================================
    // Nonlinear ElecPDE Section
    // =======================================================================
    /*
     * Common variables for elecPDE
     */
    //stores the material parameter
    shared_ptr<CoefFunctionMulti> epsilon_;

    //! Coefficient function for the multiharmonic material adaptions
    shared_ptr<CoefFunction> multiHarmCoef_;

    //stores the flux for hystersis and nonlinear models
    //std::map<RegionIdType, shared_ptr<CoefFunctionMulti> >  nlFluxCoefm_;

    //stores the flux for hystersis and nonlinear models
    shared_ptr<CoefFunctionMulti> nlFluxCoef_;
  };
  
#ifdef DOXYGEN_DETAILED_DOC
  
  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================
  
  //! \class ElecPDE
  //! 
  //! \purpose   
  //! This class is derived from class SinglePDE. It is used for solving
  //! electrostatic equation in 3D. 
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
  
} // end of namespace

#endif
