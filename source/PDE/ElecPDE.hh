#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "SinglePDE.hh" 


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
  
  //! Class for electrostatic equation (no adaptivity)
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
    ElecPDE( Grid* grid, PtrParamNode paramNode );

    //! Destructor
    virtual ~ElecPDE(){};

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );
    
    void DefineSurfaceIntegrators(){};

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators();
    
    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Turn the piezo coupling on

    //! Triggers the correct assembly of the electrostatic block in a 
    //! piezo-coupled simulation, because the coupled electrostatic block
    //! is negative compared to the normal one
    void SetPiezoCoupling();

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return ELEC_POTENTIAL; }

    /** @see virtual SinglePDE::GetNativeDOF() */
    virtual UInt GetNativeDOF() const { return 1; }


  protected:

    //! SubType of electrostatic section
    std::string subType_;

    //! Return linear stiffness integrator for a given region
    BaseBDBInt * GetStiffIntegrator( BaseMaterial* actSDMat,
                                     SubTensorType tensorType,
                                     RegionIdType regionId );
    
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

    //! Read definitions for electric impedances
    void ReadImpedances();

    //! Define integrators for electric impedances
    void DefineImpedanceIntegrators();

    //! flag for piezo-coupling
    bool isPiezoCoupled_;



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
