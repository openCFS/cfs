#ifndef FILE_FLUIDMECHPERTURBEDPDE
#define FILE_FLUIDMECHPERTURBEDPDE

#include "SinglePDE.hh" 

namespace CoupledField
{

  // forward class declaration
  class BaseResult;
  class LinearFormContext;
  
  //! Class for linearized perturbed formulation of Navier-Stoke's equations.
  class FluidMechPerturbedPDE : public SinglePDE {

  public:

    //! Constructor
    /*!
      \param grid pointer to grid
      \param paramNode pointer to the corresponding parameter node
    */
    FluidMechPerturbedPDE( Grid* grid, PtrParamNode paramNode );

    //! Destructor
    virtual ~FluidMechPerturbedPDE(){};

    //! Initialize NonLinearities
    void InitNonLin();

    //! Define all (bilinearform) integrators needed for this PDE
    void DefineIntegrators( );

    void DefineSurfaceIntegrators(){};

    //! Define the SolveStep-Driver
    void DefineSolveStep();

    //! Init the time stepping: nothing to do
    void InitTimeStepping() {;};

    //! Nothing to do
    void SetTimeStep(const Double dt) {;};

    //! Read special boundary conditions
    void ReadSpecialBCs();

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> result );

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! Calculate coupling terms
    void CalcOutputCoupling();

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);
  
    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return ELEC_POTENTIAL; }

    /** @see virtual SinglePDE::GetNativeDOF() */
    virtual UInt GetNativeDOF() const { return 1; }
    
    //! Calculate field variables at arbitrary points
    void CalcField( SolutionType solType, StdVector<const Elem*>& elems,
                    StdVector<LocPoint>& points, SingleVector& values );



  protected:

    //! SubType of electrostatic section
    std::string subType_;

    // *****************
    //  POSTPROCESSING
    // *****************
    
    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    virtual LinearFormContext* CreateRhsLinearForm( SolutionType rhsType,
                                                    shared_ptr<CoefFunction > rhsCoef );

  private:

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);

    shared_ptr<BaseFieldFunctor> meanFlowFunctor_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class FluidMechPerturbedPDE
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
