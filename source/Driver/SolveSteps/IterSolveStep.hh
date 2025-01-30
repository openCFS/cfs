#ifndef ITERSOLVESTEP_HH
#define ITERSOLVESTEP_HH

#include "BaseSolveStep.hh"

#include "PDE/BasePDE.hh"
#include "MatVec/SingleVector.hh"
#include "FeBasis/FeFunctions.hh"
#include "Domain/CoefFunction/CoefFunctionAccumulator.hh"

namespace CoupledField
{

  // forward class declarations
  class IterCoupledPDE;

  
  // ======================================================================
  //  Classes for Convergence Criterions
  // ======================================================================
  
  //! Base class for handling convergence criterions
  
  //! This class handles a convergence criterion defined for a single
  //! quantity, maybe on several entitylists / regions. This defined
  //! just an interface for a state-dependent approach, i.e. 
  class ConvCriterion {
  public:
    
    //! Enumeration for types of norms
    //! L2ABS = absolute L2-norm
    //! L2REL = relative L2 norm: (|val| - |oldval|) / |val|
    typedef enum {NO_NORM, L2ABS, L2REL} NormType;
    static Enum<NormType> NormTypeEnum;
    

    //! Constructor
    
    //! Constructor for convergence criterion
    //! \param type Type of norm to be applied (none, absolute, relative)
    //! \param value Final value of norm to be reached
    ConvCriterion( NormType type, Double value );
    
    //! Destructor
    virtual ~ConvCriterion() {};
    
    //! Reset values
    virtual void ResetValues() = 0;
    
    //! Initiate phase of convergence sampling
    virtual void StartSampling() = 0;
    
    //! Stop phase of convergence sampling
    virtual void StopSampling() = 0;

    //! Return current norm value
    virtual Double GetNorm() = 0;
    
    //! Return final norm value to be reached
    Double GetFinalNorm() {
      return finalNorm_;
    }
    
    //! Return flag, if criterion is already converged
    virtual bool Converged() = 0; 
    
    //! Get support (= region) of convergence criterion
    virtual StdVector<shared_ptr<EntityList> > GetSupport() = 0;
    
  protected:
    
    //! Calculate norm according to norm type
    Double CalcNorm( Double newVal, Double oldVal );
    
    //! Norm value to be reached
    Double finalNorm_;
    
    //! NormType
    NormType normType_; 
  };
  
  //! Calculate convergence by utilizing accumulated values
  class ConvCriterionAccu : public ConvCriterion {
  public:
    //! Constructor
    
    //! Constructor for convergence criterion
    //! \param type Type of norm to be applied (none, absolute, relative)
    //! \param value Final value of norm to be reached
    ConvCriterionAccu(NormType type, Double value, bool overrideNumInt = false );
    
    //! Destructor
    ~ConvCriterionAccu();
    
    //! Add a coefficient function accumulator for a on a given list
    void AddCoefFct( shared_ptr<EntityList> list,
                     shared_ptr<CoefFunctionAccumulator> coefFct );
    
    //! \copydoc ConvCriterion::ResetValues
    virtual void ResetValues() override;
    
    //! \copydoc ConvCriterion::StartSampling
    virtual void StartSampling() override;
    
    //! \copydoc ConvCriterion::StopSampling
    virtual void StopSampling() override;
    
    //! \copydoc ConvCriterion::GetNorm
    Double GetNorm() override;
    
    //! \copydoc ConvCriterion::Converged
    bool Converged() override;
    
    //! \copydoc ConvCriterion::GetSupport
    StdVector<shared_ptr<EntityList> > GetSupport() override;
    
  private:
    //! Map for each entitylist / region (key) an coeffunction accumulator 
    std::map<shared_ptr<EntityList>, 
            shared_ptr<CoefFunctionAccumulator> > coefs_;
    
    //! Actual norm value for quantity
    Double actNorm_;
    
    //! Old norm value for quantity
    Double oldNorm_;

    //! Bool to override numerical integration if necessary (e.g. when lpm is not fully defined)
    bool overrideNumInt_;
  };
  
  //! Special convergence criterion for displacement based values
  
  //! This convergence criterion internally holds a FeFcuntion for the
  //! displacement, which will be passed to the Grid as coordinate
  //! offset in the end.
  class ConvCriterionDisplacement : public ConvCriterion {
  public:
    //! Constructor
    
    //! Constructor for convergence criterion
    //! \param type Type of norm to be applied (none, absolute, relative)
    //! \param value Final value of norm to be reached
    ConvCriterionDisplacement(NormType type, Double value );
    
    //! Destructor
    virtual ~ConvCriterionDisplacement();
    
    //! Set FeFunction for displacement
    void SetDispFct( shared_ptr<FeFunction<Double> > disp);

    //! Set FeFunction for complex displacement
    void SetDispFctComplex( shared_ptr<FeFunction<Complex> > dispComplex);
    
    //! Set FeFunction for velocity (needed for geometry update)
    void SetVelFct( shared_ptr<FeFunction<Double> > vel);

    //! Set FeFunction for displacement (needed for geometry update)
    void SetAccFct( shared_ptr<FeFunction<Double> > acc);

    //! Add a new region to 
    void AddRegion(RegionIdType region );
    
    void SetNormFlag( bool justNorm);
    
    //! \copydoc ConvCriterion::ResetValues
    virtual void ResetValues() override;
    
    //! \copydoc ConvCriterion::StartSampling
    virtual void StartSampling() override;
    
    //! \copydoc ConvCriterion::StopSampling
    virtual void StopSampling() override;
    
    //! \copydoc ConvCriterion::GetNorm
    Double GetNorm() override;
    
    //! \copydoc ConvCriterion::Converged
    bool Converged() override;
    
    //! \copydoc ConvCriterion::GetSupport
    StdVector<shared_ptr<EntityList> > GetSupport() override;

    //! Get if PDE is complex valued
    bool GetIsComplex() { return isComplex_; };

    //! Set if PDE is complex valued
    void SetIsComplex( bool isComplex ) { isComplex_ = isComplex; };

    //! Set function to set either activate or deactivate the usage of a reference node for the phase correction in harmonic predeformed simulations
    void SetEnableRefNode( bool refNodeEnabled ) { refNodeEnabled_ = refNodeEnabled;};

    //! Set function to set the reference node name
    void SetRefNodeName( std::string refNodeName ) { refNodeName_ = refNodeName;};

    //! Set function to set the reference node DOF
    void SetRefNodeDOF( UInt refNodeDOF ) { refNodeDOF_ = refNodeDOF;};
    
    //! Set function to set the phase offset
    void SetPhaseOffset( double phaseOffset ) { phaseOffset_ = phaseOffset;};
    
  protected:
    
    //! Pointer to displacement FeFunction
    shared_ptr<FeFunction<Double> > disp_mech_;

    //! Pointer to displacement FeFunction
    //shared_ptr<FeFunction<Double> > disp_smooth_;

    //! Pointer to displacement FeFunction
    shared_ptr<FeFunction<Double> > disp_;

    //! Pointer to complex valued displacement FeFunction
    shared_ptr<FeFunction<Complex> > dispComplex_;

    //! Pointer to velocity FeFunction
    shared_ptr<FeFunction<Double> > vel_;

    //! Pointer to acceleration FeFunction
    shared_ptr<FeFunction<Double> > acc_;
   
    //! Set containing updated regions
    std::set<RegionIdType> updatedRegions_;
    
    //! Current norm of updated displacement values
    Double actNorm_;
    
    //! Old norm of updated displacement values
    Double oldNorm_;
    
    //! if just norm of mechanical displacement is of interest, but no updated geometry is needed
    bool justNorm_;

    //! bool to check if the feFunction is complex valued
    bool isComplex_ = false;

    //! enables using a reference node to which the phase is corrected to
    bool refNodeEnabled_ = false;

    //! name of the reference node used to compute the phase correction
    std::string refNodeName_;

    //! DOF of the reference node to which we map the phase to
    UInt refNodeDOF_;

    //! phase offset that is used in addition to the node based phase correction
    double phaseOffset_;

    //! phase correction including phase of DOF of node and user specified correction
    double phaseCorrection_;

    //! phase correction multiplier
    Complex phaseCorrMult_ = 1.0;
  };
  
  // ======================================================================
  
  //! Derived class for step-wise solving of iterative coupled StdPDEs
  class IterSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    IterSolveStep(IterCoupledPDE& apde, PtrParamNode paramNode, PtrParamNode infoNode);

    //! Destructor
    virtual ~IterSolveStep();

    //! Initialize struct
    void Init();
    
    //! Obtain coupling quantity

    //! This method returns a given coefficient function
    //! from a contained SinglePDE. Internally, it creates an additional
    //! surrounding struct to calculate some norm for determining an
    //! stopping criterion when evaluating the CoefFunction.
    PtrCoefFct GetCouplingCoefFct( SolutionType type,
                                   shared_ptr<EntityList>  list,
                                   const std::string& pdeName,
                                   bool& updatedGeo );

    //! Manually trigger the finalize procedure if it has not been done yet
    void TriggerFinalize();
    
    //! Function similar to GetCouplingCoefFct but only returns updatedGeo
    void GetUpdateGeoForPDE( SolutionType type,
                             shared_ptr<EntityList>  list,
                             const std::string& pdeName,
                             bool& updatedGeo );

    //----------------------- STATIC---------------------------------------

    /** routine for initializations before execution the SolveStep-method.
     * Is abstract in the BaseClass - it is poor design to enforce abstract methods,
     * when they are not really necessary */
    virtual void PreStepStatic() override {}

    //! forward to all contained PDEs (see BaseSolveStep::RefreshPrescribed())
    virtual void RefreshPrescribed() override;
 
    /** base method for solving one static step */
    virtual void SolveStepStatic() override;

    //! routine for actions after the SolveStep-method
    virtual void PostStepStatic() override {}

    //----------------------- TRANSIENT---------------------------------------

    //! Initialize additional data-structures as needed for the glm
    virtual void InitTimeStepping() override;
    
    //! routine for initializations before execution the SolveStep-method
    virtual void PreStepTrans() override;

    //! base method for solving one transient step 
    virtual void SolveStepTrans() override;
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() override {}

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initializations before execution the SolveStep-method
    virtual void PreStepHarmonic() override;


    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic() override;


    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() override {}


    //----------------------- GetRidOfZeros-----------------------------------
    void SetupGetRidOfZerosActive() override;

    //----------------------- SET/GET METHODS--------------------------------
    
    //! Set actual time
    virtual void SetActTime( const Double actTime ) override;

    //! Set actual frequency
    virtual void SetActFreq( const Double actFreq ) override;

    //! Set actual time / frequency step
    virtual void SetActStep( const UInt actStep ) override;

    //! Set the current time step value
    void SetTimeStep( Double dt );

    //! Set number of time steps
    virtual void SetNumTimeSteps (UInt numTimeStep) override;
    
    //! Set restart time / frequency step
    virtual void SetStartStep (const UInt startStep) override;

  protected:

    //! Finalize structure
    void Finalize();
    
    //! Resort order of single and coupledPDEs
    
    //! This method arranges the SinglePDEs and CoupledPDE(s) in such a way,
    //! that a meaningful solution process is possible. Currently this is
    //! is mostly hard-coded, but in the future we could incorporate
    //! information about coupling quantities and inter-PDE dependencies.
    void ResortPDEOrder();
    
    //! reference to PDE
    IterCoupledPDE &rPDE_;

    //! analysis type of all iteratively coupled PDEs is retrieved
    BasePDE::AnalysisType actAnalysisType_;
    
    //! Paramnode of <iterative>-element of the xml
    PtrParamNode param_;
    
    //! Infonode 
    PtrParamNode info_;
    
    //! Specific info node for convergence
    PtrParamNode convNode_;
    
    //! Density value (key) to integrated value (value)
    std::map<SolutionType, SolutionType> solutionMap_;

    //! Flag, if object is finalized
    bool isFinalized_;
    
    // use user defined PDE order
    bool customReorderPDE_;

    // either start or end with coupled PDEs
    bool endWithCoupledPDEs_;

    // custom PDE order
    std::string PDEorder_;

    // ----------------------------------------------------------------------
    //  Convergence related data
    // ----------------------------------------------------------------------
    
    //! Flag for nonlinear logging
    bool nonLinLogging_;
    
    //! Flag if simulation should be aborted in case of diveregence
    bool stopOnDivergence_;
    
    //! Flag indicating if mechanical displacement shall be treated as a simple norm for convergence or
    //! if an actual change in geometry shall be calculated
    bool justNorm_;
    
    //! Maximum number of iterations per time step
    UInt maxiter_;      

    //! Map, associating solution types with coupling criterions
    std::map<SolutionType, shared_ptr<ConvCriterion> > criterions_;
    
  };

} // end of namespace

#endif

