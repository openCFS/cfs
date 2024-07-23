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
    virtual void ResetValues();
    
    //! \copydoc ConvCriterion::StartSampling
    virtual void StartSampling();
    
    //! \copydoc ConvCriterion::StopSampling
    virtual void StopSampling();
    
    //! \copydoc ConvCriterion::GetNorm
    Double GetNorm();
    
    //! \copydoc ConvCriterion::Converged
    bool Converged();
    
    //! \copydoc ConvCriterion::GetSupport
    StdVector<shared_ptr<EntityList> > GetSupport();
    
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
    
    //! Set FeFunction for velocity (needed for geometry update)
    void SetVelFct( shared_ptr<FeFunction<Double> > vel);

    //! Set FeFunction for displacement (needed for geometry update)
    void SetAccFct( shared_ptr<FeFunction<Double> > acc);

    //! Add a new region to 
    void AddRegion(RegionIdType region );
    
    void SetNormFlag( bool justNorm);
    
    //! \copydoc ConvCriterion::ResetValues
    virtual void ResetValues();
    
    //! \copydoc ConvCriterion::StartSampling
    virtual void StartSampling();
    
    //! \copydoc ConvCriterion::StopSampling
    virtual void StopSampling();
    
    //! \copydoc ConvCriterion::GetNorm
    Double GetNorm();
    
    //! \copydoc ConvCriterion::Converged
    bool Converged();
    
    //! \copydoc ConvCriterion::GetSupport
    StdVector<shared_ptr<EntityList> > GetSupport();
    
  protected:
    
    //! Pointer to displacement FeFunction
    shared_ptr<FeFunction<Double> > disp_mech_;

    //! Pointer to displacement FeFunction
    shared_ptr<FeFunction<Double> > disp_smooth_;

    //! Pointer to displacement FeFunction
    shared_ptr<FeFunction<Double> > disp_;

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

    //! Function similar to GetCouplingCoefFct but only returns updatedGeo
    void GetUpdateGeoForPDE( SolutionType type,
                             shared_ptr<EntityList>  list,
                             const std::string& pdeName,
                             bool& updatedGeo );

    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic()  {;};
 
    /** base method for solving one static step */
    virtual void SolveStepStatic();

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic() {;}


    //----------------------- TRANSIENT---------------------------------------

    //! Initialize additional data-structures as needed for the glm
    virtual void InitTimeStepping();
    
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepTrans();


    //! base method for solving one transient step 
    virtual void SolveStepTrans();
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic();


    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic();


    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic(  ) {;};


    //----------------------- GetRidOfZeros-----------------------------------
    void SetupGetRidOfZeros() override;

    //----------------------- SET/GET METHODS--------------------------------
    
    //! Set actual time
    virtual void SetActTime( const Double actTime );

    //! Set actual frequency
    virtual void SetActFreq( const Double actFreq );

    //! Set actual time / frequency step
    virtual void SetActStep( const UInt actStep );

    //! Set the current time step value
    void SetTimeStep( Double dt );

    //! Set number of time steps
    virtual void SetNumTimeSteps (UInt numTimeStep);
    
    //! Set restart time / frequency step
    virtual void SetStartStep (const UInt startStep);

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

