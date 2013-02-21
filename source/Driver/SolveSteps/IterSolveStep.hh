#ifndef ITERSOLVESTEP_HH
#define ITERSOLVESTEP_HH

#include "BaseSolveStep.hh"

#include "PDE/BasePDE.hh"
#include "MatVec/SingleVector.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField
{
  // forward class declarations
  class IterCoupledPDE;

  //! Derived class for step-wise solving of iterative coupled StdPDEs
  class IterSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    IterSolveStep(IterCoupledPDE& apde, PtrParamNode );

    //! Destructor
    virtual ~IterSolveStep();

    //! Initialize struct
    void Init();
    
    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic()  {;};
 
    /** base method for solving one static step */
    virtual void SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic() {;}


    //----------------------- TRANSIENT---------------------------------------

    //! Initialize additional data-structures as needed for the glm
    virtual void InitTimeStepping();
    
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepTrans();


    //! base method for solving one transient step 
    virtual void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic();


    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic(PtrParamNode analysis_id);


    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic(  ) {;};


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

    //! calculates the norm of a vector
    Double CalcNorm(NormType normtype, SingleVector & val, SingleVector & oldval); 
    
    //! Update grid coordinates
    void UpdateGeometry();
    
    //! reference to PDE
    IterCoupledPDE &rPDE_;

    //! analysis type of all iteratively coupled PDEs is retrieved
    BasePDE::AnalysisType actAnalysisType_;
    
    //! Paramnode of <iterative>-element of the xml
    PtrParamNode param_;
    
    //! FeFunction for nodal displacements
    
    //! Pointer to nodal displacements, used for grid update at 
    //! the end of an iteration.
    shared_ptr<FeFunction<Double> > disp_;
    
    //! Set containig updated regions
    std::set<RegionIdType> updatedRegions_;

  };

} // end of namespace

#endif

