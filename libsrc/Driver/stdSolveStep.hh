#ifndef FILE_STDSOLVESTEP
#define FILE_STDSOLVESTEP

#include <map>

#include "baseSolveStep.hh"
#include "Utils/vector.hh"
#include "Utils/hysteresis.hh"

namespace CoupledField
{
  // forward class declarations
  class StdPDE;
  class Assemble;
  class BaseNodeStoreSol;
  class NodeEQN;
  class TimeStepping;
  
  //! Derived class for step-wise solving of StdPDEs

  class StdSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    StdSolveStep(StdPDE& apde);

    //! Destructor
    virtual ~StdSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    //!\param reset TRUE: perfrom new assembly, etc
    virtual void PreStepStatic( const Boolean reset ){;}
 
    //! base method for solving one static step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepStatic( const Boolean reset );

    //! solves for one linear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepStaticLin( const Boolean reset );

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepStaticNonLin( const Boolean reset )
    {Error("StepStaticNonLin not implemented!",__FILE__,__LINE__);};
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepStatic() {;};



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepTrans( const Boolean reset );
    
    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepTrans( const Boolean updatesysmat );

    //! base method for solving one transient step with slicing method 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepTrans4Slice(const Boolean reset);

    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepTransLin( const Boolean updatesysmat );

    //! solves for one nonlinear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepTransNonLin( const Boolean updatesysmat )
    {Error("Nonlinear Transient Step not implemented!",__FILE__,__LINE__);};
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans();

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepHarmonic( const Boolean reset );

    //!  base method for solving one harmonic step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepHarmonic( const Boolean reset );
    
    //! solves for one linear frequency step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepHarmonicLin( const Boolean reset );

    //! solves for one nonlinear frequency step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepHarmonicNonLin( const Boolean reset )
    {Error("Harmonic step not implemented!",__FILE__,__LINE__);};
    
    //!  routine for actions after the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PostStepHarmonic( const Boolean reset ) {;};


    //! 
    virtual void WriteResults(UInt actStep, Double actTime) {;};

    //----------------------- helpfull methods--------------------------------------

    //! Set the current time step
    void SetTimeStep( Double dt );

    //! computes linear part of RHS
    Double SetLinRHS(Double loadFactor);

    //! stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                      Double& etaLineSearch, Boolean trans=FALSE);

    //! calculates L2-norm of RHS regarding entries due to penalty formulation
    Double RhsL2Norm(Vector<Double>& stdVec);

    //! returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);

    //! set the identification tag of the PDE
    void SetPDEId( const PdeIdType pdeId );

    //! Write nonlin iteration norms to info-file
    void WriteClaNlNorms(const UInt iterationCounter,
                         const Double residualL2Norm,
                         const Double extForcesL2Norm, const Double residualErr, 
                         const Double solIncrL2Norm, const Double actSolL2Norm, 
                         const Double incrementalErr);

    //! returns the hysteresis operator
    Hysteresis * GetHystOperator(UInt iSD) {
      return hyst_[iSD];
    };

    //! transform solution and derivatives due to slicing technique
    virtual void TransformSol4Slice(UInt & shiftFactor, UInt & nodeShift,
		UInt & elemgrid, Double &  meshsize, const UInt flag);

    //! save solution of special nodes
    virtual void SaveNodes(const UInt shiftFactor, const Double timeStep,
		   const UInt numShift, const Integer nodeShift, 
		   const UInt maxnumelemz_);

  protected:


    //------------- storage vectors for nonlinear analysis --------------
    Vector<Double> RhsLinVal_; //!< external forces (for nonlin simulations)


    //---------------------------- get StdPDE-methods----------------------
    //! set Dirichlet BCs
    void SetBCs( const Double time );


    //! fetches coordinates to element nodes
    void GetElemCoords(const StdVector< UInt > connect,
                       Matrix< Double > &coordMat );

    //! returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);

    //! returns the vector of time derivative of the solution belonging to all nodes 
    //!  of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);

    //! returns the vector of 2nd time derivative of the solution belonging to all nodes 
    //! of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);

    //-------------------------------- Pointers to (Copies of) StdPDE -------------------

    StdPDE& PDE_;                   //!< reference to PDE
    std::string pdename_;            //!< name of PDE 
    UInt numPDENodes_;            //!< number of nodes belonging to the PDE
    UInt numPDEElems_;            //!< number of elements belonging to PDE
    Boolean isaxi_;                  //!< TRUE: axisymmetric problem
    StdVector<RegionIdType> subdoms_;//!< subdomain-levels belonging to PDE
    MaterialData * materialData_;    //!< pointer to material data of PDE
    Grid * ptgrid_;                  //!< pointer to grid object
    BaseSystem* algsys_;             //!< pointer to algsys object
    BaseNodeStoreSol * sol_;         //!< pointer to solution object
    NodeEQN * eqnData_;              //!< pointer to equation object
    Assemble * assemble_;            //!< pointer to assemble object  

    //! factors for computingn effective system matrix
    std::map<FEMatrixType,Double> matrix_factor_;   
    
    TimeStepping * TS_alg_;        //!< pointer to time-stepping object
                                   //!< our solution
    Boolean recalc_;               //!< flag indicating reassembling of system matrix

    Boolean isIterCoupled_;         //!< TRUE, if PDE is coupled to other ones iteratively
    Boolean firstTimeStepStatic_;  //!< needed for coupled, iterative methods
    UInt* iterCoupledCounter_;  //!< iteration counter for coupled PDE solution process

    std::string lineSearch_;   //!< switch for lineSearch
    Boolean nonLin_;           //!< flag for nonlinear calculations
    Boolean isHyst_;           //!< flag for hystersis modeling
    Boolean geoUpdate_;        //!< flag for geometric update
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    Boolean nonLinLogging_;    //!< log progress of non-linear iterations
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)

    Vector<Double> solIncr_;   //! needed in iterative coupled computation 
    Vector<Double> actSol_;    //! needed in iterative coupled computation 
    Boolean isIncrFormulation_;//! checks, if we have for the coupling a incremental solution


    //hysteresis operator;    
    StdVector<Hysteresis *> hyst_;

    //! \todo Currently only two pdes can couple. This has to be extended
    //! for the general case
    //! Identification tag for first PDE
    PdeIdType pdeId1_;

    //! Identification tag for second PDE (coupled case)
    PdeIdType pdeId2_;

    //! counts the number of resets
    UInt numReset_;

  };

} // end of namespace

#endif

