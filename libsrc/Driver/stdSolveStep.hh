#ifndef FILE_STDSOLVESTEP
#define FILE_STDSOLVESTEP

#include "baseSolveStep.hh"
#include "PDE/StdPDE.hh"
#include "Driver/assemble.hh"

namespace CoupledField
{
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
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset){;}
 
    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic(const Integer kstep, const Double asteptime,
				 const Integer level, const Boolean reset);

    //! solves for one linear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticLin(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset);

    //! solves for one nonlinear static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepStaticNonLin(const Integer kstep, const Double asteptime,
				  const Integer level, const Boolean reset)
    {Error("StepStaticNonLin not implemented!",__FILE__,__LINE__);};

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level) {;};



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean reset);

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean updatesysmat);

    //! solves for one linear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransLin(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean updatesysmat);

    //! solves for one nonlinear transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepTransNonLin(const Integer kstep, const Double asteptime,
				 const Integer level, const Boolean updatesysmat)
    {Error("Nonlinear Transient Step not implemented!",__FILE__,__LINE__);};

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
    */  
    virtual void PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level);

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset);

    //!  base method for solving one harmonic step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepHarmonic(const Integer freqStep, const Double frequency, 
				   Integer level, const Boolean reset);

    //! solves for one linear frequency step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepHarmonicLin(const Integer freqStep, const Double frequency, 
				 Integer level, const Boolean reset);

    //! solves for one nonlinear frequency step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void StepHarmonicNonLin(const Integer freqStep, const Double frequency, 
				    Integer level, const Boolean reset)
    {Error("Harmonic step not implemented!",__FILE__,__LINE__);};


    //!  routine for actions after the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const Integer freqStep, const Double frequency, 
				  Integer level, const Boolean reset) {;};


    //----------------------- helpfull methods--------------------------------------

    //! computes linear part of RHS
    Double SetLinRHS(const Integer level);

    //! stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
		      Double& etaLineSearch, Integer level, Boolean trans=FALSE);

    //! calculates L2-norm of RHS regarding entries due to penalty formulation
    Double RhsL2Norm(Vector<Double>& stdVec);

    //! returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);


    //! Write nonlin iteration norms to info-file
    void WriteClaNlNorms(const Integer iterationCounter,
			 const Double residualL2Norm,
			 const Double extForcesL2Norm, const Double residualErr, 
			 const Double solIncrL2Norm, const Double actSolL2Norm, 
			 const Double incrementalErr);

  protected:


    //------------- storage vectors for nonlinear analysis --------------
    Vector<Double> RhsLinVal_; //!< external forces (for nonlin simulations)


    //---------------------------- get StdPDE-methods----------------------
    //! set Dirichlet BCs
    void SetBCs( const Integer level,  const Double time );


    //! fetches coordinates to element nodes
    void GetElemCoords(const StdVector< Integer > connect,
			       Matrix< Double > &coordMat,
			       const Integer level);

    //! returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    //! returns the vector of time derivative of the solution belonging to all nodes 
    //!  of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    //! returns the vector of 2nd time derivative of the solution belonging to all nodes 
    //! of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);


    //-------------------------------- Pointers to (Copies of) StdPDE -------------------

    StdPDE& PDE_;                   //!< reference to PDE
    std::string pdename_;            //!< name of PDE 
    Integer numPDENodes_;            //!< number of nodes belonging to the PDE
    Boolean isaxi_;                  //!< TRUE: axisymmetric problem
    StdVector<std::string> subdoms_; //!< subdomain-levels belonging to PDE
    MaterialData * materialData_;    //!< pointer to material data of PDE
    Grid * ptgrid_;                  //!< pointer to grid object
    BaseSystem* algsys_;             //!< pointer to algsys object
    BaseNodeStoreSol * sol_;         //!< pointer to solution object
    NodeEQN * eqnData_;              //!< pointer to equation object
    Assemble * assemble_;            //!< pointer to assemble object  
    Integer actlevel_;               //!< current level 


    Double * matrix_factor_;       //!< factors for compution effective system matrix
    TimeStepping * TS_alg_;        //!< pointer to time-stepping object
    Double lasttimecalc_;          //!< last time on which we have calculated solution
    Integer laststepcalc_;         //!< Number of last timestep on which we have calculated 
                                   //!< our solution
    Boolean recalc_;               //!< flag indicating reassembling of system matrix

    Boolean pdeIsCoupled_;         //!< TRUE, if PDE is coupled to other ones
    Boolean firstTimeStepStatic_;  //!< needed for coupled, iterative methods
    Integer* iterCoupledCounter_;  //!< iteration counter for coupled PDE solution process

    std::string lineSearch_;   //!< switch for lineSearch
    Boolean nonLin_;           //!< flag for nonlinear calculations
    Boolean geoUpdate_;        //!< flag for geometric update
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    Integer nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    Boolean nonLinLogging_;    //!< log progress of non-linear iterations
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)

    DampingType dampingType_;  //!< damping type of PDE
    Vector<Double> solIncr_;      //! needed in iterative coupled computation 
    Vector<Double> actSol_;       //! needed in iterative coupled computation 
    Boolean isIncrFormulation_;   //! checks, if we have for the coupling a incremental solution

  };

} // end of namespace

#endif

