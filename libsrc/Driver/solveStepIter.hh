#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/environment.hh"
#include "PDE/basePDE.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class BaseSolveStep
  {

  public:

    //! Constructor
    BaseSolveStep(BasePDE& apde);

    //! Destructor
    virtual ~BaseSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param level level of grid
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset) {;};
 
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


    //---------------------------- get BasePDE-methods----------------------
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
  };





} // end of namespace
