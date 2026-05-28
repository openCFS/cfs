#ifndef FILE_STDSOLVESTEP
#define FILE_STDSOLVESTEP

#include <map>
#include <fstream>

#include "BaseSolveStep.hh"

#include "Driver/Assemble.hh"
#include "MatVec/SBM_Vector.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/BaseMaterial.hh"
#include "DataInOut/ResultHandler.hh"
#include "Utils/Timer.hh"
#include "Domain/Domain.hh"


namespace CoupledField
{
  // forward class declarations
  class StdPDE;
  class WriteResults;
  struct ResultInfo;
  class SingleDriver;
  class IDBC_Handler;
  class BaseIDBC_Handler;
  class FeSpace;
  class SolStrategy;
  class MHTimeFreqResult;
  class MathParser;

  //  class Domain;
  
  //! Derived class for step-wise solving of StdPDEs
  class StdSolveStep : public BaseSolveStep {

  public:

    // public typedefs
    typedef StdVector<shared_ptr<ResultInfo> > ResultList;

    //! Constructor
    StdSolveStep(StdPDE& apde);

    //! Destructor
    virtual ~StdSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic();
 
    /** base method for solving one static step */
    virtual void SolveStepStatic();

    /** @see SolveStepStatic() */ 
    virtual void StepStaticLin();

    //! solves for one nonlinear static step: incremental formulation 
    virtual void StepStaticNonLin();

    /** Solve a simple fixed point nonlinear static problem for nonlinear RHS */
    void StepStaticNonLinFixedPointSimple();

    //! routine for actions after the SolveStep-method
    virtual void PostStepStatic();

    //----------------------- TRANSIENT---------------------------------------

    //! Initialize additional data-structures as needed for the glm
    virtual void InitTimeStepping();

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepTrans();

    //! routine for computing a predictor step
    // neede in case of FSI-Iterative-Coupling
    //virtual void PredictorStep(){;};

    //! base method for solving one transient step 
    virtual void SolveStepTrans();

    //! solves for one linear transient step
    virtual void StepTransLin();

    //! solves for one nonlinear transient step: incremental formulation 
    virtual void StepTransNonLin();

    //! solves for one nonlinear transient step: total formulation 
    virtual void StepTransNonLinTotal();

    //! solves for one hysteretic transient step: direct quasi-Newton formulation
    virtual void StepTransHyst();

    //! solves for one nonlinear transient step 
    //! consideres material nonlinearities in direct coupled PDEs
    void StepTransNonLinMaterial() {REFACTOR;};
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() {};

    //----------------------- HARMONIC AND MULTIHARMONIC -------------------------
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic();

    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic();

    void SolveStepHarmonic25D(Double ExcitationFreq, Double WaveNumFreq);
    
    //! solves for one linear frequency step 
    virtual void StepHarmonicLin();

    //! solves for one nonlinear frequency step 
    virtual void StepHarmonicNonLin();
    
    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() {;};
    
    //! same as GetSoltionVal and GetRHSVal but only in the
    //! multiharmonic case and it's triggered by the MultiHarmonicDriver
    //! in the SolveProblem() method
    virtual void GetSolutionValMultHarm(const UInt& h);
    virtual void GetRHSValMultHarm(const UInt& h);

    //----------------------- EIGENFREQUENCY ----------------------------------

    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Double> & frequencies, Vector<Double> & errBounds,
                               UInt numFreq, double shift, bool sort);

    //! Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Complex> & frequencies, Vector<Double> & errBounds,
                               UInt numFreq, double shift, bool sort, bool bloch);

    //! Calculate the Eigenfrequencies in an interval [minVal,maxVal]
    UInt CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds, Double minVal, Double maxVal);

    void CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal );

    //! Calculate the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    void GetEigenMode( UInt numMode, bool right=true );
    
    //----------------------- helpfull methods--------------------------------------

    /** The Assemle object contains the bilinear forms */
    Assemble* GetAssemble() { return assemble_; }

    AlgebraicSys * GetAlgSys() { return algsys_; }

    //! Set the current time step
    void SetTimeStep( Double dt );

    //! computes linear part of RHS
    Double SetLinRHS(Double loadFactor,bool nonlin = false, bool multiharmonic = false);

    //! computes ldelta inear part of RHS; in case of sub stepping
    UInt SetDeltaLinRHS();

    //! does a line search and returns the optimal residual norm
    Double LineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol, 
                      Double& etaLineSearch, bool trans=false);

    //! does a line search: It gives the optimal line search parameter. 
    //! This minimizes the Energy functional in the Newton direction.
    //! To accelerate the process just a set of eta \in Eta is tried out.
    Double LineSearchMinEnergyHeuristic(SBM_Vector& solIncrement, SBM_Vector& actSol, 
                      Double& etaLineSearch, bool trans=false);

    //! does a line search: It gives the optimal line search parameter. 
    //! This minimizes the Energy functional in the Newton direction.
    //! Here the minimization is exactly solved by Brent's method.
    Double LineSearchMinEnergy(SBM_Vector& solIncrement, SBM_Vector& actSol, 
                      Double& etaLineSearch, bool trans=false);

    //! computes the derivative of the energy functional w.r.t the line search parameter eta.
    Double GetEnergyDerivativeEta(SBM_Vector& solIncrement, SBM_Vector& actSol, Double eta);

    //! solves a one-dimensional root-finding problem by Brent's method
    Double BrentMethod(SBM_Vector& solIncrement, SBM_Vector& actSol, double a, double b);

    //! does a line search for multiharmonic analysis and returns the optimal residual norm
    Double LineSearchMultHarm(const SBM_Vector& solIncrement, SBM_Vector& actSol,
                      Double& etaLineSearch, MHTimeFreqResult& ftRes);


    //! does a line search and returns the optimal residual norm
    Double LineSearchMag(SBM_Vector& solIncrement, SBM_Vector& actSol,
                      Double& etaLineSearch, bool trans=false);

    //! does a line search and returns the optimal residual norm
    Double LineSearchMaterial(SBM_Vector& solIncrement, 
                              SBM_Vector& actSol, 
                              Double& etaLineSearch, Double& RHSLin2Norm,
                              bool trans=false);

    //! checks if getRidOfZeros should be actually used and if yes, defines everything accordingly and warns the user about it
    void SetupGetRidOfZerosActive();
    
    void SetSolveVecZero(){
      solVec_.Init();
    }

  protected:
    
    // ========================================================
    //  Helper Methods
    // ========================================================
    
    //! Read nonlinear data from pdenode 
    virtual void ReadNonLinData();
    
    /** Write to info.xml and optional non-lin logging file.
     *  Checks programOpt->DoDetail()  */
    void OutputNonLinIterInfo(const std::string& pdeName, UInt solStep, UInt iterationCounter,
        Double residualErr, Double incrementalErr, double etaLineSearch, int coupledIterStep = -1);
    
    //------------- storage vectors for nonlinear analysis --------------
    //Vector<Double> RhsLinVal_; //!< external forces (for nonlin simulations)
    SBM_Vector oldRhsLinVal_; //!< external forces (for nonlin simulations)
    SBM_Vector tmpOldRhsLinVal_; //!< external forces (for nonlin simulations)
    SBM_Vector DeltaRhsLinVal_; //!< external forces (for nonlin simulations)

    //------------- storage vectors for nonlinear analysis --------------
    SBM_Vector RhsLinVal_; //!< external forces (for nonlin simulations)


    //-------------------------------- Pointers to (Copies of) StdPDE -------------------
    StdPDE& PDE_;                   //!< reference to PDE
    std::string pdename_;            //!< name of PDE 
    bool isaxi_;                  //!< true: axisymmetric problem
    StdVector<RegionIdType> subdoms_;//!< subdomain-levels belonging to PDE

    /** to summarize .info.xml output about nonlineariy */
    PtrParamNode nonLinInfo_;

    //! Pointer to material data of PDE
    std::map<RegionIdType, BaseMaterial*> materialData_;  

    Grid * ptgrid_;                  //!< pointer to grid object
    AlgebraicSys* algsys_;             //!< pointer to algsys object
    ResultList results_;
    Assemble * assemble_;            //!< pointer to assemble object

    
    //! Pointer to solution strategy object
    shared_ptr<SolStrategy> solStrat_;
    
    //! factors for computingn effective system matrix for each feFunction
    std::map<FeFctIdType, std::map<FEMatrixType,Double> > matrix_factor_;   
    
    //! factors for computingn effective system matrix for each feFunction
    std::map<FeFctIdType, std::map<FEMatrixType,Complex> > matrix_factor_Complex_;

                                   //!< our solution
    bool recalc_;               //!< flag indicating reassembling of system matrix

    std::string lineSearch_;   //!< switch for lineSearch
    Double lineSearchTolerance_;  //!< tolerance for line search convergence (default: 1e-3)
    UInt lineSearchMaxIter_;      //!< max iterations for line search (default: 1000)
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isHyst_;           //!< flag for hystersis modeling
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    Double lastError_; // So it stops if residum gets bigger!
    UInt consideredH_=1; // So it stops if residum gets bigger!

    Double minValidValue_;     //! stopping if any value in the region exceeds value
    Double maxValidValue_;     //! stopping if any value in the region exceeds value
    SolutionType solutionLimit_; //! solution type for which a limit is set
    RegionIdType solutionLimitReg_; //! region in which to check the min/max values for non convergence

    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    bool nonLinLogging_;    //!< log progress of non-linear iterations
    UInt minLoggingToTerminal_;
    bool nonLinTotalFormulation_;   //!< flag for total or incremental NL formulation
    bool abortOnMaxIter_; //!< flag for aborting simulation if maximum number of iterations is hit

    bool useStagnationDetection_ = false; //!< enable stagnation detection for nonlinear iteration
    UInt stagnationWindowSize_ = 5; //!< number of iterations used to inspect stagnation
    Double stagnationTolerance_ = 0.1; //!< relative improvement threshold for stagnation

    //! map for each region the type of nonlinearity
    std::map<RegionIdType, StdVector<NonLinType> > regionNonLinTypes_;

    SBM_Vector solIncr_;   //! needed in iterative coupled computation 
    SBM_Vector actSol_;    //! needed in iterative coupled computation 

    //! Vector containing all solution vectors of the FE-functions
    SBM_Vector solVec_;
    
    //! Vector containing rhs
    SBM_Vector rhsVec_;

    //! Vector containing the rhs for the current stage based on the scheme
    //! Vector containing the rhs for the current stage based on the scheme
    //! TODO: This can be obtimized if the time schemes write their rhs parts directly to the Algebraic system
    SBM_Vector stageRHS_;

    //! Map Storing FeSpaces for each solution type of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;
    
    //! Map Storing FeSpaces for each solution type of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions_;

    boost::shared_ptr<Timer> static_non_lin_step_timer_;

    std::ofstream logFile_;
    unsigned int mHandle_;
    MathParser* mParser_;

    //! Bool to check if getRidOfZeros should be used
    // we assume that we should not perform GetRidOfZeros(), if this is not the case, this will be set afterwards
    bool useGetRidOfZeros_ = false;

    
    //! Tolerance used to determine if an entry is zero or not in getRidOfZeros (the same is defined in XML schema)
    Double getRidOfZerosTol_ = 1e-20;

    //! Bool to check the setup of the evaluation of SetupGetRidOfZerosActive() only once
    bool setupGetRidOfZerosDone_ = true;

private:
  void AssembleMH(const UInt& N, const UInt& M, const bool onlyDiagBlocks = false);

  void EvaluateNonlinearity(MHTimeFreqResult& ftRes, const SBM_Vector& actSol);


  //! Vector containing all solution vectors for all harmonics
  //! in a multiharmonic analysis. We need this vector because
  //! solVec_ is only used to pass certain harmonics back to
  //! the PDE
  SBM_Vector solVecMH_;
};

} // end of namespace

#endif
