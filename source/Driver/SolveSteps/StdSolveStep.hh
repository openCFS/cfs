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
#include "Utils/mathParser/mathParser.hh"
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

    //! solves for one nonlinear transient step 
    //! consideres material nonlinearities in direct coupled PDEs
    void StepTransNonLinMaterial();
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans();

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic();

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepMultiHarmonic(StdVector<Double> harmFreq);

    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic();
    
    //! solves for one linear frequency step 
    virtual void StepHarmonicLin();

    //! solves for one nonlinear frequency step 
    virtual void StepHarmonicNonLin();
    
    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() {;};
    
    //----------------------- HARMONIC ---------------------------------------

    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Double> & frequencies, Vector<Double> & errBounds,
                               UInt numFreq, double shift, bool sort);

    //! Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Complex> & frequencies, Vector<Double> & errBounds,
                               UInt numFreq, double shift, bool sort, bool bloch);

    //! Calculate the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    void GetEigenMode( UInt numMode );
    
    //----------------------- HYSTERESIS -------------------------------------
    //! solves for one nonlinear transient step
    //! consideres hystreresis nonlinearities in direct coupled PDEs
    virtual void StepTransNonLinHysteresis();
    virtual void StepTransNonLinHysteresisTotal();

    void SetLastItOrLastTSSBMVectors(bool lastTS);
    /*!
     * Helper funciton for setting up the equation system during
     *            StepTransNonLinHysteresis()
     * Background: During the solve step, the matrices and the rhs have to be
     *  assembled multiple times during linesearch, for the calculation of the
     *  residual error and of course to get a system to be solved;
     *  for simplification, encapsulate that sequence of function calls
     *  in a separate function
     */
    /*
     * for residual computation we need a slightly different version -> see .cc file
     */
    virtual void CalcResidualAndConfigSystemForHysteresis(SBM_Vector& oldSolution,SBM_Vector& solIncrement, SBM_Vector& stageSol, Double usedEta,
                                                          UInt stage, UInt callingCnt, UInt evalVersion, bool trans, bool forceReevaluation,
                                                          bool skipReassembly, bool debugOutput, bool reset);

    //! does a line search and returns the optimal residual norm
    Double LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion, UInt callingCnt,
                      bool trans=false, bool performLineSearch=true, bool forceReevaluation=false, bool debugOutput=false, bool reset=false,UInt allowedSteps=5);
    
    //----------------------- helpfull methods--------------------------------------

    /** The Assemle opject contains the bilinear forms */
    Assemble* GetAssemble() { return assemble_; } 

    //! Set the current time step
    void SetTimeStep( Double dt );

    //! computes linear part of RHS
    Double SetLinRHS(Double loadFactor,bool nonlin = false);

    //! computes ldelta inear part of RHS; in case of sub stepping
    UInt SetDeltaLinRHS();

    //! does a line search and returns the optimal residual norm
    Double LineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol, 
                      Double& etaLineSearch, bool trans=false);

    //! does a line search and returns the optimal residual norm
    Double LineSearchMag(SBM_Vector& solIncrement, SBM_Vector& actSol,
                      Double& etaLineSearch, bool trans=false);

    //! does a line search and returns the optimal residual norm
    Double LineSearchMaterial(SBM_Vector& solIncrement, 
                              SBM_Vector& actSol, 
                              Double& etaLineSearch, Double& RHSLin2Norm,
                              bool trans=false);


    //! returns the hysteresis operator
    Hysteresis * GetHystOperator(UInt iSD) {
      return hyst_[iSD];
    };
    


  protected:
    
    // ========================================================
    //  Helper Methods
    // ========================================================
    
    //! Read nonlinear data from pdenode 
    virtual void ReadNonLinData();
    
    virtual void WriteNonLinIterToInfoXML(const std::string& pdeName, 
                                          const UInt solStep,
                                          const UInt iterationCounter,
                                          const Double residualErr, 
                                          const Double incrementalErr, 
                                          double etaLineSearch=0.0);

    virtual void WriteNonLinIterToInfoXML(const std::string& pdeName, 
                                          const UInt coupledIterStep,
                                          const UInt solStep,
                                          const UInt iterationCounter,
                                          const Double residualErr, 
                                          const Double incrementalErr, 
                                          double etaLineSearch=0.0);
    

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
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isHyst_;           //!< flag for hystersis modeling
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    Double minValidValue_;     //! stopping if any value in the region exceeds value
    Double maxValidValue_;     //! stopping if any value in the region exceeds value
    SolutionType solutionLimit_; //! solution type for which a limit is set
    RegionIdType solutionLimitReg_; //! region in which to check the min/max values for non convergence

    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    bool nonLinLogging_;    //!< log progress of non-linear iterations
    bool nonLinTotalFormulation_;   //!< flag for total or incremental NL formulation
    bool abortOnMaxIter_; //!< flag for aborting simulation if maximum number of iterations is hit

    //! map for each region the type of nonlinearity
    std::map<RegionIdType, StdVector<NonLinType> > regionNonLinTypes_;

    SBM_Vector solIncr_;   //! needed in iterative coupled computation 
    SBM_Vector actSol_;    //! needed in iterative coupled computation 

    //! Vector containing all solution vectors of the FE-functions
    SBM_Vector solVec_;
    
    //! Vector containing rhs
    SBM_Vector rhsVec_;

    // zero-vector, used in multiharmonic analysis
    Vector<Complex> zVec_;

    //! Vectors used for NonLinHysteresis
    // current > current timestep and iteration
    SBM_Vector currentLinRhsVec_;
    SBM_Vector currentNonLinRhsVec_;
    SBM_Vector currentResVec_;
    SBM_Vector currentRHSload_partial_; // for full stepping only;
    // stores f_currentTS - f_lastTS - f_nonlin_lastTS; needs to be evaluated only during first iteration
    // during each iteration f_nonlin_currentTS has to be added to get full RHSload

    // + solVec which is defined above

    // oldTS > values after last iteration of previous TS
    // to be stored after iteration suceeded
    SBM_Vector oldTSLinRhsVec_;
    SBM_Vector oldTSNonLinRhsVec_;
    SBM_Vector oldTSSolVec_;

    // oldIt > values of the current TS but from previous It
    // during first iterartion of a new TS, these vectors contain the values
    // after the last iteration of the previous TS (similar as oldTS...)
    SBM_Vector oldItNonLinRhsVec_;
    SBM_Vector oldItResVec_;

    //! Additional flags and parameter for hyst
    UInt evalVersion_;
    bool forceReevaluation_;
    bool debugOutput_;
    UInt remainingEvalParameter_;

    //! Vector containing the harmonic-frequencies (negative, 0, positive)
    StdVector<Double> harmFreq_;

    //! Vector containing the rhs for the current stage based on the scheme
    //! TODO: This can be obtimized if the time schemes write their rhs parts directly to the Algebraic system
    SBM_Vector stageRHS_;

    //! Map Storing FeSpaces for each solution type of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;
    
    //! Map Storing FeSpaces for each solution type of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions_;
    //hysteresis operator;    
    StdVector<Hysteresis *> hyst_;


    std::ofstream logFile_;
    MathParser::HandleType mHandle_;
    MathParser* mParser_;
  };

} // end of namespace

#endif

