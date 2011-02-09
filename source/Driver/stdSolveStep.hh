// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STDSOLVESTEP
#define FILE_STDSOLVESTEP

#include <map>

#include "assemble.hh"
#include "baseSolveStep.hh"
#include "MatVec/vector.hh"
#include "Utils/hysteresis.hh"
#include "Materials/baseMaterial.hh"
#include "DataInOut/resultHandler.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/domain.hh"

namespace CoupledField
{
  // forward class declarations
  class StdPDE;
  class BaseNodeStoreSol;
  class TimeStepping;
  class WriteResults;
  class EqnMap;
  struct ResultInfo;
  class SingleDriver;
  class IDBC_Handler;
  class BaseIDBC_Handler;

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
    virtual void SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    /** @see SolveStepStatic() */ 
    virtual void StepStaticLin(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! solves for one nonlinear static step 
    virtual void StepStaticNonLin(PtrParamNode analysis_id);
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepStatic();

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepTrans();

    //! routine for computing a predictor step
    // neede in case of FSI-Iterative-Coupling
    //virtual void PredictorStep(){;};

    //! base method for solving one transient step 
    virtual void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! solves for one linear transient step 
    virtual void StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! solves for one nonlinear transient step 
    virtual void StepTransNonLin(PtrParamNode analysis_id);

    //! solves for one nonlinear transient step 
    //! consideres material nonlinearities in direct coupled PDEs
    void StepTransNonLinMaterial(PtrParamNode analysis_id);

    //! solves for one nonlinear transient step 
    //! consideres hystreresis nonlinearities in direct coupled PDEs
    virtual void StepTransNonLinHysteresis(PtrParamNode analysis_id);
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans();

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic();

    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic(PtrParamNode analysis_id);
    
    //! solves for one linear frequency step 
    virtual void StepHarmonicLin(PtrParamNode analysis_id);

    //! solves for one nonlinear frequency step 
    virtual void StepHarmonicNonLin(PtrParamNode analysis_id)
    {EXCEPTION("Harmonic step not implemented!");};
    
    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() {;};
    

    //----------------------- HARMONIC ---------------------------------------

    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Double> & frequencies,
                               Vector<Double> & errBounds,
                               UInt numFreq, Double shift );

    //! Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    UInt CalcEigenFrequencies( Vector<Complex> & frequencies,
                               Vector<Double> & errBounds,
                               UInt numFreq, Double shift );

    //! Calculate the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    void CalcEigenMode( UInt numMode );
    
    
    //----------------------- helpfull methods--------------------------------------

    /** The Assemle opject contains the bilinear forms */
    Assemble* GetAssemble() { return assemble_; } 

    //! Set the current time step
    void SetTimeStep( Double dt );

    //! computes linear part of RHS
    Double SetLinRHS(Double loadFactor);

    //! computes ldelta inear part of RHS; in case of sub stepping
    UInt SetDeltaLinRHS();

    //! stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                      Double& etaLineSearch, bool trans=false);

    //! does a line search and returns the optimal residual norm
    Double LineSearchMaterial(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                              Double& etaLineSearch, Double& RHSLin2Norm,
                              bool trans=false);


    //! returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);

    //! set the identification tag of the PDE
    void SetPDEId( const PdeIdType pdeId )
    { pdeId1_ = pdeId;};

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
    
    void ReInit();

  protected:


    //! Read nonlinear data from pdenode 
    virtual void ReadNonLinData();

    //------------- storage vectors for nonlinear analysis --------------
    //Vector<Double> RhsLinVal_; //!< external forces (for nonlin simulations)
    Vector<Double> oldRhsLinVal_; //!< external forces (for nonlin simulations)
    Vector<Double> tmpOldRhsLinVal_; //!< external forces (for nonlin simulations)
    Vector<Double> DeltaRhsLinVal_; //!< external forces (for nonlin simulations)

    //------------- storage vectors for nonlinear analysis --------------
    Vector<Double> RhsLinVal_; //!< external forces (for nonlin simulations)


    //-------------------------------- Pointers to (Copies of) StdPDE -------------------
    StdPDE& PDE_;                   //!< reference to PDE
    std::string pdename_;            //!< name of PDE 
    UInt numPDENodes_;            //!< number of nodes belonging to the PDE
    UInt numPDEElems_;            //!< number of elements belonging to PDE
    bool isaxi_;                  //!< true: axisymmetric problem
    StdVector<RegionIdType> subdoms_;//!< subdomain-levels belonging to PDE

    //! Pointer to material data of PDE
    std::map<RegionIdType, BaseMaterial*> materialData_;  

    Grid * ptgrid_;                  //!< pointer to grid object
    BaseSystem* algsys_;             //!< pointer to algsys object
    BaseNodeStoreSol * sol_;         //!< pointer to solution object
    shared_ptr<EqnMap> eqnMap_;
    ResultList results_;
    Assemble * assemble_;            //!< pointer to assemble object  
    //! factors for computingn effective system matrix
    std::map<FEMatrixType,Double> matrix_factor_;   
    
    TimeStepping * TS_alg_;        //!< pointer to time-stepping object
                                   //!< our solution
    bool recalc_;               //!< flag indicating reassembling of system matrix

    std::string lineSearch_;   //!< switch for lineSearch
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isHyst_;           //!< flag for hystersis modeling
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    bool nonLinLogging_;    //!< log progress of non-linear iterations
    //! map for each region the type of nonlinearity
    std::map<RegionIdType, NonLinType> regionNonLinType_;

    //! size of rhs and algsys vector
    UInt numEqns_;
    
    Vector<Double> solIncr_;   //! needed in iterative coupled computation 
    Vector<Double> actSol_;    //! needed in iterative coupled computation 

    //hysteresis operator;    
    StdVector<Hysteresis *> hyst_;

    //! \todo Currently only two pdes can couple. This has to be extended
    //! for the general case
    //! Identification tag for first PDE
    PdeIdType pdeId1_;

    //! Identification tag for second PDE (coupled case)
    PdeIdType pdeId2_;

    MathParser::HandleType mHandle_;
    MathParser* mParser_;
  };

} // end of namespace

#endif

