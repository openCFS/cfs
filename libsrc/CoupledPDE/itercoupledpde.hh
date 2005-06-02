#ifndef FILE_ITERCOUPLEDPDE_2003
#define FILE_ITERCOUPLEDPDE_2003

#include "PDE/basePDE.hh"

namespace CoupledField
{

  // forward class declaration
  class IterSolveStep;
  class StdPDE;
  class PDECoupling;

  //! This class iteratively solve a list of given SinglePDEs 
  class IterCoupledPDE : public BasePDE
  {

    // friend declaration
    friend class IterSolveStep;

  public:

    //! Constructor
    IterCoupledPDE(StdVector<StdPDE*> & PDEs,
                   StdVector<PDECoupling*> & Couplings,
                   std::string sequenceTag); 

    //! Destructor
    ~IterCoupledPDE();

    //! calculates coupling interfaces
    void InitCoupling();
  
    //! write general defines (BCs, loads, etc.) to info-file
    void WriteGeneralPDEdefines();

    Assemble * getPDE_assemble() {
      Error("Get Assemble-Object makes no sense for itercoupledPDE"); 
      return NULL;
    };

    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();
  
    //! Defines, which of the coupled PDEs are currently solved
    //! and which are neglected. This method is mainly needed
    //! for multiSequenceAnalysis, where in different steps
    //! only a subset of all coupled PDEs is solved.
    //! \param pdes (input) Unsorted list of PDEs which
    //! are currently solved 
    void DefineSolvingPDEs(StdVector<StdPDE*> & pdes);

    //! Solve static step
    //void SolveStepStatic(const UInt kstep, const Double asteptime, 
    //           const Boolean updatesysmat);
  
    //! solve transient step
    //void SolveStepTrans(const UInt kstep, const Double asteptime,
    //                  const Boolean updatesysmat);
  
    //! write results in file
    void WriteResultsInFile(const UInt kstep = 0,
                            const Double asteptime = 0.0,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);
  
    //! set time step
    //! \param dt Current time step
    void SetTimeStep(const Double dt);

    //! Do Postprocessing as descriped in conf file
    void PostProcess();
  
  protected:

    //! Write information about coupling interfaces
    //! and coupling setup into ostream
    void WriteCouplingInfo(std::ostream &out);

    //! calculates the norm of a vector
    //Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval);

    UInt maxiter_;                        //!< maximum number of iterations per time step
    StdVector<Double> norms_;              //!< norm of coupling values

    //! pointer to SolveStep classes
    IterSolveStep * solveStep_;

    //! Flag for nonlinear logging
    Boolean nonLinLogging_;
  
    // general PDE parameters
    AnalysisType analysistype_;         //!< type of analysis
    StdVector<StdPDE *> PDEs_;         //!< list of belonging PDEs
    StdVector<PDECoupling*> Couplings_; //!< vector of coupling objects
    UInt NumPDEs_;                   //!< number of PDEs 
    std::string sequenceTag_;           //!< tag for current multisequence step
  
    // 
  
    //! vector of flags indicating if specified
    //! PDE gets solved. The ordering corresponds
    //! to that of PDEs_ in the base class.
    StdVector<Boolean> solvePDE_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class IterCoupledPDE
  //! 
  //! \purpose 
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

} // end of namespace

#endif
