#ifndef FILE_MULTISEQUENCEDRIVER_2004
#define FILE_MULTISEQUENCEDRIVER_2004

#include <map>
#include "BaseDriver.hh"
#include "SingleDriver.hh"
#include "Utils/StdVector.hh"
#include "DataInOut/SimState.hh"

namespace CoupledField
{

  // forward class declaration
  class BasePDE;
  class SinglePDE;

  //! Driver class for multi sequence simulations

  //! Driver class for multi sequence simulations.
  //! Multo sequence here means, that either different types
  //! of simulation (static, transient) can be computed one after
  //! another (e.g. prestressed materials with addtitional source
  //! terms) or that mutliple instances of the same anlysis-type can be
  //! computed (e.g. pule-echo simulation in transient case)

  class MultiSequenceDriver : public BaseDriver
  {
  public:
  
    /** @param keep in the optimization case we store the drivers and pdes instead of deleting them for a new sequence step */
    MultiSequenceDriver(shared_ptr<SimState> state, Domain* domain,
                        PtrParamNode paramNode, PtrParamNode infoNode, bool keep = false);

    virtual ~MultiSequenceDriver();

    //! Initialization method
    void Init(bool restart);
    
    //! Change the current sequenceStep. Calls SetupStep() */
    void SetSequenceStep(UInt sequenceStep);

    //! main method, where time-stepping is implemented. 
    //! it is for transient and static problem
    void SolveProblem();

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename );

    //! Return current singleDriver 
    SingleDriver* GetSingleDriver() { return actDriver_; }

    /** implement abstract identification class */ 
    DriverClass GetDriverClass() { return MULTI_SEQUENCE_DRIVER; };

    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() { return actDriver_->IsComplex(); };

    unsigned int GetNumberOfSequenceSteps() const { return numSteps_; }

    /** required for optimization of multi sequence steps. Creates a copy otherwise we would have to deal with const maps :( */
    std::map<UInt, BasePDE::AnalysisType> GetAnalyisPerStep() const { return analysisPerStep_; }

    /** small service for Optimization. The ParamNodes for the steps. CFS reads the steps one after another
     * but we need the basic information already when we initialze optimization to setup multiple excitations correctly */
    std::map<unsigned int, PtrParamNode> paramPerStep;

  private:
    
    //! Print out information about multisequence steps
    void WriteMultiSequenceStep(const UInt sequenceStep, const BasePDE::AnalysisType analysis);
    
    //! Set state to given step
    void SetupStep(UInt sequenceStep);

    //! Read restart information
    void ReadRestart();
    
    //! number of sequence steps
    UInt numSteps_;

    //! accumulated time
    Double accumulatedTime_;

    //! Flag, if analysis is restarted
    bool isRestarted_;
    
    //! current singleDriver object
    SingleDriver* actDriver_;
  
    /** shall we keep driver and pdes for optimization?
     * In that case the sequence step does not delete old stuff and reactivate existing stuff.
     * Triggered by the existence of an optimization element in xml
     * @see keptDrivers_  */
    bool keep_;

    /** For optimization we reuse the driver with keep_ set. */
    StdVector<SingleDriver*> keptDrivers_;

    /** For optimization we reuse the pdes from Domain with keep_ set.
     * Coudpled PDEs not yet supported */
    StdVector<StdVector<SinglePDE*> > keptPDEs_;


    //! stores for each step the participating pdes as name
    std::map<UInt, StdVector<std::string> > pdesPerStep_;

    //! stores for each step the analisystype of each pde
    std::map<UInt, BasePDE::AnalysisType > analysisPerStep_;
  };

}

#endif 
