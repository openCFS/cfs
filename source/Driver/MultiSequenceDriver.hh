#ifndef FILE_MULTISEQUENCEDRIVER_2004
#define FILE_MULTISEQUENCEDRIVER_2004

#include "BaseDriver.hh"
#include "SingleDriver.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  // forward class declaration
  class BasePDE;

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
  
    //! constructor
    MultiSequenceDriver(shared_ptr<SimState> state, Domain* domain );

    //! destructir
    virtual ~MultiSequenceDriver();

    //! Initialization method
    void Init();
  
    //! main method, where time-stepping is implemented. 
    //! it is for transient and static problem
    void SolveProblem(bool write_results = true, PtrParamNode given_analysis_id = PtrParamNode());

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename );

    //! Return current singleDriver 
    SingleDriver* GetSingleDriver() { return actDriver_; }

    /** implement abstract identification class */ 
    DriverClass GetDriverClass() { return MULTI_SEQUENCE_DRIVER; };

    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() { return actDriver_->IsComplex(); };

  private:
    
    //! Print out information about multisequence steps
    void WriteMultiSequenceStep(const UInt sequenceStep, 
                                const BasePDE::AnalysisType analysis);

    //! number of sequence steps
    UInt numSteps_;

    //! current sequence step
    UInt curSequenceStep_;

    //! accumulated time
    Double accumulatedTime_;

    //! current singleDriver object
    SingleDriver * actDriver_;
  
    //! stores for each step the participating pdes as name
    StdVector<StdVector<std::string> > pdesPerStep_;

    //! stores for each step the analisystype of each pde
    StdVector<BasePDE::AnalysisType > analysisPerStep_;

  };

}

#endif 
