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
    MultiSequenceDriver(shared_ptr<SimState> state, Domain* domain,
                        PtrParamNode paramNode, PtrParamNode infoNode );

    //! destructir
    virtual ~MultiSequenceDriver();

    //! Initialization method
    void Init(bool restart);
    
    //! Change the current sequenceStep
    void SetSequenceStep(UInt sequenceStep);
  
    //! main method, where time-stepping is implemented. 
    //! it is for transient and static problem
    void SolveProblem(bool write_results);

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
    SingleDriver * actDriver_;
  
    //! stores for each step the participating pdes as name
    std::map<UInt, StdVector<std::string> > pdesPerStep_;

    //! stores for each step the analisystype of each pde
    std::map<UInt, BasePDE::AnalysisType > analysisPerStep_;

  };

}

#endif 
