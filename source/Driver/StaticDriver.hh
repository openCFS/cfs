#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "SingleDriver.hh"

namespace CoupledField {

  //! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public SingleDriver {

  public:
    //! constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param true, if driver is part of  multiSequence
    StaticDriver( UInt sequenceStep,
                  bool isPartOfSequence = false );

    //! Destructor 
    ~StaticDriver();
  
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

    //! Initialization method
    void Init();

    /** @see BaseDriver::SolveProblem(double) */  
    void SolveProblem(bool write_results = true, PtrParamNode analysis_id = PtrParamNode(), AdjointParameters* adjointParams = NULL);
        
    /** @see BaseDriver::StoreResults(double) */  
    void StoreResults(UInt stepNum, double step_val);

    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() { return false; };

  protected:
    //! Number of steps before a restart file is stored
    UInt restartIncr_;

  };

}

#endif // FILE_STATICDRIVER
