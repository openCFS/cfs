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
                  bool isPartOfSequence,
                  shared_ptr<SimState> state, Domain* domain,
                  PtrParamNode paramNode, PtrParamNode infoNode );

    //! Destructor 
    ~StaticDriver();
  
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

    //! Initialization method
    void Init(bool restart);

    /** @see BaseDriver::SolveProblem(double) */  
    void SolveProblem(bool write_results = true);
        
    /** @see BaseDriver::StoreResults(double) */  
    void StoreResults(UInt stepNum, double step_val);

    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() { return false; };

  protected:

    //! \copydoc SingleDriver::SetToStepValue
    virtual void SetToStepValue(UInt stepNum, Double stepVal );
  };

}

#endif // FILE_STATICDRIVER
