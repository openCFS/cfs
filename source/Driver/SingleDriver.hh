#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "BaseDriver.hh"

namespace CoupledField {

  // forward class declaration
  class BasePDE;
  class MathParser;
  
  //! Abstract base class for single driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver {
    
  public:
    
    //! Constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    SingleDriver( UInt sequenceStep, bool isPartOfSequence,
                  shared_ptr<SimState> state, Domain* domain,
                  PtrParamNode paramNode, PtrParamNode infoNode );
    
    //! Default destructor
    virtual ~SingleDriver() {};
    
    //! set the pdes, which gets to be solved
    void SetPDE(BasePDE * pde);

    /** implement abstract identification class */ 
    DriverClass GetDriverClass() { return SINGLE_DRIVER; };

    //! Set step value from outside
    
    //! This methods allows to set the time / frequency step from outside.
    //! It is only allowed, if the SimSate object provides additional input.
    virtual void SetToStepValue(UInt stepNum, Double stepVal ) {
      EXCEPTION( "Not implemented" );
    }

    bool IsPartOfSequence() {return isPartOfSequence_;}

  protected:
  
    //! Initialize PDEs
    void InitializePDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    /** Flag if internal state of PDE (=FeFunctions) get written each step.
     * Is false and only true if set via allowPostProc in the xml file is set */
    bool writeAllSteps_;
    
    //! true, if driver is part of a multiSequence, false if first run or single run 
    bool isPartOfSequence_;

    //! Pointer to math parser instance
    MathParser * mathParser_;
    
    //! approximate source with delta functions
    bool approxSourceWithDeltaFnc_;

  };

}

#endif 
