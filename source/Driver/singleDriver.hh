// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "basedriver.hh"

namespace CoupledField {

  // forward class declaration
  class BasePDE;

  //! Abstract base class for sinlge driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver {
    
  public:
    
    //! Constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    SingleDriver( UInt sequenceStep,
                  bool isPartOfSequence = false );
    
    //! Default destructor
    virtual ~SingleDriver();
    
    //! set the pdes, which gets to be solved
    void SetPDE(BasePDE * pde);

    //! main method, where time-stepping is implemented. it is for transient and static problem
    //! There is an internal state consecutiveRun/isPartOfSequence which is controlled in the
    //! optimization case via SetConsecutiveRun(). No internal state switch!
    virtual void SolveProblem();

    /** allows the manipulation of the driver "state" */
    void SetConsecutiveRun(bool value) {
        consecutiveRun_ = value;
    };

    /** implement abstract identification class */ 
    DriverClass GetDriverClass() { return SINGLE_DRIVER; };

  protected:
  
    //! Trigger reading of restart
    virtual void ReadRestart( ) {};
    
    //! Initialize PDEs
    void InitializePDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    //! true, if driver is part of a multiSequence, false if first run or single run 
    bool isPartOfSequence_;
    
    //! Also used for tagging runs after the first run in optimization case. See SolveStep() 
    bool consecutiveRun_;

    //! current sequences step in multiSequence simulation
    UInt sequenceStep_;
    
  };

}

#endif 
