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
    virtual void SolveProblem() = 0;
  
  protected:
  
    //! Trigger reading of restart
    virtual void ReadRestart( ) {};
    
    //! Initialize PDEs
    void InitializePDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    //! true, if driver is part of a multiSequence
    bool isPartOfSequence_;

    //! current sequences step in multiSequence simulation
    UInt sequenceStep_;
    
  };

}

#endif 
