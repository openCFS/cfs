// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "singleDriver.hh"

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
    void SolveProblem(bool write_results = true, const std::string& comment = "");
        
    /** @see BaseDriver::StoreResults(double) */  
    void StoreResults(double step_val);

  protected:
    //! Number of steps before a restart file is stored
    UInt restartIncr_;

  };

}

#endif // FILE_STATICDRIVER
