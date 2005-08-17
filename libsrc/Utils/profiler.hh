#ifndef PROFILER_HH
#define PROFILER_HH

#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h> 

#include "General/environment.hh"


namespace CoupledField {

  // forward class declarations
  class MyClock;

  //! Class for gathering information about memory consumption and time
  class Profiler {

  public:
    
    //! Constructor
    Profiler();

    //! Destructor
    virtual ~Profiler();

    //! Write tracing marker into file
    void Trace(std::string name);
    
  private:
    
    //! Process id of this program
    pid_t myPid_;

    //! Clock object
    MyClock * clock_;
    
    //! Name of output file
    std::string fileName_;
    
    //! File stream for output
    std::ofstream *memOut_;

    //! Previous timing values
    double wTimeLast_, cTimeLast_;

    //! Flag for enabling/disabling profiling

    //! In this flag we store the information whether profiling is to be
    //! performed, or whether the object should only perform a dummy
    //! behaviour, i.e. we can call it, but the methods simply do nothing.
    Boolean doProfiling_;

  };


} // end of namespace


#endif
