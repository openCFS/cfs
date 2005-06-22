#ifndef PROFILER_HH
#define PROFILER_HH

#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h> 



namespace CoupledField {

  class Profiler {

  public:
    
    //!
    Profiler();

    virtual ~Profiler();

    //!
    void Trace(std::string name);
    
  private:
    
    //! process id of this program
    pid_t myPid_;
    
    //! ostream for tracing
    std::ofstream *memOut_;
  };


} // end of namespace


#endif
