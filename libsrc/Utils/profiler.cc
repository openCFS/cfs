#include "profiler.hh"

#include <iostream>
#include <sstream>

namespace CoupledField{


  Profiler::Profiler() {
    myPid_ = getpid();
    memOut_ = new std::ofstream("cfs.mem");
  }
  Profiler::~Profiler() {
    if ( memOut_ )
      delete memOut_;
  }
  
  
  void Profiler::Trace(std::string name) {
    std::ostringstream command;
    command << "ps -F " << myPid_ 
            << " | grep " << myPid_ 
            << " | gawk '{printf $6}' >> cfs.mem";
    memOut_->close();
    std::system( (command.str()).c_str() );
    memOut_->open("cfs.mem", std::ofstream::out | std::ofstream::app);
    (*memOut_)  << '\t' << name << std::endl;

    
  }
  
  
} // end of namespace
