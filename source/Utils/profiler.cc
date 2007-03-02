#include "profiler.hh"

#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <iomanip>

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "myclock.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  Profiler::Profiler() {

    // Determine whether profiling should be performed
    doProfiling_ = commandLine->GetDoProfile() == true;

    if ( doProfiling_ == true ) {

      // get process id
      myPid_ = getpid();

      // intialize profile file
      fileName_ = commandLine->GetSimName() + ".profile";
      memOut_ = new std::ofstream(fileName_.c_str());

      (*memOut_).setf(std::ios::fixed,std::ios::floatfield);
      (*memOut_).precision(2);

      (*memOut_) << "# Profiling output\n"
                 << "# Memory | Wall time | User Time "
                 << "| Rel. wall time | Rel. user time | Name\n";
      (*memOut_) << std::setw(80) << std::setfill('=') << "" 
                 << std::setfill(' ') << std::endl;

      // Initialise attributes
      wTimeLast_ = 0.0;
      cTimeLast_ = 0.0;

      // Start timing
      clock_ = new MyClock();
      clock_->Start();
    }

    // We do not do profiling
    else {
      memOut_ = NULL;
      clock_  = NULL;
    }
  }


  // **************
  //   Destructor
  // **************
  Profiler::~Profiler() {

    delete memOut_;
    memOut_ = NULL;

    delete clock_;
    clock_ = NULL;
  }
  
  
  // *********
  //   Trace
  // *********
  void Profiler::Trace( std::string name ) {

    if ( doProfiling_ == true ) {

      std::ostringstream command;
      Double wTime, cTime;

      // Get memory consumption
      command << "ps -F " << myPid_ 
              << " | grep " << myPid_ 
              << " | gawk '{printf \"% 5d MB\",$6/1024}' >> " << fileName_;
      memOut_->close();
      std::system( (command.str()).c_str() );
      memOut_->open( fileName_.c_str(), std::ofstream::out |
                     std::ofstream::app);

      // Get time
      clock_->GetTime(wTime, cTime);
      (*memOut_) << "\t" << wTime << "\t" << cTime;

      // Write relative time
      (*memOut_) << "\t" <<  wTime-wTimeLast_
                 << "\t" <<  cTime-cTimeLast_;
      wTimeLast_ = wTime;
      cTimeLast_ = cTime;

      // Write name
      (*memOut_)  << '\t' << name << std::endl;
    }
  }

}
