#include "myclock.hh"

namespace CoupledField
{

  MyClock::MyClock () {
    ENTER_FCN( "MyClock::MyClock" );
  }
 
  MyClock::~MyClock() {
    ENTER_FCN( "MyClock::~MyClock" );
  }

  void MyClock::Start() {
    ENTER_FCN( "MyClock::Start" );

    wTimeStart_ = std::time(NULL);
    cTimeStart_ = std::clock();
  }

  void MyClock::Reset() {
    ENTER_FCN( "MyClock::Reset" );
    wTimeStart_ = std::time(NULL);
    cTimeStart_ = std::clock();
  }

  void MyClock::GetTime(Double & wallTime, Double & userTime) {
    
    std::time_t wTimeCur, cTimeCur;

    // get current time and clock counts
    wTimeCur = std::time(NULL);
    cTimeCur = std::clock();
    
    wallTime = (Double) std::difftime(wTimeCur, wTimeStart_);
    userTime = ( (Double) (cTimeCur - cTimeStart_)) / CLOCKS_PER_SEC;
  }
  
} // end of namespace  
