#ifndef CLOCK_2001
#define CLOCK_2001


#include <string>
#include <ctime>
#include "General/environment.hh"


namespace CoupledField 
{

  //! Count wall and CPU time
  class MyClock {
    
  public:
    
    //! Constructor
    MyClock();

    //! Destructor
    virtual ~MyClock();

    //! Start clock counter
    void Start();

    //! End clock counter
    void Reset();

    //! Get time since Start() command
    void GetTime(Double & wallTime, Double & userTime );

  private:

    //! Start value of wall time
    std::time_t wTimeStart_;

    //! Start value of cpu / clock time
    std::clock_t cTimeStart_;
 
  };

} // end of namespace
#endif
