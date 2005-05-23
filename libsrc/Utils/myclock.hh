#ifndef CLOCK_2001
#define CLOCK_2001

#include <string>
#include <fstream>
#include "General/environment.hh"


namespace CoupledField 
{

  //! Count wall and CPU time
  /*!
    This class provides tools for timing program
  */
  class MyClock
  {
  public:

    //! Constructor
    /*!
      \param atitle name for the output. it can be omitted.
    */
    MyClock(char * atitle=NULL);

    //! enumeration for the status of a clock. 
    /*!
      \param beg start to count the time
      \param finish end to count the time
    */
    enum status{beg, end};
 
    //! Function for measure time of running programm
    /*!
      \param status beg or end
      \param astring name of the process, that is timed. it can be omitted.
    */
    void ClockCount(enum status, const std::string astring="");
 
    //! Deconstructor(close file for output time).
    ~MyClock();
 
  private:
    //!
    time_t tm, tm_tmp;
    //!
    clock_t ck;
    //!
    clock_t ck_tmp;
    //!
    std::ofstream filetime;
    //!
    Boolean InFile;
 
  };

} // end of namespace
#endif
