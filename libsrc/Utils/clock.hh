#ifndef CLOCK_2001
#define CLOCK_2001

using namespace std;
namespace CoupledField 
{

  //! Count wall and CPU time
  /*!
    This class provides tools for timing program
  */

class Clock
{
public:

  //! Constructor
  Clock(char * atitle=NULL);

  //!
  enum status{beg, end};
 
  //! Function for measure time of running programm
   void ClockCount(enum status, const string astring="");
 
  //! Deconstructor(close file for output time).
   ~Clock();
 
private:
  time_t tm, tm_tmp;
  clock_t ck;
  clock_t ck_tmp;
  ofstream filetime;
  Boolean InFile;
 
};

} // end of namespace
#endif
