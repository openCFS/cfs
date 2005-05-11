#ifndef FILE_TIMEFUNC_2001
#define FILE_TIMEFUNC_2001

#include "filetype.hh"
#include <list>

namespace CoupledField
{

  //! Class we store information about given time function
  /*!
   This class stores information about given time function. 
   */

class TimeFunc
  {
public:
    //! Constructor
    TimeFunc(FileType *);

    //! Deconstructor
    ~TimeFunc();

    //! return value of time function with number <em>num</em> at the time <em>time</em>
    /*!
      \param time in: time
      \param fncname in: name of time function
    */
    Double TimeFuncAtTime(const Double time,  const std::string fncname);

    //! Print values of time function in stream outfile
    void Print(std::ostream * outfile) const;

    //! return the number of time functions
    Integer GetmaxTimeFnc() 
    {
      return maxnumTF_;
    }

private:

    //! pointer to input file. needed only for <em>datfile</em>
    FileType * ptFileType;

    //! read time functions from different dat-file
    void ReadTimeFuncs();
 
    //! 
    Integer maxnumTF_; //<! number of time functions

    //!
    std::list<Double> * timeTF_; //<! time 

    //!
    std::list<Double> * valTF_; //!< val

    //!
    StdVector<std::string> fnc_names_;

    //!
    Boolean timeFncDatFiles_;
};
} // end of namespace
#endif
