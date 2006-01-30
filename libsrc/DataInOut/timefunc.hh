#ifndef FILE_TIMEFUNC_2001
#define FILE_TIMEFUNC_2001

#include <list>
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"

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
    TimeFunc();

    //! Deconstructor
    ~TimeFunc();

    //! return value of time function with number <em>num</em> at the time <em>time</em>
    /*!
      \param time in: time
      \param fncname in: name of time function
    */
    Double TimeFuncAtTime(const Double time,  const std::string fncname);

    //! return values of time and space dependent function with number <em>num</em>
    //! at the time <em>time</em>.
    //! At the moment only one is allowed ('spc_dependent_fnc').
    //! That means that this function is called if 'spc_dependent_fnc' 
    //! is given as name of function data file in XML parameter file.
    /*!
      \param time in: time
      \param fncname in: name of time-space function ('spc_dependent_fnc')
      \param nodes in: vector containing nodes where IDBC is to be applied
      \param ptgrid in: pointer to grid to get coordinates of nodes
    */
    StdVector<Double> TimeSpcFuncAtTime(const Double time, 
                                               const std::string fncname,
                                               StdVector<UInt> nodes, Grid * ptgrid);
    

    //! Print values of time function in stream outfile
    void Print(std::ostream * outfile) const;

    //Set StartVector for structured Grid
    void SetStartTimeVector(Integer numBcs);

    //! return the number of time functions
    UInt GetmaxTimeFnc() 
    {
      return maxnumTF_;
    }

  private:

    //! read time functions from different dat-file
    void ReadTimeFuncs();
 
    //! 
    UInt maxnumTF_; //<! number of time functions

    //!
    std::list<Double> * timeTF_; //<! time 

    //!
    std::list<Double> * valTF_; //!< val

    //!
    StdVector<std::string> fnc_names_;

    //!
    Boolean timeFncDatFiles_;

    //!
    StdVector<Double> startTime_;
  };
} // end of namespace
#endif
