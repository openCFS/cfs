#ifndef FILE_TIMEFUNC_2001
#define FILE_TIMEFUNC_2001

#include "filetype.hh"

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

    //! return value of time function with number {\tt num} at the time {\tt time}.
    /*!
      \param time in: time
      \param num in: number of the time function
    */
    Double TimeFuncAtTime(const Double time, const Integer num);

    //! Print values of time function in stream outfile
    void Print(std::ostream * outfile) const;

    //! return the number of time functions
    Integer GetmaxTimeFnc() 
    {
      return maxnumTimeFunc_;
    }

private:

    //! pointer to input file. needed only for {\tt datfile}.
    FileType * ptFileType;

    //! read time func from dat-file
    void ReadTimeFunc(const std::string nametf);

    //! read time func old style (needs header line with func nr & nr of entries)
    void ReadTimeFuncOldType(const std::string nametf);
    
    //! value of time func through interpolation values from data-file
    Double ValTimeFuncDatFile(const Double time, const Integer num);

    //! 
    Integer maxnumTimeFunc_;

    //!
    Integer * maxvalTimeFunc;

    //!
    Double ** timeTimeFunc;

    //!
    Double ** valTimeFunc;

    //!
    Boolean timeFncDatFile_;

    //!
    Integer argTimeFnc;

    //!
    pfn1var ptTimeFnc;

    //! interval of time func 
    Double intervalTF_a;
    Double intervalTF_b;

};
} // end of namespace
#endif
