#ifndef FILE_TIMEFUNC_2001
#define FILE_TIMEFUNC_2001

#include "filetype.hh"

namespace CoupledField
{

  //! Class we store information about given time function
  /*!
   This class stores information about given time function. Later it should provide methods for handling time functions, interpolation (1D, 2D, 3D) of initial and
load (nodal-, surface-, volume-type) data.
   */

class TimeFunc
  {
public:
    //! Constructor
    TimeFunc(FileType *);

    //! Deconstructor
    ~TimeFunc();

    //! Print 
    void Print(std::ostream * outfileDat) const;

    //! 
    Double TimeFuncAtTime(const Double, const Integer num);

    //!
    //    Double TimeFuncWaveSt(const Double time);

private:

    //!
    FileType * ptFileType;

    //! read time func from dat-file
    void ReadTimeFunc(const std::string nametf);

    //! value of time func through interpolation values from data-file
    Double ValTimeFuncDatFile(const Double time, const Integer num);

    //!
    Integer maxnumTimeFunc;

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
