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

private:

    //!
    FileType * ptFileType;

    //!
    Integer maxnumTimeFunc;

    //!
    Integer * maxvalTimeFunc;

    //!
    Double ** timeTimeFunc;

    //!
    Double ** valTimeFunc;
};
} // end of namespace
#endif
