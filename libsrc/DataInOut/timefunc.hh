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
    TimeFunc(const Char * const name);

    //! Deconstructor
    ~TimeFunc();

    //! Print 
    void Print(std::ostream * outfileDat) const;

    //! 
    Double TimeFuncAtTime(const Double, const Integer num);

    //!
    Double TimeFuncWaveSt(const Double time);

private:

    //!
    std::ifstream timefncfile;

    //!
    Integer maxnumTimeFunc;

    //!
    Integer * maxvalTimeFunc;

    //!
    Double ** timeTimeFunc;

    //!
    Double ** valTimeFunc;

    //!
    ofstream testtf;
};
} // end of namespace
#endif
