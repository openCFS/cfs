#ifndef FILE_FREQFUNC_2001
#define FILE_FREQFUNC_2001

#include <list>
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"

namespace CoupledField
{

  //! Class interpolates and retrieves information from given nodal freq functions
  /*!
    This class interpolates and retrieves information from given nodal freq functions. 
  */

  class FreqFunc
  {
  public:
    //! Constructor
    FreqFunc();

    //! Deconstructor
    ~FreqFunc();

    //! return ampl. and phase of nodal freq function with number <em>num</em> at the freq <em>freq</em>
    /*!
      \param freq in: freq
      \param fncname in: name of freq function
      \param nodeNumber in: node from which its ampl. and phase are to be retrieved from its nodal freq. file
    */
    StdVector<Double> NodalFreqFuncAtFreq(const Double freq,  const std::string fncname, UInt nodeNumber);

    //! return the number of freq functions
    UInt GetmaxFreqFnc() 
    {
      return maxnumFF_;
    }

  private:

    //! read a nodal freq function file from a set files under ./nodalSrcs/
    //! directory having common root names and <em>nodeNumber</em> as suffix
    void ReadFreqFuncs(UInt nodeNumber);
 
    //! 
    UInt maxnumFF_; //<! number of freq functions

    //!
    std::list<Double> * freqFF_; //<! freq 

    //!
    std::list<Double> * valFF_; //!< ampplitude

    //!
    std::list<Double> * phaseFF_; //!< phase

    //!
    StdVector<std::string> fnc_names_;

    //!
    StdVector<Double> startFreq_;
  };
} // end of namespace
#endif
