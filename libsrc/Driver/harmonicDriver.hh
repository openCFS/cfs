#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField
{

//! driver for harmonic problems. it is derived from BaseDriver
class HarmonicDriver : public SingleDriver
{

public:
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  HarmonicDriver(Domain * adomain);

   //! deconstructor 
  virtual ~HarmonicDriver();

  //!  main method, where harmonic analysus is implemented.
  virtual void SolveProblem();


  private:
 
  Double  startFreq_;
  Double  stopFreq_;
  Integer numFreq_;
  Integer saveType_;
};

}

#endif // FILE_HARMONICDRIVER
