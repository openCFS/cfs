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
  //! \param adomain pointer to class Domain
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time
  //! \param driverTag tag for current driver
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  HarmonicDriver(Domain * adomain,
		 Integer stepOffset = 0,
		 Double timeOffset = 0.0,
		 std::string driverTag = "anyTag",
		 Boolean isPartOfSequence = FALSE);

   //! deconstructor 
  virtual ~HarmonicDriver();

  //!  main method, where harmonic analysus is implemented.
  virtual void SolveProblem();


  private:
 
  Double  startFreq_;
  Double  stopFreq_;
  Integer numFreq_;
  Integer saveType_;
  Boolean adjustDamping_;
};

}

#endif // FILE_HARMONICDRIVER
