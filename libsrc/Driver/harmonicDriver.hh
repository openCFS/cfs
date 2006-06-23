#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField
{

  //! driver for harmonic problems. it is derived from BaseDriver
  class HarmonicDriver : public SingleDriver {

  public:

    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    HarmonicDriver(Domain * adomain,
                   UInt stepOffset = 0,
                   Double timeOffset = 0.0,
                   std::string driverTag = "anyTag",
                   bool isPartOfSequence = false);

    //! Detructor 
    virtual ~HarmonicDriver();

    //! Main method, where harmonic analysis is implemented.
    virtual void SolveProblem();

  private:

    //! Method for computing the actual frequency for a given frequency number

    //! This method implements the different sampling variants for the
    //! frequency domain. We currently support linear, logarithmic and
    //! reverse logarithmic sampling.
    //! \param freqIndex index of the frequency that shall be computed, i.e.
    //!                  an integral value from [1:numFreq_]
    Double ComputeNextFrequency( UInt freqIndex );

    //! First frequency for which a simulation is performed
    Double startFreq_;

    //! Last frequency for which a simulation is performed
    Double stopFreq_;

    //! Number of frequencies for which a simulation is performed
    UInt numFreq_;

    //! A attribute storing the type of algorithm used for frequency sampling
    FreqSamplingType samplingType_;

    //! Flag for signalling frequency dependend damping

    //! If we have frequency dependend damping, then for each new frequency
    //! the damping must be adapted2
    bool adjustDamping_;
  };

}

#endif // FILE_HARMONICDRIVER
