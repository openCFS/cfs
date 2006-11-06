#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField
{

  //! driver for harmonic problems. it is derived from BaseDriver
  class HarmonicDriver : public virtual SingleDriver {

  public:

    //! Constructor
    //! \param stepOffset offset for starting (frequency)step
    //! \param reqOffset offset for starting frequency
    //! \param driverTag tag for current driver
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    HarmonicDriver( UInt stepOffset = 0,
                    Double freqOffset = 0.0,
                    std::string driverTag = "anyTag",
                    bool isPartOfSequence = false);

    //! Detructor 
    virtual ~HarmonicDriver();

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { 
      return actFreqStep_; }

    //! Return current frequency
    Double GetActFreq( ) { return actFreq_; }

    //! Initialization method
    void Init();
    
    //! Main method, where harmonic analysis is implemented.
    void SolveProblem();

  protected:

    //! Method for computing the actual frequency for a given frequency number

    //! This method implements the different sampling variants for the
    //! frequency domain. We currently support linear, logarithmic and
    //! reverse logarithmic sampling.
    //! \param freqIndex index of the frequency that shall be computed, i.e.
    //!                  an integral value from [1:numFreq_]
    Double ComputeNextFrequency( UInt freqIndex );

    //! Current frequency value
    Double actFreq_;
    
    //! First frequency for which a simulation is performed
    Double startFreq_;

    //! Last frequency for which a simulation is performed
    Double stopFreq_;

    //! Number of frequencies for which a simulation is performed
    UInt numFreq_;

    //! Current frequency step
    UInt actFreqStep_;

    //! A attribute storing the type of algorithm used for frequency sampling
    FreqSamplingType samplingType_;

    //! Offset for first frequency step
    UInt stepOffset_;

    //! offset for first time
    Double freqOffset_;

    //! Flag for signalling frequency dependend damping

    //! If we have frequency dependend damping, then for each new frequency
    //! the damping must be adapted2
    bool adjustDamping_;
  };

}

#endif // FILE_HARMONICDRIVER
