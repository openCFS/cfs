// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField
{

class ParamNode;
class Timer;

//! driver for harmonic problems. it is derived from BaseDriver
class HarmonicDriver : public virtual SingleDriver
{

public:

  //! Constructor
  //! \param sequenceStep current step in multisequence simulation
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  HarmonicDriver(UInt sequenceStep, bool isPartOfSequence = false);

  //! Detructor
  virtual ~HarmonicDriver();

  //! Return current time / frequency step of simulation
  UInt GetActStep(const std::string& pdename)
  {
    return actFreqStep_;
  }

  //! Return current frequency
  Double GetActFreq()
  {
    return actFreq_;
  }

  //! Return start frequency
  Double GetStartFreq()
  {
    return startFreq_;
  }

  //! Return stop frequency
  Double GetStopFreq()
  {
    return stopFreq_;
  }

  //! Initialization method
  void Init();

  //! Main method, where harmonic analysis is implemented.
  void SolveProblem(bool write_results = true, PtrParamNode analysis_id = PtrParamNode(), AdjointParameters* adjointParams = NULL);

  /** This allows optimization to handle the individual frequency steps, e.g. to compute
   * objective values. Internally this is is a service function for SolveProblem()
   * @param actFreqStep sets the actFreq_ attribute, to start with 1 and not to exceed numFreq_ */
   Double ComputeFrequencyStep(UInt actFreqStep, PtrParamNode analysis_id);

   /** This StoreResults meant for Optimization only */
  void StoreResults(UInt stepNum,
                    double step_val );

  /** This is the list of all frequencies. As long as we have no adaptive
   * frequeceny steps this makes no problem. */
  struct Frequency
  {
    /** the frequencys are 1-based such freqs[i].step should be i+1 */
    unsigned int step;

    double freq;

    /** The weight is only used in multiple frequency optimization */
    double weight;
  };

  StdVector<Frequency> freqs;

protected:


  /** parse the frequencyList if it exists
   * @return if the frquencyList exists*/
  bool ReadFrequencyList();

  /** read the numFreq/startFreq/stopFreq */
  bool ReadParametrizedFrequencies();

  //! Method for computing the actual frequency for a given frequency number

  //! This method implements the different sampling variants for the
  //! frequency domain. We currently support linear, logarithmic and
  //! reverse logarithmic sampling.
  //! \param freqIndex index of the frequency that shall be computed, i.e.
  //!                  an integral value from [1:numFreq_]
  Double ComputeNextFrequency(UInt freqIndex) const;

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

  /** This is the pointer to our analysis description */
  PtrParamNode pn_;
  
  // =======================================================================
  //  Timing estimation
  // =======================================================================

  //! Timer for estimating remaining runtime 
  Timer * timer_;

};

}

#endif // FILE_HARMONICDRIVER
