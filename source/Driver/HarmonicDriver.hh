#ifndef FILE_HARMONICDRIVER_2001
#define FILE_HARMONICDRIVER_2001

#include "SingleDriver.hh"

#include <memory>


namespace CoupledField
{


class Timer;

//! driver for harmonic problems. it is derived from BaseDriver
class HarmonicDriver : public virtual SingleDriver
{

public:

  //! Constructor
  //! \param sequenceStep current step in multisequence simulation
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  HarmonicDriver(UInt sequenceStep, bool isPartOfSequence,
                 shared_ptr<SimState> state,
                 Domain* domain, PtrParamNode paramNode, PtrParamNode infoNode );

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

  /** helper for ContextManager in case of multiple sequence optimization.  */
  static unsigned int GetNumFreq(PtrParamNode node);

  //! Initialization method
  void Init(bool restart);

  //! Main method, where harmonic analysis is implemented.
  void SolveProblem();

   /** This StoreResults meant for Optimization only */
   unsigned int StoreResults(UInt stepNum, double step_val);

  //! \copydoc SingleDriver::SetToStepValue
  virtual void SetToStepValue(UInt stepNum, Double stepVal );

  /** This is the list of all frequencies. As long as we have no adaptive
   * frequeceny steps this makes no problem. */
  struct Frequency
  {
    /** the frequencys are 1-based such freqs[i].step should be i+1 */
    unsigned int step;

    double freq;

    /** The weight is only used in multiple frequency optimization */
    double weight;

    friend std::ostream &
    operator<<(std::ostream & os, Frequency const & freqStp)
    {
       return os << "Frequency: " << freqStp.freq << ", Step: " 
        << freqStp.step << ", Weight: " << freqStp.weight;
    }
  };

  StdVector<Frequency> freqs;

  /** This allows optimization to handle the individual frequency steps, e.g. to compute
   * objective values. Internally this is is a service function for SolveProblem()
   * @param freq sets the actFreq_ and actFreqStep attributes, to start with 1 and not to exceed numFreq_ */
   Double ComputeFrequencyStep(Frequency const& freqStp);

  /** Helper method which determines if an AnalyisType is complex. */
  virtual bool IsComplex() { return true; };

  //! Static method being called in the case of a Ctr-C signal
  static void SignalHandler( int sig);

  
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

  //! Frequency step to proceed from when performing restarted simulation
  UInt restartStep_; 
  
  //! Current frequency value
  Double actFreq_;

  //! First frequency for which a simulation is performed
  Double startFreq_;

  //! Last frequency for which a simulation is performed
  Double stopFreq_;
  
  //! Last frequency step
  UInt stopFreqStep_;

  //! Number of frequencies for which a simulation is performed
  UInt numFreq_;

  //! Current frequency step
  UInt actFreqStep_;

  //! A attribute storing the type of algorithm used for frequency sampling
  FreqSamplingType samplingType_;
  
  //! Flag whether stepping strategy is standard or tree based
  bool treeStepping_;

  //! Number of tree levels in case of tree stepping
  UInt maxTreeLevel_;

  //! Number of frequency steps, that should be included in the output
  UInt numSteps_;

  // =======================================================================
  //  Restart related data
  // =======================================================================
  
  //! Read restart information
  void ReadRestart();
  
  //! Flag, if analysis is restarted
  bool isRestarted_;
  
  //! Flag, if restart file is to be written
  bool writeRestart_;
  
  // =======================================================================
  //  Timing estimation
  // =======================================================================

  //! Estimated time per step
  Double timePerStep_;
  
  //! Timer for estimating remaining runtime 
  shared_ptr<Timer> timer_;

};

}

#endif // FILE_HARMONICDRIVER
