// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
/*!
 *       \file     CoefFunctionGridNodalDefault.hh
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *
 *       \date     May 2015
 *       \author   Manfred Kaltenbbacher
 */
//================================================================================================
#ifndef FILE_INVERSESOURCEDRIVER_2014
#define FILE_INVERSESOURCEDRIVER_2014

#include "SingleDriver.hh"

#include <memory>


namespace CoupledField
{

class ParamNode;
class Timer;

//! driver for identifying RHS of a PDE. it is derived from BaseDriver
class InverseSourceDriver : public virtual SingleDriver
{

public:

  //! Constructor
  //! \param sequenceStep current step in multisequence simulation
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  InverseSourceDriver(UInt sequenceStep, bool isPartOfSequence,
                 shared_ptr<SimState> state,
                 Domain* domain, PtrParamNode paramNode, PtrParamNode infoNode );

  //! Detructor
  virtual ~InverseSourceDriver();

  //! Return current time / frequency step of simulation
  UInt GetActStep(const std::string& pdename)
  {
    return actFreqStep_;
  }

  //! Return current frequency
  Double GetActFreq()
  {
    return freq_;
  }

  //! Initialization method
  void Init(bool restart);

  //! Main method, where harmonic analysis is implemented.
  void SolveProblem();

  /** This allows optimization to handle the individual frequency steps, e.g. to compute
   * objective values. Internally this is is a service function for SolveProblem()
   * @param actFreqStep sets the actFreq_ attribute, to start with 1 and not to exceed numFreq_ */
   Double ComputeFrequencyStep(UInt actFreqStep);

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
  };

  StdVector<Frequency> freqs;

  /** Helper method which determines if an AnalyisType is complex. */
  virtual bool IsComplex() { return true; };

  //! Static method being called in the case of a Ctr-C signal
  static void SignalHandler( int sig);

  
protected:

  //! Frequency step to proceed from when performing restarted simulation
  UInt restartStep_; 
  
  //! Current frequency value
  Double actFreq_;

  //! First frequency for which a simulation is performed
  Double freq_;

  //! Number of frequencies for which a simulation is performed
  UInt numFreq_;

  //! Current frequency step
  UInt actFreqStep_;

  //! Last frequency step
  UInt stopFreqStep_;
  
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

  //! pointer to loadRHS: source
  PtrCoefFct rhsSource_;

  //! pointer to loadRHS: measure
  PtrCoefFct rhsMeas_;

  //!file name containing the measured data
  std::string fileNameMeasdata_;

  //! regularization parameter 1
  Double alpha_;

  //! regularization parameter 2
  Double beta_;

  //! regularization parameter 3
  Double rho_;

  //! regularization parameter 4
  Double qExp_;

  //! relative stopping criterion for outer iteration
  Double resStopCritRel_;

  //! maximal number of iterations reducing alpha, beta, rho (outer loop)
  UInt maxReduceParSteps_;

  //! maximal number of gradient steps (inner loop)
  UInt maxGradSteps_;

  //! maximal number of iterations
  UInt maxNumLineSearch_;

  //!square of L2-norm of measured pressure at mic-positions
  Double measL2squared_;

  //! minimal allowed distance of two micophones
  Double minMicDistant_;

  //! level of log writing
  std::string logLevel_;

  //! scales the measured data to this value
  Double scale2Val_;
};

}

#endif // FILE_INVERSESOURCERIVER
