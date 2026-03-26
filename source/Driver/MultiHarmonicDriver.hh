// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MultiHarmonicDriver.hh
 *       \brief    Driver for nonlinear frequency domain analysis
 *
 *       \date     Jan 2, 2017
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_MULTIHARMONICDRIVER
#define FILE_MULTIHARMONICDRIVER

#include "SingleDriver.hh"

#include <memory>


namespace CoupledField
{


class Timer;

//! driver for multiharmonic problems (nonlinear in frequency-domain). it is derived from BaseDriver
class MultiHarmonicDriver : public virtual SingleDriver
{

public:

  //! Constructor
  //! \param sequenceStep current step in multisequence simulation
  //! \param isPartOfSequence true, if driver is part of  multiSequence
  MultiHarmonicDriver(UInt sequenceStep, bool isPartOfSequence,
                 shared_ptr<SimState> state,
                 Domain* domain, PtrParamNode paramNode, PtrParamNode infoNode );


  //! Destructor
  virtual ~MultiHarmonicDriver();

  //! Return current time / frequency step of simulation
  UInt GetActStep( const std::string& pdename ) { return 1;}

  //! Initialization method
  void Init(bool restart);

  //! Main method, where harmonic analysis is implemented.
  void SolveProblem();

  //! \copydoc SingleDriver::SetToStepValue
  virtual void SetToStepValue(UInt stepNum, Double stepVal );

  /** Helper method which determines if an AnalyisType is complex. */
  virtual bool IsComplex() { return true; };

  //! Compute the index of a given harmonic
  UInt IndexOfHarmonic(const Integer& harmonic);

  //! Compute the harmonic of a given harmonic
  Integer HarmonicOfIndex(const UInt& Index);

  //! Get the number of the considered frequencies (positive and negative)
  UInt GetNumFreq(){return numFreq_;}

  //! True, if we are considering the full system with all harmonics (odd and even ones)
  bool IsFullSystem(){ return fullSystem_; }

  //! Base frequency for which a simulation is performed
  Double baseFreq_;
  
  //! Number of considered harmonics for the solution quantity
  UInt numHarmonics_N_;
  
  //! Number of considered harmonics for the nonlinearity, e.g. nonlinear BH-curve in electromagnetics
  UInt numHarmonics_M_;

  //! Number of considered time evaluation points for FFT and iFFT
  UInt numFFT_;

  //! Harmonic, we are currently considering
  UInt actHarm_;

  //! Boolean, which tells us if we need to incorporate the zero harmonic
  bool fullSystem_;

  //! Vector containing the harmonic-frequencies (negative, 0, positive)
  StdVector<Double> harmFreq_;

  //! Number of considered frequencies (positive and negative)
  UInt numFreq_;

protected:

  //! Timer for estimating remaining runtime
  shared_ptr<Timer> timer_;
};

}

#endif // FILE_MULTIHARMONICDRIVER
