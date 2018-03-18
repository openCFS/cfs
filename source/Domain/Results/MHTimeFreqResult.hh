// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MHTimeFreqResult.hh
 *       \brief    Storage container for multiharmonic results
 *
 *       \date     Mar 9, 2018
 *       \author   kroppert
 */
//================================================================================================

#ifndef FILE_CFS_MHTIMEFREQRESULT_HH
#define FILE_CFS_MHTIMEFREQRESULT_HH

#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include "General/Environment.hh"
#include "MatVec/SBM_Vector.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField {


//! Class for storing a result in time and frequency domain

//! This class gets handed over a time- or frequency-domain
//! result from a multiharmonic analysis and stores the
//! SBM result vectors. Furthermore it can transform between
//! the two solution domains, via FFT.
//! The number of considered time evaluation points is given
//! by nFFT_
class MHTimeFreqResult{

public:

  //! Constructor
  MHTimeFreqResult(const UInt& N,
                   const UInt& M,
                   const Double& baseFreq,
                   const UInt& nFFT);

  //! Destructor
  ~MHTimeFreqResult();


  // ========================================================================
  //  GETTERS / SETTERS
  // ========================================================================

  void SetFrequencyResult(const SBM_Vector& freqResult);


  //! Method for evaluating a Fourier series with given Fourier-coefficients.
  //! Evaluate Fourier series at given discrete times in timeVec_
  void FourierToTime();

  //! Method for computing the Fourier series from a given time-signal.
  void TimeToFourier();



private:

  //! Called by constructor
  //! Create the time-vector, based on the base-frequency
  //! and the considered number of evaluation points in
  //! one period of the base frequency
  void CreateTimeVec();

  //! Result-vector in time-domain. Outer vector is time-step,
  //! inner one contains the spatial dof's
  StdVector< Vector<Complex> > timeResult_;

  //! Result-vector in frequency domain. Outer vector is frequency (harmonic),
  //! inner one contains the spatial dof's
  StdVector< Vector<Complex> > freqResult_;

  //! Number of time points in one period of base-frequency.
  //! The FFT of the timeResult into freqResult is based on this
  //! number of evaluation points.
  //! Also the other way round (transformation from frequency- into
  //! time-domain uses these number of steps)
  //! It is provided by the xml-parameter "numFFTPoints"
  UInt nFFT_;

  //! Number of frequencies in result
  UInt numFreqs_;


  //! Base angular frequency (1st harmonic)
  Double omega0_;

  //! Number of (spatial) dof's in computational domain
  UInt spatialSize_;

  //! Time-vector, based on the base-frequency and the considered
  //! number of evaluation points in one period of the base period
  StdVector<Double> timeVec_;

  //! Number of considered harmonics
  UInt N_;

  //! Buffer for all input data of the FFT (one time series per row)
  //! Will be overwritten with the result of the FFT.
  Matrix<Complex> data_;

};// end of class MHTimeFreqResult

}// end of namespace
#endif
