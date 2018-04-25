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
#include "MatVec/Matrix.hh"
#include "General/Environment.hh"
#include "MatVec/SBM_Vector.hh"


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

  //! Default Constructor
  MHTimeFreqResult();

  //! Constructor
  MHTimeFreqResult(const UInt& N,
                   const UInt& M,
                   const Double& baseFreq,
                   const UInt& nFFT);

  //! Destructor
  ~MHTimeFreqResult();


  //! Initialization, does the same as the constructor
  void Init(const UInt& N,
       const UInt& M,
       const Double& baseFreq,
       const UInt& nFFT);

  // ========================================================================
  //  GETTERS / SETTERS
  // ========================================================================

  void SetFrequencyResult(const SBM_Vector& freqResult);

  Matrix<Complex> GetFullFrequencyResult(){ return freqResult_;}

  Matrix<Complex> GetFullTimeResult(){ return timeResult_;}

  Vector<Complex> GetTimeResult(int tStep){
    Vector<Complex> ret;
    timeResult_.GetCol(ret, tStep);
    return ret;
  }

  Vector<Complex> GetFreqResult(int tStep){
    Vector<Complex> ret;
    freqResult_.GetCol(ret, tStep);
    return ret;
  }

  UInt GetNumTimeSteps(){ return timeResult_.GetNumRows(); }

  UInt GetNumFreqSteps(){ return freqResult_.GetNumRows(); }


  // ========================================================================
  //  SHIFT BETWEEN TIME AND FREQUENCY DOMAIN
  // ========================================================================

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

  //! Result-vector in time-domain.
  //! One time series per row
  Matrix<Complex> timeResult_;

  //! Result-vector in frequency domain.
  //! One time series per row
  Matrix<Complex> freqResult_;

  //! We need this result to fill the correct freqResult_ matrix
  //! after the MKL FFT computation
  Matrix<Complex> deleteMeResult_;



  //! Number of time points in one period of base-frequency.
  //! The FFT of the timeResult into freqResult is based on this
  //! number of evaluation points.
  //! Also the other way round (transformation from frequency- into
  //! time-domain uses these number of steps)
  //! It is provided by the xml-parameter "numFFTPoints"
  UInt nFFT_;

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
