// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MHTimeFreqResult.cc
 *       \brief    Storage container for multiharmonic results
 *
 *       \date     Mar 9, 2018
 *       \author   kroppert
 */
//================================================================================================

#include "MHTimeFreqResult.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include <complex>
#include <math.h>

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif


namespace CoupledField {
  // declare logging stream
  DECLARE_LOG(mhtimefreqresult)
  DEFINE_LOG(mhtimefreqresult, "mhtimefreqresult")

  MHTimeFreqResult::MHTimeFreqResult(const UInt& N,
                                     const UInt& M,
                                     const Double& baseFreq,
                                     const UInt& nFFT){
    N_ = N;
    // time results can only be real valued
    timeResult_.Resize(0);
    // frequency results can only be complex
    freqResult_.Resize(0);
    //freqResult_.ResetEntryType(BaseMatrix::COMPLEX);
    nFFT_ = nFFT;
    omega0_ = 2 * M_PI * baseFreq;
    spatialSize_ = 0;

    // Create the time vector
    CreateTimeVec();
  }

  MHTimeFreqResult::~MHTimeFreqResult(){}



  void MHTimeFreqResult::SetFrequencyResult(const SBM_Vector & freqRes){
    UInt s = freqRes.GetSize();
    freqResult_.Resize(s);
    for(UInt i = 0; i < s; ++i){
      freqResult_[i] = freqRes(i); // SBM is one-based
    }


    // Set the number of spatial dof's as the size of the first frequency
    // result. It doesn't matter which one we consider, since they all
    // have the same size
    spatialSize_ = freqResult_[0].GetSize();
  }



  void MHTimeFreqResult::FourierToTime(){
    // First of all, check if we have the number of spatial dof's
    if(spatialSize_ == 0) EXCEPTION("MHTimeFreqResult::FourierToTime no result was set!!")

    // Evaluate Fourier series at given discrete times in timeVec_

    // Since the timeVec_ size is much smaller than the number of
    // dof's in the domain, our outer loop will be the timeVec_ entries
    // and the inner loop will be the parallel loop over dof's
    timeResult_.Resize(timeVec_.GetSize());

    // Outer loop over discrete times
    for(UInt i = 0; i < timeVec_.GetSize(); ++i ){
      Double t = timeVec_[i];
      Vector<Complex>& A = timeResult_[i];
      A.Resize( spatialSize_ );

      //TODO: Performance of the openMP and singlethread version was not measured yet
#ifdef USE_OPENMP
      StdVector< Complex* > local(CFS_NUM_THREADS);
#pragma omp parallel num_threads(CFS_NUM_THREADS)
      {
        UInt np = omp_get_num_threads();
        std::vector< Complex > localA(spatialSize_);
        local[omp_get_thread_num()] = localA.data();

        // add values to local array
        for(UInt k = 0; k < 2* N_ + 1; ++k){
          // SBMVector index is 1-based!!
          Vector<Complex> & fR = freqResult_[k];
          // harmonic number
          int h = k - N_;
          Complex c = (cos(h * omega0_ * t) + Complex(0.0,1.0)*sin(h * omega0_ *t));
#pragma omp for
          for(UInt j = 0; j < spatialSize_; ++j){
            localA[j] += fR[j] * c;
          }
        }
        // implicit barrier ensures all local copies are ready for aggregation

        // aggregate local copies into global array
#pragma omp for
        for(UInt k = 0; k < spatialSize_; ++k){
          for(UInt p = 0; p < np; ++p){
            A[k] += local[p][k];
          }
        }
        // implicit barrier ensures no local copy is deleted before aggregation is done

      }// end of omp parallel
#else
      for(UInt k = 0; k < 2 * N_ + 1; ++k){
        Vector<Complex> & fR = freqResult_[k];
        // harmonic number
        int h = k - N_;
        for(UInt j = 0; j < spatialSize_; ++j){
          A[j] += fR[j] * (cos(h * omega0_ * t) + Complex(0.0,1.0)*sin(h * omega0_ *t));
        }

      }
#endif

    }// loop over time array


    if(IS_LOG_ENABLED(mhtimefreqresult, dbg)){
      UInt i = 0;
      for(auto& t : timeResult_){
        LOG_DBG3(mhtimefreqresult) << "timeResult[" << i << "]:REAL: " << t.GetPart(Global::REAL).ToString();
        LOG_DBG3(mhtimefreqresult) << "timeResult[" << i << "]:IMAG: " << t.GetPart(Global::IMAG).ToString();
        ++i;
      }
    }

  }




  void MHTimeFreqResult::CreateTimeVec(){
    timeVec_.Resize(nFFT_ * N_ + 1, 0.0);

    // Period length of base frequency
    Double baseFreq = omega0_ / (2 * M_PI);
    Double tBase = 1.0 / baseFreq;
    Double dT = tBase / nFFT_;
    Double tAc = 0.0;
    // Fill vector with equally spaced timesteps over [0, tBase]
    for(auto& t : timeVec_){
      t = tAc;
      tAc += dT;
    }

    LOG_DBG3(mhtimefreqresult) << "MHTimeFreqResult: Constructed timeline = [" << timeVec_.ToString()<<" ] \n" ;
  }

}
