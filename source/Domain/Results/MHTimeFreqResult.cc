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

#include <def_use_blas.hh>
#ifdef USE_MKL
#include "mkl_dfti.h"
#endif

// Copy from cfsdat FFtFilter.cc
#define DFTI_CHECK_STATUS(handle, status) \
  if (status != 0) { \
    if (handle) DftiFreeDescriptor(&handle); \
    char *msg = DftiErrorMessage(status); \
    EXCEPTION("FFT failed: " << msg); \
  }


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
    nFFT_ = nFFT + 1 ;
    omega0_ = 2 * M_PI * baseFreq;
    spatialSize_ = 0;
    //numFreqs_ = nFFT/2 - 1;
    numFreqs_ = N_;

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


  void MHTimeFreqResult::TimeToFourier(){
    if(timeResult_.GetSize() == 0) EXCEPTION("MHTimeFreqResult::TimeToFourier No time signal was set!");

    // Fill the data matrix array with time values
    // Row i contains the time-results for point i
    data_.Resize(spatialSize_, nFFT_);
    data_.Init();

    //TODO Window function could be applied here
    for(UInt tStep = 0; tStep < timeVec_.GetSize(); ++tStep){
      for (UInt i = 0; i < spatialSize_; ++i) {
        data_[i][tStep] = timeResult_[tStep][i];
        std::cout<<"data_["<<i<<"]["<<tStep<<"]="<<data_[i][tStep]<<std::endl;
      }
    }


#ifdef USE_MKL
    /* ========================================================
     * 1) Create Descriptor
     * ======================================================== */

    /*
     * Allocate a fresh descriptor for the problem with a call to the
     * DftiCreateDescriptor or DftiCreateDescriptorDM function.
     * The descriptor captures the configuration of the transform,
     * such as the dimensionality (or rank), sizes, number of transforms,
     * memory layout of the input/output data (defined by strides)
     * and scaling factors. Many of the configuration settings are assigned
     * default values in this call which you might need to modify in your application.
     */

    // create MKL DFT descriptor
    DFTI_DESCRIPTOR_HANDLE dftHandle = 0;
    // Dimensionality of the time signal
    UInt dim = 1;
    // Double or single precision of the fft
    DFTI_CONFIG_VALUE precision = DFTI_DOUBLE;
    // We have a complex valued time signal but only the real value should be present
    DFTI_CONFIG_VALUE domain = DFTI_COMPLEX;

    MKL_LONG status = DftiCreateDescriptor(&dftHandle, precision, domain, dim, nFFT_);
    DFTI_CHECK_STATUS(dftHandle, status);

    /* ========================================================
     * 2) Set Value
     * ======================================================== */

    /*
     * Optionally adjust the descriptor configuration with a call to the
     * DftiSetValue or DftiSetValueDM function as needed. Typically, you must
     * carefully define the data storage layout for an FFT or the data distribution
     * among processes for a Cluster FFT. The configuration settings of the descriptor,
     * such as the default values, can be obtained with the DftiGetValue or DftiGetValueDM function.
     */

    status = DftiSetValue(dftHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DFTI_CHECK_STATUS(dftHandle, status);

    // Using a packed storage format
    //status = DftiSetValue(dftHandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
    //DFTI_CHECK_STATUS(dftHandle, status);


    // Set how many spatial dof's we consider (number of Fourier transforms)
    //status = DftiSetValue(dftHandle, DFTI_NUMBER_OF_TRANSFORMS, spatialSize_);
    //DFTI_CHECK_STATUS(dftHandle, status);

    // Set number of time steps
    //status = DftiSetValue(dftHandle, DFTI_INPUT_DISTANCE, nFFT_);
    //DFTI_CHECK_STATUS(dftHandle, status);

    // Set how many harmonics we want to consider
    //status = DftiSetValue(dftHandle, DFTI_OUTPUT_DISTANCE, numFreqs_);
    //DFTI_CHECK_STATUS(dftHandle, status);

    /* ========================================================
     * 3) Commit Descriptor
     * ======================================================== */

    /*Commit the descriptor with a call to the DftiCommitDescriptor or
     * DftiCommitDescriptorDM function, that is, make the descriptor ready
     * for the transform computation. Once the descriptor is committed,
     * the parameters of the transform, such as the type and number of
     * transforms, strides and distances, the type and storage layout of
     * the data, and so on, are "frozen" in the descriptor.
     */
    status = DftiCommitDescriptor(dftHandle);
    DFTI_CHECK_STATUS(dftHandle, status);



    /* ========================================================
     * 4) Perform FFT
     * ======================================================== */

    /*
     * Compute the transform with a call to the DftiComputeForward/DftiComputeBackward
     * or DftiComputeForwardDM/DftiComputeBackwardDM functions as many times as needed.
     * Because the descriptor is defined and committed separately, all that the
     * compute functions do is take the input and output data and compute the
     * transform as defined. To modify any configuration parameters for another
     * call to a compute function, use DftiSetValue followed by DftiCommitDescriptor
     * (DftiSetValueDM followed by DftiCommitDescriptorDM) or create and commit another descriptor.
     */

      // No averaging, let MKL handle the loop


    for(UInt i = 0; i < spatialSize_; ++i){
      status = DftiComputeForward(dftHandle, data_[i]);
      DFTI_CHECK_STATUS(dftHandle, status);
    }

    status = DftiFreeDescriptor(&dftHandle); // ignore errors from Free function



#endif






    // fill the freqResult
    freqResult_.Resize(2 * N_ + 1);
    for(UInt i = 0; i < freqResult_.GetSize(); ++i){
      Vector<Complex>& fRes = freqResult_[i];
      fRes.Resize(spatialSize_);
      if(i != N_){
        for(UInt t = 0; t < spatialSize_; ++t){
          fRes[i] = data_[t][i+1];
        }
      }else{
        // The static component is at position 0
        for(UInt t = 0; t < spatialSize_; ++t){
          fRes[i] = data_[t][0];
        }
      }
    }



    std::cout<<"freqResult_[";
    for(UInt i = 0; i < freqResult_.GetSize(); ++i){
      std::cout<<freqResult_[i][0];
    }
    std::cout<<"]"<<std::endl;

    std::cout<<"timeResult_[";
    for(UInt i = 0; i < timeResult_.GetSize(); ++i){
      std::cout<<timeResult_[i][0];
    }
    std::cout<<"]"<<std::endl;




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
    timeVec_.Resize(nFFT_, 0.0);

    // Period length of base frequency
    Double baseFreq = omega0_ / (2 * M_PI);
    Double tBase = 1.0 / baseFreq;
    Double dT = tBase / (nFFT_ - 1);
    Double tAc = 0.0;
    // Fill vector with equally spaced timesteps over [0, tBase]
    for(auto& t : timeVec_){
      t = tAc;
      tAc += dT;
    }

    LOG_DBG3(mhtimefreqresult) << "MHTimeFreqResult: Constructed timeline = [" << timeVec_.ToString()<<" ] \n" ;
  }

}
