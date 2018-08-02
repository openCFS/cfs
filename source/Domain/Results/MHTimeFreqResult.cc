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
    spatialSize_ = 0;
    omega0_ = 0.0;
    timeResInit_ = false;
    nFFT_ = 0;
    N_ = 0;

    this->Init(N, M, baseFreq, nFFT);
  }

  // Default constructor
  MHTimeFreqResult::MHTimeFreqResult(){
    spatialSize_ = 0;
    omega0_ = 0.0;
    timeResInit_ = false;
    nFFT_ = 0;
    N_ = 0;
  }



  MHTimeFreqResult::~MHTimeFreqResult(){}


  void MHTimeFreqResult::Init(const UInt& N,
                              const UInt& M,
                              const Double& baseFreq,
                              const UInt& nFFT){
    N_ = N;
    timeResult_.Resize(0,0);
    freqResult_.Resize(0,0);
    nFFT_ = nFFT  ;
    omega0_ = 2 * M_PI * baseFreq;
    spatialSize_ = 0;
    // Create the time vector
    this->CreateTimeVec();
  }


  void MHTimeFreqResult::PrintTimeVector(){
    LOG_DBG3(mhtimefreqresult) << "CoefFunctionHarmBalance: Time vector: " << timeVec_.ToString();
  }


  void MHTimeFreqResult::PrintTimeResults(const UInt& index){
    std::string logString = "CoefFunctionHarmBalance: Time result for element with index ";
    logString.append( boost::lexical_cast<std::string>(index) );
    logString.append("\n \t and timesteps 0 to ");
    logString.append( boost::lexical_cast<std::string>(this->GetNumTimeSteps() - 1) );
    logString.append(" : \n");

    // construct logging output
    for(UInt i = 0; i < this->GetNumTimeSteps(); ++i){
      logString.append( boost::lexical_cast<std::string>(this->GetTimeResult(i)[index].real()) );
      logString.append(",\n");
    }
    LOG_DBG3(mhtimefreqresult) << logString <<"\n";
  }

  void MHTimeFreqResult::PrintFreqResults(const UInt& index){
    std::string logString = "CoefFunctionHarmBalance: Frequency result for element with index ";
    logString.append( boost::lexical_cast<std::string>(index) );
    logString.append("\n \t and frequency steps 0 to ");
    logString.append( boost::lexical_cast<std::string>(this->GetNumFreqSteps() - 1) );
    logString.append(" : \n");

    // construct logging output
    for(UInt i = 0; i < this->GetNumFreqSteps(); ++i){
      logString.append( boost::lexical_cast<std::string>(this->GetFreqResult(i)[index]) );
      logString.append(",\n");
    }
    LOG_DBG3(mhtimefreqresult) << logString <<"\n";
  }


  void MHTimeFreqResult::SetFrequencyResult(const SBM_Vector & freqRes){
    // If this is not the first calculation, clean stuff up
    if(spatialSize_ != 0){
      timeResult_.Resize(0,0);
      freqResult_.Resize(0,0);
    }

    // Set the number of spatial dof's as the size of the first frequency
    // result. It doesn't matter which one we consider, since they all
    // have the same size
    spatialSize_ = freqRes(0).GetSize();
    UInt numFreq = freqRes.GetSize();

    // Sanity check
    if(numFreq != N_ + 1){
      EXCEPTION("MHTimeFreqResult::SetFrequencyResult This should not happen!");
    }

    freqResult_.Resize(spatialSize_, 2 * N_ + 1);
    Complex in = 0.0;
    freqResult_.InitValue(in);
    Integer h;
    UInt ind;
    // Handle harmonics
    for(UInt i = 0; i < N_ + 1; ++i){
      // which harmonic are we considering (in the optimized version)
      h = -N_ + 2 * i;
      // where is this harmonic in the full harmonics vector (including
      // all harmonics, also the even ones)
      ind = N_ + h;
      // insert it
      for(UInt j = 0; j < spatialSize_; ++j){
        freqRes(i).GetEntry(j, freqResult_[j][ind]);
      }
    }

//    for(UInt i = 0; i < freqResult_.GetNumCols(); ++i){
//      Vector<Complex> tmp;
//      freqResult_.GetCol(tmp, i);
//      std::cout<<"Initial freqResult_["<<i<<"] = "<< freqResult_[5][i]<<std::endl;
//    }
  }

  void MHTimeFreqResult::SetTimeResult(const Vector<Double>& timeRes,
                                       const UInt& timeStep){
    // check if result was already initialized
    if(!timeResInit_){
      EXCEPTION("MHTimeFreqResult::SetTimeResult You need to initialize the time result before use via InitTimeResult(spatialSize)");
    }

    // check the spatial size
    if(timeRes.GetSize() != timeResult_.GetNumRows()){
      EXCEPTION("MHTimeFreqResult::SetTimeResult You want to set a vector with size "<< timeRes.GetSize()<<
          " into the timeResult_ matrix with only "<< timeResult_.GetNumRows()<< " entries!!");
    }

    // check the time size
    if(timeStep >= timeResult_.GetNumCols()){
      EXCEPTION("MHTimeFreqResult::SetTimeResult You want to set the time vector for step "<< timeStep<<
          " but the result allows only  "<< timeResult_.GetNumCols()<< " steps!!");
    }

    UInt i = timeStep;
    Double t;
    for(UInt j = 0; j < spatialSize_; ++j){
      timeRes.GetEntry(j, t);
      timeResult_[j][i] = (Complex)t;
    }
    //Vector<Complex>tmp;
    //timeResult_.GetCol(tmp,timeStep);
    //std::cout<<"timeResult_ row "<< timeStep<< " : "<<tmp.ToString()<<std::endl;
  }



  void MHTimeFreqResult::TimeToFourier(){
    if(timeResult_.GetNumCols() == 0) EXCEPTION("MHTimeFreqResult::TimeToFourier No time signal was set!");

    // MKL overwrites the provided time array with frequency results,
    // therefore copy the timeResult into freqResult and perform the
    // FFT on freqResult_
    freqResult_.Resize(timeResult_.GetNumRows(), 2 * N_ + 1);
    freqResult_.Init();

    deleteMeResult_.Resize(timeResult_.GetNumRows(), timeResult_.GetNumCols());
    deleteMeResult_ = timeResult_;
    //for(UInt i = 0; i < timeResult_.GetNumRows(); ++i){
    //  for(UInt j = 0; j < timeResult_.GetNumCols(); ++j ){
    //    deleteMeResult_[i][j] = timeResult_[i][j];
    //  }
    //}

    Complex tmp;

#ifdef USE_MKL
    /*========================================================
     * 1) Create Descriptor
     * ========================================================
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

    /*========================================================
     * 2) Set Value
     * ========================================================
     * Optionally adjust the descriptor configuration with a call to the
     * DftiSetValue or DftiSetValueDM function as needed. Typically, you must
     * carefully define the data storage layout for an FFT or the data distribution
     * among processes for a Cluster FFT. The configuration settings of the descriptor,
     * such as the default values, can be obtained with the DftiGetValue or DftiGetValueDM function.
     */

    status = DftiSetValue(dftHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DFTI_CHECK_STATUS(dftHandle, status);

    /* ========================================================
     * 3) Commit Descriptor
     * ========================================================
     * Commit the descriptor with a call to the DftiCommitDescriptor or
     * DftiCommitDescriptorDM function, that is, make the descriptor ready
     * for the transform computation. Once the descriptor is committed,
     * the parameters of the transform, such as the type and number of
     * transforms, strides and distances, the type and storage layout of
     * the data, and so on, are "frozen" in the descriptor.
     */

    status = DftiCommitDescriptor(dftHandle);
    DFTI_CHECK_STATUS(dftHandle, status);

    /*========================================================
     * 4) Perform FFT
     * ========================================================
     * Compute the transform with a call to the DftiComputeForward/DftiComputeBackward
     * or DftiComputeForwardDM/DftiComputeBackwardDM functions as many times as needed.
     * Because the descriptor is defined and committed separately, all that the
     * compute functions do is take the input and output data and compute the
     * transform as defined. To modify any configuration parameters for another
     * call to a compute function, use DftiSetValue followed by DftiCommitDescriptor
     * (DftiSetValueDM followed by DftiCommitDescriptorDM) or create and commit another descriptor.
     */

    // No averaging, let MKL handle the loop and
    // reorder the harmonics because the static harmonic 0
    // is located at position [0] and the output from MKL looks
    // like
    // [0   1   2   ...   N-1   N   -N    -N+1   ...   -2   -1]
    // but we need the following order
    // [-N   -N+1   ...   -1   0    1   ...    N-1    N]
    // Furthermore the size of freqResult_-columns (output from
    // the FFT) is nFFT_, therefore we have to extract the
    // static harmonic and N_ positive and N_ negative ones from
    // the matrix ... or simply set the unwanted to zero
    for(UInt i = 0; i < spatialSize_; ++i){
      //status = DftiComputeForward(dftHandle, freqResult_[i]);
      //DFTI_CHECK_STATUS(dftHandle, status);

      status = DftiComputeForward(dftHandle, deleteMeResult_[i]);
      DFTI_CHECK_STATUS(dftHandle, status);

      // Static entry for harmonic 0
      tmp = deleteMeResult_[i][0];
      // Fill the correct frequency result matrix and also scale the FFT correctly
      for(UInt k = 1; k < N_ + 1; ++k){
        freqResult_[i][N_ + k] = deleteMeResult_[i][k];
      }
      for(UInt k = 1; k < N_ + 1; ++k){
        freqResult_[i][N_ - k] = deleteMeResult_[i][nFFT_ - k];
      }
      freqResult_[i][N_] = tmp;

    }

    status = DftiFreeDescriptor(&dftHandle); // ignore errors from Free function
#endif

    freqResult_ = freqResult_ * (1.0 / ((nFFT_ )/2));

//    for(UInt i = 0; i < deleteMeResult_.GetNumCols(); ++i){
//      std::cout<<"deleteMeResult_["<<i<<"] = "<<deleteMeResult_[10][i]<<std::endl;
//    }
//    for(UInt i = 0; i < freqResult_.GetNumCols(); ++i){
//      std::cout<<"freqResult_["<<i<<"] = "<<freqResult_[10][i]<<std::endl;
//    }
 }


  void MHTimeFreqResult::FourierToTime(){
    if(timeResInit_){
      EXCEPTION("MHTimeFreqResult::FourierToTime It seems that a timeresult was already set");
    }

    // First of all, check if we have the number of spatial dof's
    if(spatialSize_ == 0) EXCEPTION("MHTimeFreqResult::FourierToTime no result was set!!")

    if( freqResult_.GetNumCols() != 2 * N_ + 1){
      EXCEPTION("MHTimeFreqResult::FourierToTime There are " << freqResult_.GetNumCols()
          << " frequencies given but \n  the number of harmonics is N_ = " << N_);
    }

    // Evaluate Fourier series at given discrete times in timeVec_

    // Since the timeVec_ size is much smaller than the number of
    // dof's in the domain, our outer loop will be the timeVec_ entries
    // and the inner loop will be the parallel loop over dof's
    timeResult_.Resize(spatialSize_ ,nFFT_ );
    timeResult_.InitValue(0.0);

    // Outer loop over discrete times
    Double t = 0.0;
    int h = 0;
    for(UInt i = 0; i < timeVec_.GetSize(); ++i ){
      t = timeVec_[i];
      for(UInt k = 0; k < 2 * N_ + 1; ++k){
        // harmonic number
        h = k - N_;
        for(UInt j = 0; j < spatialSize_; ++j){
          // multiplication with 0.5 due to double sided spectrum
          timeResult_[j][i] += freqResult_[j][k] * (cos(h * omega0_ * t) + Complex(0.0,1.0)*sin(h * omega0_ *t)) * 0.5;
        }
      }
    }// loop over time array

//    for(UInt i = 0; i < timeResult_.GetNumCols(); ++i){
//      std::cout<<"timeResult_["<<i<<"] = "<<timeResult_[10][i]<<std::endl;
//    }

  }


  void MHTimeFreqResult::CreateTimeVec(){
    // -1 because the last "signal point" in the time signal
    // must not be equal to the first one
    timeVec_.Resize(nFFT_ , 0.0);

    // Period length of base frequency
    Double baseFreq = omega0_ / (2 * M_PI);
    Double tBase = 1.0 / baseFreq;
    Double dT = tBase / (nFFT_);
    Double tAc = 0.0;
    // Fill vector with equally spaced timesteps over [0, tBase]
    for(auto& t : timeVec_){
      t = tAc;
      tAc += dT;
    }

    LOG_DBG3(mhtimefreqresult) << "MHTimeFreqResult: Constructed timeline = [" << timeVec_.ToString()<<" ] \n" ;
  }


}
