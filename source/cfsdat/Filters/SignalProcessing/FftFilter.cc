/*
 * FftFilter.cc
 *
 *  Created on: Dec 5, 2016
 *      Author: jens
 */

#include <cfsdat/Filters/SignalProcessing/FftFilter.hh>
#define _USE_MATH_DEFINES
#include <cmath>

// check if Intel MKL is available
#include <def_use_blas.hh>
#ifdef USE_MKL
# include <mkl.h>
#endif

#define DFTI_CHECK_STATUS(handle, status) \
  if (status != 0) { \
    if (handle) DftiFreeDescriptor(&handle); \
    char *msg = DftiErrorMessage(status); \
    EXCEPTION("FFT failed: " << msg); \
  }

namespace CFSDat {

CF::Enum<FftFilter::WindowType> FftFilter::WindowTypeEnum = FftFilter::InitWindowTypeEnum();

CF::Enum<FftFilter::WindowType> FftFilter::InitWindowTypeEnum() {
  CF::Enum<FftFilter::WindowType> ret;
  ret.SetName("WindowType");
  ret.Add(WIN_RECT, "rect");
  ret.Add(WIN_BLACKMAN, "blackman");
  ret.Add(WIN_BLACKMAN_HARRIS, "blackmanHarris");
  ret.Add(WIN_HAMMING, "hamming");
  ret.Add(WIN_HANN, "hann");
  ret.Add(WIN_NUTTAL, "nuttal");
  return ret;
}

void FftFilter::GetWindow(Vector<Double> & winVec, WindowType winType, UInt N,
                          Double param)
{
  Double p;
  winVec.Resize(N);
  
  switch (winType) {
    case WIN_RECT:
      for (UInt i = 0; i < N; ++i) {
        winVec[i] = 1.0;
      }
      break;
      
    case WIN_BLACKMAN:
      for (UInt i = 0; i < N; ++i) {
        p = 2.0 * M_PI * (Double) i / (Double) N;
        winVec[i] = 0.42 - 0.5 * std::cos(p) + 0.08 * std::cos(2.0*p);
      }
      break;
      
    case WIN_BLACKMAN_HARRIS:
      for (UInt i = 0; i < N; ++i) {
        p = 2.0 * M_PI * (Double) i / (Double) N;
        winVec[i] = 0.35875 - 0.48829 * std::cos(p)
                            + 0.14128 * std::cos(2.0*p)
                            - 0.01168 * std::cos(3.0*p);
      }
      break;
      
    case WIN_HAMMING:
      for (UInt i = 0; i < N; ++i) {
        p = 2.0 * M_PI * (Double) i / (Double) N;
        winVec[i] = 0.54 - 0.46 * std::cos(p);
      }
      break;
      
    case WIN_HANN:
      for (UInt i = 0; i < N; ++i) {
        p = 2.0 * M_PI * (Double) i / (Double) N;
        winVec[i] = 0.5 - 0.5 * std::cos(p);
      }
      break;
      
    case WIN_NUTTAL:
      for (UInt i = 0; i < N; ++i) {
        p = 2.0 * M_PI * (Double) i / (Double) N;
        winVec[i] = 0.355768 - 0.487396 * std::cos(p)
                             + 0.144232 * std::cos(2.0*p)
                             - 0.012604 * std::cos(3.0*p);
      }
      break;
  }
}

FftFilter::FftFilter(UInt numWorkers, PtrParamNode config, PtrResultManager resMan)
: BaseFilter(numWorkers, config, resMan)
{
#ifndef USE_MKL
  EXCEPTION("Compile with Intel MKL for FFT support");
#endif
  filtStreamType_ = ACCUMULATIVE_FILTER;
  
  // Parse config
  nFFT_ = params_->Get("numSamples")->As<UInt>();
  windowType_ = WindowTypeEnum.Parse(params_->Get("window")->Get("type")->As<std::string>());
  windowParam_ = params_->Get("window")->Get("param")->As<Double>();
  overlap_ = params_->Get("overlap")->As<Double>();
  numFreqs_ = nFFT_/2 + 1;
  
  numSamples_ = 0;
  windowLength_ = 0;
  deltaT_ = 0.0;
  deltaF_ = 0.0;
  
  // parse input and output quantities
  PtrParamNode singleResult = params_->Get("singleResult");
  std::string inRes = singleResult->Get("inputQuantity")->Get("resultName")->As<std::string>();
  std::string outRes = singleResult->Get("outputQuantity")->Get("resultName")->As<std::string>();
  upResNames.insert(inRes);
  filtResNames.insert(outRes);
}

ResultIdList FftFilter::SetUpstreamResults() {
  return SetDefaultUpstreamResults();
}

void FftFilter::AdaptFilterResults() {
  uuids::uuid inId = *upResIds.Begin();
  uuids::uuid outId = *filterResIds.Begin();
  
  // copy result info from input and switch to complex 
  resultManager_->CopyResultData(inId, outId);
  resultManager_->SetDType(outId, ExtendedResultInfo::COMPLEX);
  resultManager_->SetValid(outId);

  // get input timeline and store it
  StdVector<Double> timeline = resultManager_->GetTimeLine(inId);
  numSamples_ = timeline.GetSize();
  if (numSamples_ < 2) {
    EXCEPTION("FFT filter needs more than one time step to operate on");
  }
  // re-initialize step-to-value map from input
  globalStepValueMap_.clear();
  for (UInt i = 0; i < numSamples_; ++i) {
    globalStepValueMap_[i+1] = timeline[i];
  }
  myStepIter_ = globalStepValueMap_.begin();
  
  // assume equidistant time steps
  deltaT_ = timeline[1] - timeline[0];
  
  // compute frequencies from timeline of input
  deltaF_ = 1.0 / (deltaT_ * nFFT_);
  StdVector<Double> freqline(numFreqs_);
  for (UInt i = 0; i < numFreqs_; ++i) {
    freqline[i] = (Double) i * deltaF_;
  }
  resultManager_->SetTimeLine(outId, freqline);
  
  // initialize window function vector
  windowLength_ = std::min(numSamples_, nFFT_);
}

bool FftFilter::Run() {
  
  // If Run() is called for the first time, we have to grind through the whole
  // timeline of the upstream pipeline to obtain the time signals.
  if (data_.GetNumRows() == 0 || data_.GetNumCols() == 0) {
    
    // compute window function
    Vector<Double> winVec;
    GetWindow(winVec, windowType_, windowLength_, windowParam_);
    Double winFactor = 1.0 / winVec.Sum();
        
    // determine if we have to do averaging
    UInt numBlocks = 1;
    Integer blockOffset = 0;
    if (windowLength_  < numSamples_) {
      blockOffset = (Integer) std::floor((1.0-overlap_)*(Double)windowLength_);
      if (blockOffset <= 0) {
        EXCEPTION("Window overlap " << overlap_ << " is too large");
      }
      numBlocks = ((UInt) std::floor((Double)(numSamples_-windowLength_) /
                                     (Double)blockOffset)) + 1;
    }
    
    // get equation numbers of input
    uuids::uuid inId = *upResIds.Begin();
    StdVector<UInt> eqnVec;
    Vector<Double> &resVec = resultManager_->GetResultVector<Double>(inId, eqnVec);
    UInt numPoints = eqnVec.GetSize() == 0 ? resVec.GetSize() : eqnVec.GetSize();
    
    // initialize data matrix
    // 2*numFreqs_ == numSamples_+1 if numSamples is an even number
    data_.Resize(numPoints, 2*numFreqs_);
    data_.Init();
    
    // loop over all time steps of input and obtain data
    //UInt numSources = sources_.GetSize();
    std::map<UInt, Double>::iterator timeIt = globalStepValueMap_.begin(),
                                     timeEnd = globalStepValueMap_.end();
    for ( ; timeIt != timeEnd; ++timeIt) {
      /**
      resultManager_->SetStepValue(inId, timeIt->second);
      resultManager_->ActivateResult(inId);
      for (UInt i = 0; i < numSources; ++i) {
        sources_[i]->Run();
      }
      resVec = resultManager_->GetResultVector<Double>(inId, eqnVec);
      **/
      resVec = GetUpstreamResultVector<Double>(inId, timeIt->second);
      assert(resVec.GetSize() == numPoints);
      
      UInt step = timeIt->first - 1;
      if (numBlocks == 1) {
        // if there is just one block of data, apply the window function now
        #pragma omp parallel for num_threads(CFS_NUM_THREADS)
        for (UInt i = 0; i < numPoints; ++i) {
          data_[i][step] = resVec[i] * winVec[step];
        }
      }
      else {
        #pragma omp parallel for num_threads(CFS_NUM_THREADS)
        for (UInt i = 0; i < numPoints; ++i) {
          data_[i][step] = resVec[i];
        }
      }
      resultManager_->DeactivateResult(inId);
    }
    
#ifdef USE_MKL
    // create MKL DFT descriptor
    DFTI_DESCRIPTOR_HANDLE dftHandle = 0;
    MKL_LONG status = DftiCreateDescriptor(&dftHandle, DFTI_DOUBLE, DFTI_REAL, 1, nFFT_);
    DFTI_CHECK_STATUS(dftHandle, status);
    status = DftiSetValue(dftHandle, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
    DFTI_CHECK_STATUS(dftHandle, status);
    status = DftiSetValue(dftHandle, DFTI_PACKED_FORMAT, DFTI_CCE_FORMAT);
    DFTI_CHECK_STATUS(dftHandle, status);
    if (numBlocks == 1) {
      // If there is just one block we don't need to do averaging.
      // Let the MKL handle the loop over all points.
      status = DftiSetValue(dftHandle, DFTI_NUMBER_OF_TRANSFORMS, numPoints);
      DFTI_CHECK_STATUS(dftHandle, status);
      status = DftiSetValue(dftHandle, DFTI_INPUT_DISTANCE, data_.GetNumCols());
      DFTI_CHECK_STATUS(dftHandle, status);
      status = DftiSetValue(dftHandle, DFTI_OUTPUT_DISTANCE, numFreqs_);
      DFTI_CHECK_STATUS(dftHandle, status);
    }
    status = DftiCommitDescriptor(dftHandle);
    DFTI_CHECK_STATUS(dftHandle, status);
    
    if (numBlocks == 1) {
      // No averaging, let MKL handle the loop
      status = DftiComputeForward(dftHandle, data_[0]);
      DFTI_CHECK_STATUS(dftHandle, status);
    }
    else {
      // Allocate input and output buffers needed for averaging
      Vector<Double> inVec(2*numFreqs_), outVec(numFreqs_);
      Double avgFactor = 1.0 / (Double) numBlocks;
      
      for (UInt point = 0; point < numPoints; ++point) {
        inVec.Init();
        outVec.Init();
        
        for (UInt block = 0; block < numBlocks; ++block) {
          // Apply window function to current block
          UInt offset = block * blockOffset;
          for (UInt i = 0; i < windowLength_; ++i) {
            inVec[i] = data_[point][i+offset] * winVec[i];
          }
          
          // Compute FFT
          status = DftiComputeForward(dftHandle, inVec.GetPointer());
          DFTI_CHECK_STATUS(dftHandle, status);
          
          // Add result of block to overall result
          for (UInt i = 0; i < numFreqs_; ++i) {
            outVec[i] += std::sqrt(inVec[2*i]*inVec[2*i] + inVec[2*i+1]*inVec[2*i+1]);
          }
        }
        
        // store accumulative result back in to data_ matrix
        for (UInt i = 0; i < numFreqs_; ++i) {
          data_[point][2*i  ] = outVec[i] * avgFactor; // store magnitude in real part
          data_[point][2*i+1] = 0.0; // imaginary part is always zero for averaged spectrum
        }
      }
    }
    
    status = DftiFreeDescriptor(&dftHandle); // ignore errors from Free function
#endif
    
    // scale result of FFT to compensate for effects of windowing and one-sided spectrum
    const Double winFactor2 = 2.0 * winFactor;
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt point = 0; point < numPoints; ++point) {
      data_[point][0] *= winFactor;
      for (UInt freq = 1; freq < 2*numFreqs_; ++freq) {
        data_[point][freq] *= winFactor2;
      }
    }
  }
  
  // for now there can be only one result
  uuids::uuid outRes = *filterResIds.Begin();
  StdVector<UInt> eqns;
  Vector<Complex> &resVec = resultManager_->GetResultVector<Complex>(outRes, eqns);
  UInt numPoints =  data_.GetNumRows();

  // copy results
  UInt stepIdx = myStepIter_->first - 1;
  for (UInt point = 0; point < numPoints; ++point) {
    resVec[point] = Complex(data_[point][2*stepIdx], data_[point][2*stepIdx+1]);
  }
  ++myStepIter_;
  
  return true;
}

} /* namespace CAATool */
