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
    timeResult_.ResetEntryType(BaseMatrix::DOUBLE);

    // frequency results can only be complex
    timeResult_.ResetEntryType(BaseMatrix::COMPLEX);
    nFFT_ = nFFT;
    baseFreq_ = baseFreq;

    // Create the time vector
    CreateTimeVec();
  }

  MHTimeFreqResult::~MHTimeFreqResult(){}



  void MHTimeFreqResult::CreateTimeVec(){
    timeVec_.Resize(nFFT_ * N_ + 1, 0.0);

    // Period length of base frequency
    Double tBase = 1.0 / baseFreq_;
    Double dT = tBase / nFFT_;
    Double tAc = 0.0;
    // Fill vector with equally spaced timesteps over [0, tBase]
    for(auto& t : timeVec_){
      t = tAc;
      tAc += dT;
    }

  }

}
