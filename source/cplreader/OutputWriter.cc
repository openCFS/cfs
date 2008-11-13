// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "OutputWriter.hh"

namespace CoupledField
{
  OutputWriter::OutputWriter() :
    dim_(3),
    doneWriting_(false)
  {
  }

  OutputWriter::~OutputWriter() {
  }

  void OutputWriter::BeginStep(UInt stepNum, Double stepValue) {
    actStepNum_ = stepNum;
    actStepValue_ = stepValue;
  }

  void OutputWriter::EndStep() {
    std::vector<std::string>::iterator it, end;

    it = actStepResults_.begin();
    end = actStepResults_.end();

    for( ; it != end; it++ ) {
      stepNums_[*it].push_back(actStepNum_);
      stepValues_[*it].push_back(actStepValue_);
    }
  }
}
