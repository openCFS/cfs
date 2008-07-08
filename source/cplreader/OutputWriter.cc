#include "OutputWriter.hh"

namespace CoupledField
{
  OutputWriter::OutputWriter() :
    dim_(3) {
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
