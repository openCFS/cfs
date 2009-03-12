#include <cstdlib>

#include "ResultCache.hh"
#include "resultHandler.hh"
#include "simInput.hh"
#include "Domain/domain.hh"

namespace CoupledField
{

bool ResultCache::resultValid_ = false;

bool ResultCache::stepsValid_ = false;

Double ResultCache::curStepVal_ = -1.0;

ResultCache::OutputType ResultCache::outType_ = OUT_REAL;

BaseMatrix::EntryType ResultCache::entryType_ = BaseMatrix::DOUBLE;

UInt ResultCache::dof_ = 0;

UInt ResultCache::index_ = 0;

UInt ResultCache::lastStepNum_ = 0;

UInt ResultCache::numDofs_ = 1;

UInt ResultCache::sequenceStep_ = 1;

std::string ResultCache::entityName_;

std::string ResultCache::readerId_ = "default";

std::map<UInt, Double> ResultCache::stepValues_;

SolutionType ResultCache::solType_ = NO_SOLUTION_TYPE;

Vector<Double> ResultCache::resCacheDouble_;

Vector<Complex> ResultCache::resCacheComplex_;

void ResultCache::SetOutputType(OutputType outType) {
  outType_ = outType;
}

void ResultCache::SetDof(UInt dof) {
  dof_ = dof;
}

void ResultCache::SetEntity(const std::string& entityName) {
  if (entityName != entityName_) {
    resultValid_ = false;
    stepsValid_ = false;
  }
  entityName_ = entityName;
}

void ResultCache::SetInfo(OutputType outType,
                          UInt dof,
                          const std::string& entityName,
                          SolutionType solType,
                          Double stepValue) {
  outType_ = outType;
  dof_ = dof;
  if (entityName != entityName_) {
    resultValid_ = false;
    stepsValid_ = false;
  }
  entityName_ = entityName;
  if (solType != solType_) {
    resultValid_ = false;
    stepsValid_ = false;
  }
  solType_ = solType;
  if (stepValue != curStepVal_)
    resultValid_ = false;
  curStepVal_ = stepValue;
}

void ResultCache::SetSolution(SolutionType solType) {
  if (solType != solType_) {
    resultValid_ = false;
    stepsValid_ = false;
  }
  solType_ = solType;
}

void ResultCache::SetStepValue(Double stepValue) {
  if (stepValue != curStepVal_)
    resultValid_ = false;
  curStepVal_ = stepValue;
}

Double ResultCache::GetResult(const char* inputId, Double sequenceStep) {
  ResultHandler* resHandler = domain->GetResultHandler();
  std::string readerId(inputId);

  // do we need to load a new result?
  if (!resultValid_ || (readerId != readerId_)
      || (sequenceStep != sequenceStep_)) {
    if (!stepsValid_) {
      LoadStepValues(readerId, (UInt) sequenceStep);
      stepsValid_ = true;
      lastStepNum_ = 0;
    }

    std::map<UInt, Double>::const_iterator it
      = stepValues_.lower_bound(lastStepNum_);
    std::map<UInt, Double>::const_iterator itEnd = stepValues_.end();
    UInt curStepNum = it->first;

    while (it != itEnd) {
      if (it->second >= curStepVal_)
        break;
      curStepNum = (++it)->first;
    }
    if (it != itEnd) {
      Double dt = it->second - curStepVal_;
      if (curStepVal_ - (--it)->second < dt)
        curStepNum = it->first;
    }
    else {
        curStepNum = (--it)->first;
    }

    shared_ptr<BaseResult> tmpRes = resHandler->GetResult(readerId,
                                                          (UInt) sequenceStep,
                                                          curStepNum,
                                                          solType_,
                                                          entityName_);
    entryType_ = tmpRes->GetEntryType();
    if (entryType_ == BaseMatrix::COMPLEX) {
      resCacheComplex_ = dynamic_cast<Result<Complex>&>(*tmpRes).GetVector();
    }
    else {
      resCacheDouble_ = dynamic_cast<Result<Double>&>(*tmpRes).GetVector();
    }
    readerId_ = readerId;
    sequenceStep_ = (UInt) sequenceStep;
    resultValid_ = true;
  }


  if (entryType_ == BaseMatrix::COMPLEX) {
#ifndef NDEBUG
    if (index_*numDofs_ + dof_ - 1 >= resCacheComplex_.GetSize()) {
      std::string resultName;
      Enum2String(solType_, resultName);
      EXCEPTION("Index of '" << resultName << "' from input file (id="
          << readerId_ << ") out of bounds.");
    }
#endif

    switch (outType_) {
    case OUT_REAL:
      return resCacheComplex_[index_*numDofs_ + dof_ -1].real();
    case OUT_IMAG:
      return resCacheComplex_[index_*numDofs_ + dof_ -1].imag();
    case OUT_AMPL:
      return std::abs(resCacheComplex_[index_*numDofs_ + dof_ -1]);
    case OUT_PHASE:
      // we need to convert radians to degrees, because this goes through
      // the math parser!!!
      return std::arg(resCacheComplex_[index_*numDofs_ + dof_ -1])*180.0/PI;
    }
  }
  else {
#ifndef NDEBUG
     if (index_*numDofs_ + dof_ - 1 >= resCacheDouble_.GetSize()) {
      std::string resultName;
      Enum2String(solType_, resultName);
      EXCEPTION("Index of '" << resultName << "' from input file (id="
          << readerId_ << ") out of bounds.");
    }
#endif

    return resCacheDouble_[index_*numDofs_ + dof_ - 1];
  }

  return 0.0;
}

void ResultCache::LoadStepValues(const std::string& readerId,
                                 UInt sequenceStep) {
  shared_ptr<SimInput> inputReader = domain->GetResultHandler()
      ->GetInputReader(readerId);
  StdVector< shared_ptr<ResultInfo> > resInfos;

  try {
    inputReader->GetResultTypes(sequenceStep, resInfos);
  }
  catch (Exception& ex) {
    EXCEPTION("Could not get data from input (id='" << readerId << "')");
  }
  UInt numResults = resInfos.GetSize();
  UInt iRes;
  for (iRes = 0; iRes < numResults; ++iRes) {
    if (resInfos[iRes]->resultType == solType_) {
      break;
    }
  }
  if (iRes >= resInfos.GetSize()) {
    std::string resultName;
    Enum2String(solType_, resultName);
    EXCEPTION("Result '" << resultName << "' not found in input file (id="
        << readerId << ")");
  }

  numDofs_ = resInfos[iRes]->dofNames.GetSize();
  inputReader->GetStepValues(sequenceStep, resInfos[iRes], stepValues_);

}

} // namespace CoupledField
