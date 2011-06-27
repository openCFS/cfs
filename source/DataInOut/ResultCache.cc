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

std::map<UInt, UInt> ResultCache::nodeNum2Index_;


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
    entityName_ = entityName;
  }
}

void ResultCache::SetInfo(OutputType outType,
                          UInt dof,
                          const std::string& entityName,
                          SolutionType solType) {
  outType_ = outType;
  dof_ = dof;
  SetEntity(entityName);
  SetSolution(solType);
}

void ResultCache::SetSolution(SolutionType solType) {
  if (solType != solType_) {
    resultValid_ = false;
    stepsValid_ = false;
    solType_ = solType;
  }
}

void ResultCache::SetStepValue(Double stepValue) {
  if (stepValue != curStepVal_) {
    resultValid_ = false;
    curStepVal_ = stepValue;
  }
}

Double ResultCache::GetResult( const char* inputId, Double sequenceStep,
                               Double dofUser ) {
  UInt dof = 0, locIndex;
  ResultHandler* resHandler = domain->GetResultHandler();
  std::string readerId(inputId);

  // do we need to load a new result?
  if (!resultValid_ || (readerId != readerId_)
      || (sequenceStep != sequenceStep_))
  {
    if (!stepsValid_)
    {
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

  // Determine dof index
  if ( dofUser <= 0.0 ) {
    dof = dof_;
  } else {
    dof = (UInt) ceil(dofUser);
  }
  // check that dof makes sense
  if ( dof <= 0 || dof > numDofs_ )
  {
    EXCEPTION( "Not a valid DoF index for result '"
        << SolutionTypeEnum.ToString(solType_) << "': " << dof );
  }

  // calculate local index in dataset
  std::map<UInt, UInt>::iterator it = nodeNum2Index_.find(index_);
  if ( it == nodeNum2Index_.end() ) {
    EXCEPTION("Not a valid index in result '"
        << SolutionTypeEnum.ToString(solType_) << "': " << index_);
  } else {
    locIndex = it->second;    
  }
  
  if (entryType_ == BaseMatrix::COMPLEX) {
#ifndef NDEBUG
    if (locIndex*numDofs_ + dof - 1 >= resCacheComplex_.GetSize()) {
      EXCEPTION("Index of '" << SolutionTypeEnum.ToString(solType_) << "' from input file (id="
          << readerId_ << ") out of bounds.");
    }
#endif

    switch (outType_) {
    case OUT_REAL:
      return resCacheComplex_[locIndex*numDofs_ + dof -1].real();
    case OUT_IMAG:
      return resCacheComplex_[locIndex*numDofs_ + dof -1].imag();
    case OUT_AMPL:
      return std::abs(resCacheComplex_[locIndex*numDofs_ + dof -1]);
    case OUT_PHASE:
      // we need to convert radians to degrees, because this goes through
      // the math parser!!!
      return std::arg(resCacheComplex_[locIndex*numDofs_ + dof -1])*180.0/PI;
    }
  }
  else {
#ifndef NDEBUG
     if (locIndex*numDofs_ + dof - 1 >= resCacheDouble_.GetSize()) {
      EXCEPTION("Index of '" << SolutionTypeEnum.ToString(solType_) << "' from input file (id="
          << readerId_ << ") out of bounds.");
    }
#endif

    return resCacheDouble_[locIndex*numDofs_ + dof - 1];
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
    EXCEPTION("Result '" << SolutionTypeEnum.ToString(solType_) << "' not found in input file (id="
        << readerId << ")");
  }

  numDofs_ = resInfos[iRes]->dofNames.GetSize();
  inputReader->GetStepValues(sequenceStep, resInfos[iRes], stepValues_);

  // rebuild map node number => index
  StdVector<UInt> entityNums;
  
  switch ( resInfos[iRes]->definedOn )
  {
    case ResultInfo::NODE:
      domain->GetGrid("default")->GetNodesByName( entityNums, entityName_ );
      break;
    case ResultInfo::ELEMENT:
    case ResultInfo::SURF_ELEM:
      domain->GetGrid("default")->GetElemNumsByName( entityNums, entityName_ );
      break;
    default:
      std::string defName;
      ResultInfo::Enum2String( resInfos[iRes]->definedOn, defName );
      EXCEPTION("Result '" << SolutionTypeEnum.ToString(solType_)
                << "' is defined on " << defName
                << "s, which cannot be handled by the input function");
  }
  
  nodeNum2Index_.clear();
  for ( UInt i=0, n=entityNums.GetSize(); i<n; ++i )
  {
    nodeNum2Index_[entityNums[i]] = i;
  }

}

} // namespace CoupledField
