// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FIRFilter.cc
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2018
 *       \author   mtautz
 */
//================================================================================================


#include "FIRFilter.hh"

#include <cstdlib>
#include <boost/tokenizer.hpp>

namespace CFSDat{

FIRFilter::FIRFilter(UInt numWorkers, CoupledField::PtrParamNode config, shared_ptr<ResultManager> resMan)
: BaseFilter(numWorkers,config,resMan){

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  
  if( params_->Has("singleResult") ){
    PtrParamNode singleNode = params_->Get("singleResult");
    std::string iRes = singleNode->Get("inputQuantity")->Get("resultName")->As<std::string>();
    std::string oRes = singleNode->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(oRes);
    upResNames.insert(iRes);
    InOutBiMap::value_type pair(oRes,iRes);
    inOutNames_.insert(pair);
  } 
  
  // checking for up- and downsampling
  upsampling_ = false;
  downsampling_ = false;
  Integer resampleFactor = 1;
  if( params_->Has("upsampling") ) {
    upsampling_ = true;
    resampleFactor = params_->Get("upsampling")->As<Integer>();
  }
  if( params_->Has("downsampling") ) {
    if (upsampling_) {
      EXCEPTION("Can not perform upsampling and downsampling at the same time");
    }
    downsampling_ = true;
    resampleFactor = params_->Get("downsampling")->As<Integer>();
  }
  if (upsampling_ || downsampling_) {
    if (resampleFactor == 0) {
      EXCEPTION("Resampling factor can not be zero");
    }
    if (resampleFactor < 0) {
      resampleFactor = -resampleFactor;
    }
    if (resampleFactor == 1) {
      upsampling_ = false;
      downsampling_ = false;
    }
    resampleFactor_ = resampleFactor;
  }
  
  // loading coefficients from nodes or using default values
  if( params_->Has("coefficients") ) {
    std::string coefString = params_->Get("coefficients")->As<std::string>();
    typedef  boost::char_separator<char>  c_sep;
    typedef  boost::tokenizer< c_sep > tok;
    c_sep sep(" ,");
    tok varTok(coefString, sep);
    Integer index = 0;
    for(tok::iterator beg=varTok.begin(); beg!=varTok.end();++beg){
      Double coef = atof(beg->c_str());
      if (coef != 0.0 && coef != -0.0) {
        coefficients_.Push_back(coef);
        coefficientOffsets_.Push_back(index);
      }
      index++;
    }
    CoupledField::PtrParamNode coefNode = params_->Get("coefficients");
    bool center = false;
    Integer offsetMinus = 0;
    if (coefNode->Has("offset")) {
      std::string offsetString = coefNode->Get("offset")->As<std::string>();
      if (!offsetString.compare("forward")) {
        offsetMinus = 0;
      } else if (!offsetString.compare("backward")) {
        offsetMinus = index - 1;
      } else if (!offsetString.compare("center")) {
        center = true;
      } else {
        offsetMinus = atoi(offsetString.c_str());
      }
    } else {
      center = true;
    }
    if (center) {
      offsetMinus = index / 2;
    }
    for (UInt iCoef = 0; iCoef < coefficients_.GetSize(); iCoef++) {
      coefficientOffsets_[iCoef] -= offsetMinus;
    }
  } else {
    SetDefaultCoefficients();
  }
  minOffset_ = coefficientOffsets_[0];
  maxOffset_ = minOffset_;
  for (UInt iCoef = 1; iCoef < coefficients_.GetSize(); iCoef++) {
    minOffset_ = std::min(minOffset_,coefficientOffsets_[iCoef]);
    maxOffset_ = std::max(maxOffset_,coefficientOffsets_[iCoef]);
  }
  if (upsampling_) {
    // multiplication for uspampling to preserve energy content
    Double factor = (Double)resampleFactor_;
    for (UInt iCoef = 0; iCoef < coefficients_.GetSize(); iCoef++) {
      coefficients_[iCoef] *= factor;
    }
  }
  /**
  for (UInt iCoef = 0; iCoef < coefficients_.GetSize(); iCoef++) {
    std::cout << "coefIdx " << coefficientOffsets_[iCoef] << "   coef " << coefficients_[iCoef] << std::endl;
  }
  **/
}

bool FIRFilter::UpdateResults(std::set<uuids::uuid>& upResults) {
  std::set<uuids::uuid>::iterator oIter = upResults.begin();
  for(; oIter != upResults.end(); ++oIter){
    std::string resName = resultManager_->GetExtInfo(*oIter)->resultName;
    uuids::uuid upRes = filtResToUpResIds_[resName];
    Vector<Double>& returnVec = GetOwnResultVector<Double>(*oIter);
    
    Integer actStepIndex = resultManager_->GetStepIndex(*oIter);
    if (upsampling_) {
      Integer downActStepIndex = actStepIndex / resampleFactor_;
      // now we select all steps from upstream, which we really need
      // some steps are omitted due to zero padding
      StdVector<Double> selectCoefficients;
      StdVector<Integer> selectCoefficientOffsets; 
      for (UInt iCoeff = 0; iCoeff < coefficients_.GetSize(); iCoeff++) {
        Integer iCoeffStep = actStepIndex + coefficientOffsets_[iCoeff];
        while (iCoeffStep < 0) {
          iCoeffStep += resampleFactor_;
        }
        if (iCoeffStep % resampleFactor_ == 0) {
          selectCoefficients.Push_back(coefficients_[iCoeff]);
          selectCoefficientOffsets.Push_back(iCoeffStep / resampleFactor_ - downActStepIndex);
        }
      }
      ApplyFIRFilter(returnVec, upRes, downActStepIndex, selectCoefficients, selectCoefficientOffsets);
    } else if (downsampling_) {
      actStepIndex *= resampleFactor_;
      ApplyFIRFilter(returnVec, upRes, actStepIndex, coefficients_, coefficientOffsets_);
    } else {
      ApplyFIRFilter(returnVec, upRes, actStepIndex, coefficients_, coefficientOffsets_);
    }
  }
  return true;
}

void FIRFilter::ApplyFIRFilter(Vector<Double>& returnVec, 
       uuids::uuid upRes, Integer actStepIndex,
       StdVector<Double>& coefficients, StdVector<Integer>& coefficientOffsets) {
  returnVec.Init(0.0);
  const UInt size = returnVec.GetSize();
  for (UInt iCoeff = 0; iCoeff < coefficients.GetSize(); iCoeff++) {
    const Double factor = coefficients[iCoeff];
    Vector<Double>& inVec = GetUpstreamResultVector<Double>(upRes, 
              actStepIndex + coefficientOffsets[iCoeff]);
    #pragma omp parallel for num_threads(CFS_NUM_THREADS) 
    for (UInt i = 0; i < size; i++) {
      returnVec[i] += inVec[i] * factor;
    }
  }
}

ResultIdList FIRFilter::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
    std::string upstreamRes = inOutNames_.left.at(filterResName);
    uuids::uuid newId;
    if (upsampling_) {
      Integer offsetDistance = maxOffset_ - minOffset_;
      Integer steps = offsetDistance / resampleFactor_ 
                      + (offsetDistance % resampleFactor_ > 0 ? 1 : 0);
      Integer maxOffset = maxOffset_ / resampleFactor_ 
                      + (maxOffset_ % resampleFactor_ > 0 ? 1 : 0);
      Integer minOffset = maxOffset - steps;
      newId = RegisterUpstreamResult(upstreamRes, minOffset, maxOffset, *aIt);
    } else {
      newId = RegisterUpstreamResult(upstreamRes, minOffset_, maxOffset_, *aIt);
    }
    generated.Push_back(newId);
    //additionally we store the uuids belonging to one upRes
    filtResToUpResIds_[filterResName] = newId;
  }
  return generated;
}

Integer FIRFilter::GetDownStreamMaxStepOffset(std::set<uuids::uuid> downStreamResultIds) {
  Integer downStreamMaxStepOffset = BaseFilter::GetDownStreamMaxStepOffset(downStreamResultIds);
  if (upsampling_) {
    downStreamMaxStepOffset /= resampleFactor_;
  } else if (downsampling_) {
    downStreamMaxStepOffset *= resampleFactor_;
  }
  return downStreamMaxStepOffset;
}

void FIRFilter::CopyTimeLineUpstream(uuids::uuid upStreamId, uuids::uuid downStreamId) {
  if (resultManager_->IsStatic(downStreamId)) {
  //if (resultManager_->IsConstant(downStreamId)) {
    BaseFilter::CopyTimeLineUpstream(upStreamId, downStreamId);
    return;
  }
  if (upsampling_) {
    StdVector<Double> dTimeLine = resultManager_->GetTimeLine(downStreamId);
    StdVector<Double> uTimeLine;
    DownsampleTimeLine(dTimeLine, uTimeLine);
    resultManager_->SetTimeLine(upStreamId, uTimeLine);
  } else if (downsampling_) {
    StdVector<Double> dTimeLine = resultManager_->GetTimeLine(downStreamId);
    StdVector<Double> uTimeLine;
    UpsampleTimeLine(dTimeLine, uTimeLine);
    resultManager_->SetTimeLine(upStreamId, uTimeLine);
  } else {
    BaseFilter::CopyTimeLineUpstream(upStreamId, downStreamId);
  }
}

void FIRFilter::UpsampleTimeLine(StdVector<Double>& inputTimeLine, StdVector<Double>& upSampleTimeLine) {
  // linear interpolation of time step values
  StdVector<Double> FactorsA;
  StdVector<Double> FactorsB;
  for (UInt j = 0; j < resampleFactor_; j++) {
    FactorsA.Push_back(Double(resampleFactor_-j) / Double(resampleFactor_));
    FactorsB.Push_back(Double(j) / Double(resampleFactor_));
  }
  Double addValue = inputTimeLine[inputTimeLine.GetSize() - 1];
  addValue += (inputTimeLine[1] - inputTimeLine[0]);
  inputTimeLine.Push_back(addValue);
  upSampleTimeLine.Clear();
  for (UInt i = 0; i < inputTimeLine.GetSize() - 1; i++) {
    for (UInt j = 0; j < resampleFactor_; j++) {
      upSampleTimeLine.Push_back(inputTimeLine[i] * FactorsA[j] + inputTimeLine[i+1] * FactorsB[j]);
    }
  }
}

void FIRFilter::DownsampleTimeLine(StdVector<Double>& inputTimeLine, StdVector<Double>& downSampleTimeLine) {
  downSampleTimeLine.Clear();
  for (UInt i = 0; i < inputTimeLine.GetSize(); i += resampleFactor_) {
    downSampleTimeLine.Push_back(inputTimeLine[i]);
  }
}

void FIRFilter::AdaptFilterResults(){
  //we would need to copy related data to our filter results
  //in our case we obtain the information from the upstream
  // we set our own results
  for(UInt aRes = 0; aRes < filterResIds.GetSize(); aRes++){
    ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(filterResIds[aRes]);
    std::string aFiltResName = aInfo->resultName;

    uuids::uuid assocId = filtResToUpResIds_[aFiltResName];
    ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(assocId);

    //got the upstream result validated?
    if(!inInfo->isValid){
      EXCEPTION("Problem in filter pipeline detected. Time derivative input result \"" <<  inInfo->resultName << "\" could not be provided.")
    }

    resultManager_->CopyResultData(assocId,filterResIds[aRes]);
    if (!resultManager_->IsStatic(assocId) && (upsampling_ || downsampling_)) {
    //if (!resultManager_->IsConstant(assocId) && (upsampling_ || downsampling_)) {
      StdVector<Double> uTimeLine = resultManager_->GetTimeLine(assocId);
      StdVector<Double> dTimeLine;
      if (upsampling_) {
        UpsampleTimeLine(uTimeLine, dTimeLine);
      } else {
        DownsampleTimeLine(uTimeLine, dTimeLine);
      }
      resultManager_->SetTimeLine(filterResIds[aRes], dTimeLine);
    }
    resultManager_->SetValid(filterResIds[aRes]);
  }

}


void FIRFilter::SetDefaultCoefficients() {
  coefficients_.Resize(1);
  coefficientOffsets_.Resize(1);
  coefficients_[0] = 1.0;
  coefficientOffsets_[0] = 0;
}

}

