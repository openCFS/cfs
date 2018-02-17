// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseFilter.cc
 *       \brief    <Description>
 *
 *       \date     Aug 28, 2015
 *       \author   ahueppe
 */
//================================================================================================



#include "BaseFilter.hh"
#include "Filters/Arithmetic/BinOpFilter.hh"
#include "Filters/Input/InputFilter.hh"
#include "Filters/Output/OutputFilter.hh"
#include "Filters/Derivatives/RotatingSubstDt.hh"
#include "Filters/Derivatives/TimeDerivFilter.hh"
#include "Filters/Derivatives/PostLighthillSource.hh"
#include "SignalProcessing/FftFilter.hh"
#include "SignalProcessing/TimeMeanFilter.hh"
#include "SignalProcessing/TemporalBlendingFilter.hh"
#include <boost/tokenizer.hpp>
#include <Filters/BaseMeshFilterType.hh>

namespace CFSDat{

FilterPtr BaseFilter::Generate(PtrParamNode filtNode,PtrResultManager resMana) {
  FilterPtr newPtr;
  if (filtNode->GetName() == "meshInput") {
    newPtr = FilterPtr(new CFSDat::InputFilter(0,filtNode,resMana));
  }
  else if (filtNode->GetName() == "meshOutput") {
    newPtr = FilterPtr(new CFSDat::OutputFilter(0,filtNode,resMana));
  }
  else if (filtNode->GetName() == "timeDeriv1") {
    newPtr = FilterPtr(new CFSDat::TimeDerivFilterD1(0,filtNode,resMana));
  }
  else if (filtNode->GetName() == "substantialDeriv") {
    newPtr = FilterPtr(new CFSDat::RotatingSubstDt(0,filtNode,resMana));
  }
  else if (filtNode->GetName() == "postLighthill") {
    newPtr = FilterPtr(new CFSDat::PostLighthillSource(0, filtNode, resMana));
  }
  else if (filtNode->GetName() == "syntheticTurbulence_SNGR") {
//    newPtr = FilterPtr(new CFSDat::SNGRFilter(0, filtNode, resMana));
  }
  else if ((filtNode->GetName() == "interpolation") ||
           (filtNode->GetName() == "differentiation") ||
           (filtNode->GetName() == "aeroacoustic") ) {
    newPtr = BaseMeshFilterType::Generate(filtNode,resMana);
  }
  else if (filtNode->GetName() == "binaryOperation") {
    newPtr = BinOpFilter::GenerateOperator(filtNode,resMana);
  }
  else if (filtNode->GetName() == "fft") {
    newPtr.reset(new FftFilter(0, filtNode, resMana));
  }
  else if (filtNode->GetName() == "temporalBlending") {
      newPtr.reset(new CFSDat::TemporalBlendFilter(0, filtNode, resMana));
    }
  else if (filtNode->GetName() == "timeMean") {
    newPtr = FilterPtr(new CFSDat::TimeMeanFilter1(0,filtNode,resMana));
  }
  return newPtr;
}


BaseFilter::BaseFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
            : numWorkers_(numWorkers),
              params_(config){

  filterId_ = config->Get("id", CF::ParamNode::EX)->As<std::string>();
  
  std::cout << "\t---> Creating Filter with ID \"" << filterId_ << "\"" << std::endl;

  filtStreamType_ = BaseFilter::NO_STREAM;
  filterTag_ = boost::uuids::random_generator()();

  std::string inIds = config->Get("inputFilterIds", CF::ParamNode::EX)->As<std::string>();

  if (inIds.length() > 0 && inIds != "none") {
    typedef  boost::char_separator<char>  c_sep;
    typedef  boost::tokenizer< c_sep > tok;
  
    c_sep sep(" ,");
    tok varTok(inIds, sep);
  
    for(tok::iterator beg=varTok.begin(); beg!=varTok.end();++beg){
      inputIds_.insert(beg->c_str());
    }
  }

  resultManager_ = resMan;

  //determine global step value map
  PtrParamNode stepNode = params_->GetParent()->Get("stepValueDefinition");
  if(stepNode->Has("startStop")){
    PtrParamNode stNode = stepNode->Get("startStop");
    UInt start    = stNode->Get("startStep")->Get("value")->As<UInt>();
    UInt numSteps = stNode->Get("numSteps")->Get("value")->As<UInt>();
    Double starttime    = stNode->Get("startTime")->Get("value")->As<Double>();
    Double delta    = stNode->Get("delta")->Get("value")->As<Double>();
    string rOffset = stNode->Get("deleteOffset")->Get("value")->As<string>();
    //first step is always 1...
    if(rOffset == "yes")
      startTime_ = starttime - delta;
    for(UInt i=0;i<numSteps;++i){
        globalStepValueMap_[start+i+1] = (starttime-startTime_) + (start+i)*delta;
    }
  }
  // initializing fields for initialization of results
  initSourceResults_ = 0;
  initSinkResults_ = 0;
}


void BaseFilter::InitResults(){
  if (filtStreamType_ != OUTPUT_FILTER) {
    EXCEPTION("Function InitResults() should only be called for output filters");
  }
  std::cout << "\t---> Initializing Output Filter id " << this->filterId_ << std::endl;
  InitResultsUpstream();
}

  
void BaseFilter::InitResultsUpstream(){
  std::cout << "\t---> Initializing Filter id " << this->filterId_ << " Upstream " << std::endl;
  // search for results provided by filter itself
  ExtractFilterResults();
  
  UInt numSinks = sinks_.GetSize();
  if (numSinks > 0) {
    initSinkResults_++;
    if (initSinkResults_ < numSinks) {
      // here we wait for all requested results to be initialized
      return;
    } else if (initSinkResults_ > numSinks) {
      EXCEPTION("Function InitResultsUpstream() called too often, should not exceed number of downstream filters, maybe a circular dependency occurred");
    }
  }
  //now as all upstream filters have arrived, we go further upstream
  
  //check if filter is needed
  if (filtStreamType_ != OUTPUT_FILTER) {
    if (!CheckNeeded()) {
      return;
    }
  }
  
  //obtain results of upstream data
  //and give the implementing filter the possiblity to modify
  upResIds = this->SetUpstreamResults();

  if (upResIds.GetSize() == 0) {
    //now we arrived at an upstream end of the pipeline
    //go downstream again
    InitResultsDownstream();
    return;
  }
  
  // deactivate remaining results needed from downstream
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    resultManager_->DeactivateResult(*aIter);
  }
  // activate upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->ActivateResult(upResIds[aRes]);
  }
  // Going further upstream
  for(UInt aSrc=0;aSrc< sources_.GetSize(); aSrc++){
    sources_[aSrc]->InitResultsUpstream();
  }
  // reactivate remaining results needed from downstream
  aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    resultManager_->ActivateResult(*aIter);
  }
}


void BaseFilter::ExtractFilterResults(){
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    uuids::uuid aId = *aIter;
    ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(aId);
    //print_ConstExtInfoFields((*aInfo.get()));

    if(filtResNames.find(aInfo->resultName) != filtResNames.end()){
      //only push back the result id once.
      if(filterResIds.Find(aId) == -1){
        for (UInt i = 0; i < filterResIds.GetSize(); i++) {
          uuids::uuid otherId = filterResIds[i];
          ResultManager::ConstInfoPtr oInfo = resultManager_->GetExtInfo(otherId);
          if (aInfo->resultName == oInfo->resultName) {
            resultManager_->CombineResults(aId,otherId);
            break;
          }
        }
        filterResIds.Push_back(aId);
      }
      resultManager_->DeactivateResult(aId);
    }
  }
}


bool BaseFilter::CheckNeeded() {
  if (filterResIds.GetSize() == 0){
    std::cerr << "WARNING: filter " << filterId_ << " is not needed for the requested results. Is this what you want?" << std::endl;
    return false;
  }
  return true;
}

ResultIdList BaseFilter::SetDefaultUpstreamResults() {
  ResultIdList resultIds;
  std::set<std::string>::iterator inItr = upResNames.begin();
  for (; inItr != upResNames.end(); inItr++) {
    std::string name = *inItr;
    uuids::uuid newId = RegisterUpstreamResult(name, filterResIds[0]);
    resultIds.Push_back(newId);
  }
  return resultIds;
}


uuids::uuid BaseFilter::RegisterUpstreamResult(std::string name, uuids::uuid downStreamResultId) {
  return RegisterUpstreamResult(name, 0, 0, downStreamResultId);
}


uuids::uuid BaseFilter::RegisterUpstreamResult(std::string name, Integer minOffset, 
                         Integer maxOffset, uuids::uuid downStreamResultId) {
  Integer downStreamMaxStepOffset = 0;
  if (downStreamResultId != uuids::nil_uuid()) {
    if (filterResIds.Find(downStreamResultId) == -1) {
      EXCEPTION("RegisterUpstreamResult not called with valid");
    }
    // here we get the offset of the newest downstream result to be evaluated
    ResultManager::ConstInfoPtr outInPtr = resultManager_->GetExtInfo(downStreamResultId);
    downStreamMaxStepOffset = outInPtr->maxStepOffset;
  }
  std::cout << " register " << name << std::endl;
  uuids::uuid newId = resultManager_->AddResult(name,this->filterTag_,
                            downStreamMaxStepOffset+minOffset,
                            downStreamMaxStepOffset+maxOffset);
  if (downStreamResultId != uuids::nil_uuid()) {
    resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(downStreamResultId)->timeLine.get()));
  }
  upResNameIds[name] = newId;
  return newId;
}


void BaseFilter::InitResultsDownstream(){
  std::cout << "\t---> Initializing Filter id " << this->filterId_ << " Downstream " << std::endl;
  
  UInt numSources = sources_.GetSize();
  if (numSources > 0) {
    initSourceResults_++;
    if (initSourceResults_ < numSources) {
      return;
    } else if (initSourceResults_ > numSources) {
      EXCEPTION("Function InitResultsDownstream() called too often, should not exceed number of upstream filters");
    }
  }
  //now as all upstream filters have arrived, we go further downstream

  //adapt the filter results according to the upstream data
  AdaptFilterResults();
  
  //loop over sinks
  for(UInt aSin=0;aSin< sinks_.GetSize(); aSin++){
    sinks_[aSin]->InitResultsDownstream();
  }
}

bool BaseFilter::Run() {
  std::set<uuids::uuid> obsoleteResults = ExtractObsoleteResults();
  if (obsoleteResults.empty()) {
    return true;
  }
  bool success = this->UpdateResults(obsoleteResults);
  DeactivateUpstreamResults();
  if (success) {
    std::set<uuids::uuid>::iterator sIter = obsoleteResults.begin();
    for (; sIter != obsoleteResults.end(); ++sIter) {
      resultManager_->SetResultVecUpToDate(*sIter,true);
    }
  } else {
    WARN("Some results needed for " << this->filterId_ << " could not be calculated");
  }
  return success;
}
  
bool BaseFilter::UpdateResults(std::set<uuids::uuid>& upResults) {
  EXCEPTION("Function UpdateResults() called but not implemented for filter");
  return false;
}

std::set<uuids::uuid> BaseFilter::ExtractObsoleteResults() {
  std::set<uuids::uuid> requestedResults;
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator actRes = activeResults.begin();
  for (;actRes != activeResults.end(); actRes++) {
    if(filterResIds.Find(*actRes) > -1) {
      if (!resultManager_->IsResultVecUpToDate(*actRes)) {
        requestedResults.insert(resultManager_->GetMasterResult(*actRes));
      }
    }
  }
  return requestedResults;
}


void BaseFilter::DeactivateUpstreamResults() {
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }
}

void BaseFilter::PrepareUpstreamResult(uuids::uuid resId) {
  if (!resultManager_->IsResultVecUpToDate(resId)) {
    resultManager_->ActivateResult(resId);
    CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
    for(; srcIter != sources_.End() ; srcIter++){
      // should we check here anything for success?
      (*srcIter)->Run();
    }
    if (!resultManager_->IsResultVecUpToDate(resId)) {
      WARN("Receiving upstream result " << resultManager_->GetResultName(resId) << " did not work");
    }
  }
}


#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template CF::Vector<Double>& BaseFilter::GetOwnResultVector(uuids::uuid, CF::StdVector<UInt>&);
  template CF::Vector<Complex>& BaseFilter::GetOwnResultVector(uuids::uuid, CF::StdVector<UInt>&);
  template CF::Vector<Double>& BaseFilter::GetOwnResultVector(uuids::uuid);
  template CF::Vector<Complex>& BaseFilter::GetOwnResultVector(uuids::uuid);

  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(uuids::uuid, Double, CF::StdVector<UInt>&);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(uuids::uuid, Double, CF::StdVector<UInt>&);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(uuids::uuid, CF::StdVector<UInt>&);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(uuids::uuid, CF::StdVector<UInt>&);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(uuids::uuid, Double);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(uuids::uuid, Double);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(uuids::uuid);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(uuids::uuid);
  
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(std::string, Double, CF::StdVector<UInt>&);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(std::string, Double, CF::StdVector<UInt>&);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(std::string, CF::StdVector<UInt>&);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(std::string, CF::StdVector<UInt>&);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(std::string, Double);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(std::string, Double);
  template CF::Vector<Double>& BaseFilter::GetUpstreamResultVector(std::string);
  template CF::Vector<Complex>& BaseFilter::GetUpstreamResultVector(std::string);
#endif
}
