// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TimeMeanFilter.cc
 *       \brief    <Description>
 *
 *       \date     Apr 1, 2017
 *       \author   sschoder
 */
//================================================================================================


#include <Filters/SignalProcessing/TimeMeanFilter.hh>

namespace CFSDat{

TimeMeanFilter::TimeMeanMap TimeMeanFilter::tMeanMap_ = TimeMeanFilter::init_mean_map();

TimeMeanFilter1::TimeMeanFilter1(UInt numWorkers, CoupledField::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
: TimeMeanFilter(numWorkers,config,resMan){
 //ok, so either the user specifies
 //results to be differentiated or we choose default possibilities from our map
 // for now this is the way to go

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  //first case we have explicitly given the single result tag
  if(params_->Has("meanSpan")){
    N_ = params_->Get("meanSpan")->As<int>();
  } else {
    WARN("No time meanSpan set for filter " << filterId_ << ", setting 20 now");
    N_ = 20;
  }
  if (N_ == 0) {
    WARN("meanSpan set to 0 for filter " << filterId_ << ", should at least be one, setting 1 now");
    N_ = 1;
  }
  if( params_->Has("singleResult") ){
    PtrParamNode singleNode = params_->Get("singleResult");
    std::string iRes = singleNode->Get("inputQuantity")->Get("resultName")->As<std::string>();
    std::string oRes = singleNode->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(oRes);
    upResNames.insert(iRes);
    InOutBiMap::value_type pair(oRes,iRes);
    inOutNames_.insert(pair);
  }
  //TODO FIR FILTER SETUP general
}

bool TimeMeanFilter1::UpdateResults(std::set<uuids::uuid>& upResults) {
  std::set<uuids::uuid>::iterator oIter = upResults.begin();
  //std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  //std::set<uuids::uuid>::iterator aIter = activeResults.begin();
  for(; oIter != upResults.end(); ++oIter){
    //we now set the time for the associated results
    std::string resName = resultManager_->GetExtInfo(*oIter)->resultName;
    //do we already have the timestep? is it constant?
    //yes for now
    uuids::uuid upRes = filtResToUpResIds_[resName];
    
    Vector<Double>& returnVec = GetOwnResultVector<Double>(*oIter);
    
    const UInt size = returnVec.GetSize();
    const Integer actStepIndex = resultManager_->GetStepIndex(*oIter);
    
    const Integer iStart = N_ > 0 ? -N_+1 : 0;
    const Integer iEnd = N_ > 0 ? 0 : -N_-1;
    
    returnVec.Init(0.0);
    for (Integer i = iStart; i <= iEnd; i++) {
      Vector<Double>& inVec = GetUpstreamResultVector<Double>(upRes, actStepIndex + i);
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt i = 0; i < size; i++) {
        returnVec[i] += inVec[i];
      }
    }
    
    const Double factor = std::abs(1.0 / (Double)N_);
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < size; i++) {
      returnVec[i] *= factor;
    }
  }
  
  return true;
}


ResultIdList TimeMeanFilter1::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
    std::string upstreamRes = inOutNames_.left.at(filterResName);
    uuids::uuid newId;
    if (N_ > 0) {
      newId = RegisterUpstreamResult(upstreamRes, -N_+1, 0, *aIt);
    } else {
      newId = RegisterUpstreamResult(upstreamRes, 0, -N_-1, *aIt);
    }
    generated.Push_back(newId);
    //additionally we store the uuids belonging to one upRes
    filtResToUpResIds_[filterResName] = newId;
  }
  return generated;
}

void TimeMeanFilter1::AdaptFilterResults(){
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
      EXCEPTION("Problem in filter pipeline detected. Time mean input result \"" <<  inInfo->resultName << "\" could not be provided.")
    }


    CF::StdVector<Double> curTime = *inInfo->timeLine.get();
    if(curTime.GetSize() == 0){
      std::stringstream s;
      s << "No upstream filter provided time information for the result "  << aFiltResName;
      CF::Exception(s.str());
    }
    //determine the timestep by subtracting first and second entry in timeline
    timeSteps_[filterResIds[aRes]] = std::abs( curTime[1] - curTime[0]);

    resultManager_->CopyResultData(assocId,filterResIds[aRes]);
    resultManager_->SetValid(filterResIds[aRes]);

  }

}
}

