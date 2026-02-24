// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DerivFilterD1.cc
 *       \brief    <Description>
 *
 *       \date     Nov 4, 2015
 *       \author   ahueppe, sschoder
 */
//================================================================================================
//http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/

#include "TimeDerivFilter.hh"

namespace CFSDat{

TimeDerivFilter::TimeDerivMap TimeDerivFilter::tDerivMap_ = TimeDerivFilter::init_deriv_map();

TimeDerivFilterD1::TimeDerivFilterD1(UInt numWorkers, CoupledField::PtrParamNode config, shared_ptr<ResultManager> resMan)
: TimeDerivFilter(numWorkers,config,resMan){
 //ok, so either the user specifies
 //results to be differentiated or we choose default possibilities from our map
 // for now this is the way to go

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  //first case we have explicitly given the single result tag
  if( params_->Has("singleResult") ){
    PtrParamNode singleNode = params_->Get("singleResult");
    std::string iRes = singleNode->Get("inputQuantity")->Get("resultName")->As<std::string>();
    std::string oRes = singleNode->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(oRes);
    upResNames.insert(iRes);
    InOutBiMap::value_type pair(oRes,iRes);
    inOutNames_.insert(pair);
  } //take standard derivatives input and output defined in TimeDerivMap (header file)
  else if (params_->Has("allResults") ){
    TimeDerivFilter::TimeDerivMap::iterator resIt = TimeDerivFilter::tDerivMap_.begin();
    for(;resIt != TimeDerivFilter::tDerivMap_.end();resIt++){
      std::string myName = CF::SolutionTypeEnum.ToString(resIt->second.d1Type);
      std::string upName = CF::SolutionTypeEnum.ToString(resIt->first);
      filtResNames.insert(myName);
      upResNames.insert(upName);
      //in this filter we have a 1-1 connection so we can use bimap
      //other filters may need different things
      InOutBiMap::value_type pair(upName,myName);
      inOutNames_.insert(pair);
    }
  }
}

bool TimeDerivFilterD1::UpdateResults(std::set<uuids::uuid>& upResults) {
  std::set<uuids::uuid>::iterator oIter = upResults.begin();
  for(; oIter != upResults.end(); ++oIter){
    std::string resName = resultManager_->GetExtInfo(*oIter)->resultName;
    uuids::uuid upRes = filtResToUpResIds_[resName];

    Vector<Double>& returnVec = GetOwnResultVector<Double>(*oIter);
    const UInt size = returnVec.GetSize();
    const Double dT = timeSteps_[*oIter];
    const Integer actStepIndex = resultManager_->GetStepIndex(*oIter);
    
    if (dT == 0.0) {
      // temporal derivative of a constant value is zero
      returnVec.Init(0.0);
    } else {
      Vector<Double>& inVec = GetUpstreamResultVector<Double>(upRes, actStepIndex - 2);
      Double factor = -(1.0/8.0)/dT;
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt i = 0; i < size; i++) {
        returnVec[i] = inVec[i] * factor;
      }
  
      inVec = GetUpstreamResultVector<Double>(upRes, actStepIndex - 1);
      factor = -(1.0/4.0)/dT;
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt i = 0; i < size; i++) {
        returnVec[i] += inVec[i] * factor;
      }
  
      inVec = GetUpstreamResultVector<Double>(upRes, actStepIndex + 1);
      factor = (1.0/4.0)/dT;
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt i = 0; i < size; i++) {
        returnVec[i] += inVec[i] * factor;
      }
    
      inVec = GetUpstreamResultVector<Double>(upRes, actStepIndex + 2);
      factor = (1.0/8.0)/dT;
      #pragma omp parallel for num_threads(CFS_NUM_THREADS)
      for (UInt i = 0; i < size; i++) {
        returnVec[i] += inVec[i] * factor;
      }
    }
  }
  return true;
}


ResultIdList TimeDerivFilterD1::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
    std::string upstreamRes = inOutNames_.left.at(filterResName);
    uuids::uuid newId = RegisterUpstreamResult(upstreamRes, -2, 2, *aIt);
    generated.Push_back(newId);
    //additionally we store the uuids belonging to one upRes
    filtResToUpResIds_[filterResName] = newId;
  }
  return generated;
}

void TimeDerivFilterD1::AdaptFilterResults(){
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

    CF::StdVector<Double> curTime = *inInfo->timeLine.get();
    if(curTime.GetSize() == 0){
      std::stringstream s;
      s << "No upstream filter provided time information for the result "  << aFiltResName;
      CF::Exception(s.str());
    } else if (curTime.GetSize() == 1) {
      timeSteps_[filterResIds[aRes]] = 0.0;
    } else {
      //determine the timestep by subtracting first and second entry in timeline
      timeSteps_[filterResIds[aRes]] = std::abs( curTime[1] - curTime[0]);
    }

    resultManager_->CopyResultData(assocId,filterResIds[aRes]);
    resultManager_->SetValid(filterResIds[aRes]);
  }

}
}

