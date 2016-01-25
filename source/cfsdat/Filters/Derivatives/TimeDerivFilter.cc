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
 *       \author   ahueppe
 */
//================================================================================================


#include "TimeDerivFilter.hh"

namespace CFSDat{

TimeDerivFilter::TimeDerivMap TimeDerivFilter::tDerivMap_ = TimeDerivFilter::init_deriv_map();

TimeDerivFilterD1::TimeDerivFilterD1(UInt numWorkers, CoupledField::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
: TimeDerivFilter(numWorkers,config,resMan){
 //ok, so either the user specifies
 //results to be differentiated or we choose default possibilities from our map
 // for now this is the way to go
  this->filtSteamType_ = FIFO_FILTER;
  //first case we have explicitly given the single result tag
  if( params_->Has("singleResult") ){
    PtrParamNode singleNode = params_->Get("singleResult");
    std::string iRes = singleNode->Get("inputQuantity")->Get("resultName")->As<std::string>();
    std::string oRes = singleNode->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(oRes);
    upResNames.insert(iRes);
    InOutBiMap::value_type pair(oRes,iRes);
    inOutNames_.insert(pair);

  }else if (params_->Has("allResults") ){
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

bool TimeDerivFilterD1::Run(){
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1)
      continue;
    //we now set the time for the associated results
    Double aTF = resultManager_->GetStepValue(*aIter);
    std::string resName = resultManager_->GetExtInfo(*aIter)->resultName;
    //do we already have the timestep? is it constant?
    //yes for now
    uuids::uuid upRes = filtResToUpResIds_[resName];

    //here is the trick, only the obsolete result gets activated
    // for the rest we just interchange filtResToUpResIds_[resName]
    // for this second order derivative we cannot do anything
    // but if we go to higher order filteres we try to pass only
    // results which are new to upstream
    resultManager_->SetTimeValue(upRes,aTF);
    // now we deactivate our own result and activate the others
    resultManager_->ActivateResult(upRes);

    resultManager_->DeactivateResult(*aIter);
  }

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }
  aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1)
      continue;
    std::string resName = resultManager_->GetExtInfo(*aIter)->resultName;
    uuids::uuid upRes = filtResToUpResIds_[resName];

    CF::StdVector<UInt> eqnNums; //they should be equal...
    Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(*aIter,eqnNums);
    Vector<Double>& r1 = resultManager_->GetResultVector<Double>(upRes,eqnNums,-2);
    Vector<Double>& r2 = resultManager_->GetResultVector<Double>(upRes,eqnNums,-1);
    Vector<Double>& r3 = resultManager_->GetResultVector<Double>(upRes,eqnNums,1);
    Vector<Double>& r4 = resultManager_->GetResultVector<Double>(upRes,eqnNums,2);
    UInt last = (eqnNums.GetSize() == 0)? returnVec.GetSize() : eqnNums.GetSize();


    for(UInt i=0;i<last;++i){
      returnVec[i] = 2*(((r3[i]-r2[i])+r4[i]-r1[i])/(8*timeSteps_[*aIter]));
    }
    resultManager_->ActivateResult(*aIter);
  }

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }
  return true;
}


ResultIdList TimeDerivFilterD1::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
    std::string upstreamRes = inOutNames_.left.at(filterResName);
    CF::StdVector<Integer> timeLine(4);
    timeLine[0] = -2;
    timeLine[1] = -1;
    timeLine[2] = 1;
    timeLine[3] = 2;
    uuids::uuid newId = resultManager_->AddResult(upstreamRes,this->filterTag_,timeLine);
    resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
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
    }
    //determine the timestep by subtracting first and second entry in timeline
    timeSteps_[filterResIds[aRes]] = std::abs( curTime[1] - curTime[0]);

    resultManager_->CopyResultData(assocId,filterResIds[aRes]);
    resultManager_->SetValid(filterResIds[aRes]);
  }

}
}

