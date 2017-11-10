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

bool TimeMeanFilter1::Run(){
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
    UInt last = (eqnNums.GetSize() == 0)? returnVec.GetSize() : eqnNums.GetSize();


    int initialTimeStep = (*(resultManager_->GetExtInfo(upRes)->timeLine.get()))[0] / timeSteps_[*aIter];
    int currentStep = resultManager_->GetStepValue(*aIter,0)/timeSteps_[*aIter] - initialTimeStep;
    if (0 == currentStep){
      tmp_ = resultManager_->GetResultVector<Double>(*aIter,eqnNums);
      tmp_.Init();
    }


    if (currentStep < N_){
      Vector<Double>& r = resultManager_->GetResultVector<Double>(upRes,eqnNums,0);

      // Calculation procedure
      for(UInt i=0;i<last;++i){
        tmp_[i] = (((currentStep)*tmp_[i] + r[i])/(currentStep + 1));
        returnVec[i] = tmp_[i];
      }
    }else {
      Vector<Double>& r1 = resultManager_->GetResultVector<Double>(upRes,eqnNums,-N_);
      Vector<Double>& r2 = resultManager_->GetResultVector<Double>(upRes,eqnNums,0);

      // Calculation procedure
      for(UInt i=0;i<last;++i){
        tmp_[i] += ((r2[i] - r1[i]  )/(N_));
        returnVec[i] = tmp_[i];
      }
    }
//    std::cout<< returnVec[0] << std::endl;
//    std::cout<< returnVec[1] << std::endl;
//    std::cout<< returnVec[2] << std::endl;
    resultManager_->ActivateResult(*aIter);
  }

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }
  return true;
}


ResultIdList TimeMeanFilter1::SetUpstreamResults(){
  ResultIdList generated;
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();

  for(;aIt!=filterResIds.End();++aIt){
    std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;
    std::string upstreamRes = inOutNames_.left.at(filterResName);

    uuids::uuid newId = resultManager_->AddResult(upstreamRes,this->filterTag_);
    resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
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

    if (N_ == 0){
      std::cout<<"\nWARNING:     \n" <<
                 "           Problem in filter input of "<<  inInfo->resultName <<" detected. Time meanSpan must be >0. \n" <<
                 "           N_ is set to 20.\n"<< std::endl;
      N_ = 20;
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

