// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     OutputFilter.cc
 *       \brief    <Description>
 *
 *       \date     Oct 26, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "OutputFilter.hh"
#include <DataInOut/DefineInOutFiles.hh>

namespace CFSDat{


OutputFilter::OutputFilter(UInt numWorkers, CoupledField::PtrParamNode config, shared_ptr<ResultManager> resMan)
             :BaseFilter(numWorkers,config,resMan){

  infoNode_ = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT));
  std::string name = "cfsDatInfoFilter-";
  name += this->filterId_;
  infoNode_->SetName(name.c_str());

  outFile_ = CoupledField::DefineInOutFiles::CreateSingleOutputFileObject(this->filterId_,params_->Get("outputFile")->GetChild(),infoNode_,false);

  this->filtStreamType_ = OUTPUT_FILTER;
}

OutputFilter::~OutputFilter(){
  outFile_->FinishMultiSequenceStep();
}

bool OutputFilter::Run(){
  /**
  ResultIdList::iterator rIter = upResIds.Begin();
  for(; rIter != upResIds.End() ; rIter++){
    resultManager_->SetStepValue(*rIter,aStepIter_->second);
    resultManager_->ActivateResult(*rIter);
    std::cout << " activate result " << std::endl;
  }

  CF::StdVector< shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }
**/
/*  // Filter Start value offset
  //determine global step value map
  PtrParamNode stepNode = params_->GetParent()->Get("stepValueDefinition");
  UInt start = 1;
  if(stepNode->Has("startStop")){
    start = stepNode->Get("startStop")->Get("startStep")->Get("value")->As<UInt>();
    if(aStepIter_->first == start){
      std::cout<<aStepIter_->first<<std::endl;
    }
  }*/

  //lets write results if they are valid
  ResultIdList::iterator rIter = upResIds.Begin();
//  rIter = upResIds.Begin();
  bool allsuccess = true;
  for(; rIter != upResIds.End() ; rIter++){
    //if(resultManager_->IsResultVecUpToDate(*rIter)){
      if(resultManager_->GetExtInfo(*rIter)->isMeshResult){
        StdVector< shared_ptr<BaseResult> > cResVec = resultManager_->GetBaseResultVector(*rIter);
        CF::StdVector<UInt> eqnVec;

        outFile_->BeginStep(aStepIter_->first,aStepIter_->second); //Add step values and time line in output file
        if(resultManager_->GetExtInfo(*rIter)->dType == ExtendedResultInfo::COMPLEX){
          Vector<Complex> & fullVec = GetUpstreamResultVector<Complex>(*rIter,aStepIter_->second,eqnVec);
          if (!resultManager_->IsResultVecUpToDate(*rIter)) {
            resultManager_->DeactivateResult(*rIter);
            std::cout << "      Error writing field " << resultManager_->GetResultName(*rIter) << " was not up to date " << std::endl; 
            continue;
          }
          std::cout << "      Writing field " << resultManager_->GetResultName(*rIter) << std::endl; 
          
          //now we loop over the result array and copy the values according to
          for(UInt aRe = 0; aRe < cResVec.GetSize(); ++aRe){

            Result<Complex>* myResult = dynamic_cast<Result<Complex>* >(cResVec[aRe].get());
            Vector<Complex> & resVec =  myResult->GetVector();
            
            if( resVec.ContainsNaN() || resVec.ContainsInf() ){
              WARN("Detected result with NAN or INF values. Setting this result to zero.");
              resVec.Init(0.0);
            }

            //obtain region equations
            std::string regName = cResVec[aRe]->GetEntityList()->GetName();
            CF::RegionIdType rId = resultManager_->GetExtInfo(*rIter)->ptGrid->GetRegion().Parse(regName);
            resultManager_->GetEqnMap(*rIter)->GetRegionEquations(eqnVec,rId);
            resVec.Resize(eqnVec.GetSize()); //TODO
            for(UInt aEq = 0; aEq<eqnVec.GetSize();++aEq){
              resVec[aEq] =  fullVec[eqnVec[aEq]];
            }
            outFile_->AddResult(cResVec[aRe]);
          }
        }else{
          Vector<Double> & fullVec = GetUpstreamResultVector<Double>(*rIter,aStepIter_->second,eqnVec);
          if (!resultManager_->IsResultVecUpToDate(*rIter)) {
            resultManager_->DeactivateResult(*rIter);
            std::cout << "      Error writing field " << resultManager_->GetResultName(*rIter) << " was not up to date " << std::endl; 
            continue;
          }
          std::cout << "      Writing field " << resultManager_->GetResultName(*rIter) << std::endl; 

          //now we loop over the result array and copy the values according to
          for(UInt aRe = 0; aRe < cResVec.GetSize(); ++aRe){

            Result<Double>* myResult = dynamic_cast<Result<Double>* >(cResVec[aRe].get());
            Vector<Double> & resVec =  myResult->GetVector();
            
            if( resVec.ContainsNaN() || resVec.ContainsInf() ){
              WARN("Detected result with NAN or INF values. Setting this result to zero.");
              resVec.Init(0.0);
            }

            //obtain region equations
            std::string regName = cResVec[aRe]->GetEntityList()->GetName();
            CF::RegionIdType rId = resultManager_->GetExtInfo(*rIter)->ptGrid->GetRegion().Parse(regName);
            resultManager_->GetEqnMap(*rIter)->GetRegionEquations(eqnVec,rId);
            resVec.Resize(eqnVec.GetSize());
            for(UInt aEq = 0; aEq<eqnVec.GetSize();++aEq){
              resVec[aEq] =  fullVec[eqnVec[aEq]];
            }
            outFile_->AddResult(cResVec[aRe]);
          }
        }

        outFile_->FinishStep();

        //now deactivate own upstream results
        DeactivateUpstreamResults();
        
        resultManager_->SetValid(*rIter);

      }else{
        allsuccess &= false;
      }

  }

  ++aStepIter_;
  //check for last step
  if(aStepIter_ == globalStepValueMap_.end())
    return allsuccess;
  else
    return false;
}

ResultIdList OutputFilter::SetUpstreamResults(){
  // check for global step values
  bool hasGlobalStepvalues = globalStepValueMap_.size() > 0;
  CF::StdVector<Double> globalStepValues;
  if (hasGlobalStepvalues) {
    globalStepValues.Reserve(globalStepValueMap_.size());
    std::map<UInt,Double>::iterator sIter = globalStepValueMap_.begin();
    for (;sIter != globalStepValueMap_.end();++sIter) {
      globalStepValues.Push_back(sIter->second);
    }
  }

  //determine requested results
  //loop over paramnode
  ResultIdList generated;
  ParamNodeList saveList = params_->Get("saveResults")->GetList("result");
  for(UInt aRes = 0; aRes < saveList.GetSize(); ++aRes){
    PtrParamNode rNode = saveList[aRes];
    std::string resultName = rNode->Get("resultName")->As<std::string>();
    upResNames.insert(resultName);
    uuids::uuid newId = RegisterUpstreamResult(resultName, uuids::nil_uuid());
    resultManager_->SetAsOutputResult(newId,true);
    if (hasGlobalStepvalues) {
      resultManager_->SetTimeLine(newId,globalStepValues);
    }
    generated.Push_back(newId);
  }
  return generated;

}

void OutputFilter::AdaptFilterResults(){
  //every result would must be valid by now. if not something is wrong
  ResultIdList::iterator iter = upResIds.Begin();
  std::set<uuids::uuid> toRemove;

  for(;iter!=upResIds.End();++iter){
    if(!resultManager_->GetExtInfo(*iter)->isValid){
      std::string s("The filter chain did not verify the result: ");
      s += resultManager_->GetExtInfo(*iter)->resultName;
      WARN(s);
      toRemove.insert(*iter);
    }
  }

  //remove invalid results
  std::set<uuids::uuid>::iterator aiter = toRemove.begin();
  for(;aiter != toRemove.end();++aiter){
    upResIds.Erase(upResIds.Find(*aiter));
  }

  //check that all upstream result pipelines are compatible
  //i.e. same time line and same result Type
  bool allSame = true;
  iter = upResIds.Begin();
  for(;iter!=upResIds.End();++iter){
    ResultIdList::iterator iter2 = iter;
    ResultManager::ConstInfoPtr res1 =resultManager_->GetExtInfo(*iter);
    const CF::StdVector<Double>& t1 = *res1->timeLine.get();
    for(;iter2!=upResIds.End();++iter2){
      ResultManager::ConstInfoPtr res2 =resultManager_->GetExtInfo(*iter2);
      const CF::StdVector<Double>& t2 = *res2->timeLine.get();
      allSame &= t1 == t2;
      allSame &= res1->dType == res2->dType;
    }
  }

  if(!allSame && globalStepValueMap_.size() == 0){
    WARN("Inconsistency in result infos from upstream filters detected. This indicates a problem in the filter pipeline definition.");
  }
  if(params_->GetParent()->Get("stepValueDefinition")->Has("dynamic")){
    std::cout << "Dynamic creation of output result timeline :" << std::endl;
    iter = upResIds.Begin();

    for(;iter!=upResIds.End();++iter){
      ResultManager::ConstInfoPtr res1 =resultManager_->GetExtInfo(*iter);
      if(res1->timeLine->GetSize() > 0){
        std::cout << "Using time line from result " << res1->resultName << "\n with "  \
                  << res1->timeLine->GetSize() << " steps from " << (*res1->timeLine.get())[0] << " to "\
                  << (*res1->timeLine.get())[res1->timeLine->GetSize()-1] << std::endl;
        if(res1->timeLine->GetSize()>1){
          std::cout << "Assumed delta " << (*res1->timeLine.get())[1]-(*res1->timeLine.get())[0] << std::endl;
        }
        CF::StdVector<Double> & timeVal = (*res1->timeLine.get());
        //CF::StdVector<UInt> & timeStep = (*res1->stepNumbers.get());
        for(UInt i=0;i<timeVal.GetSize(); ++i){
          globalStepValueMap_[i+1] = timeVal[i];
        }
        break;
      }
    }
    if(globalStepValueMap_.size() == 0)
      EXCEPTION("Could not determine timestep values for the computation. Check your pipeline inputs. Aborting...")

  }

  aStepIter_ = globalStepValueMap_.begin();
}

void OutputFilter::PrepareCalculation(){
  //now we register our results
  ResultIdList::iterator iter = upResIds.Begin();
  outFile_->Init(resultManager_->GetExtInfo(*iter)->ptGrid, false);

  for(;iter!=upResIds.End();++iter){
    StdVector< shared_ptr<BaseResult> >::iterator ResIter = resultManager_->GetBaseResultVector(*iter).Begin();
    ResultManager::ConstInfoPtr cInfo = resultManager_->GetExtInfo(*iter);
    PtrParamNode pNode = params_->Get("saveResults")->GetByVal("result","resultName",cInfo->resultName);
    UInt start = pNode->Get("saveBegin")->As<UInt>();
    UInt stop  = pNode->Get("saveEnd")->As<UInt>();
    UInt inc   = pNode->Get("saveInc")->As<UInt>();
    for(;ResIter != resultManager_->GetBaseResultVector(*iter).End();++ResIter){
      outFile_->RegisterResult(*ResIter,start,inc,stop,!cInfo->isMeshResult);
    }
  }
  iter = upResIds.Begin();
  //take first one and assume everybody is the same...
  BasePDE::AnalysisType aType;
  UInt numSteps = globalStepValueMap_.size();
  if(resultManager_->GetExtInfo(*iter)->dType == ExtendedResultInfo::COMPLEX){
    aType = BasePDE::HARMONIC;
  }else{
    aType = BasePDE::TRANSIENT;
  }

  //take the values from the first timeline.
  outFile_->BeginMultiSequenceStep(1,aType,numSteps);
}

}
