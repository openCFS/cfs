// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GenericInputFilter.cc
 *       \brief    <Description>
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================

#include <string>
#include <cstdlib>
#include <cstdio>

#include <Filters/Input/InputFilter.hh>
#include <Domain/Mesh/GridCFS/GridCFS.hh>

#include "DataInOut/DefineInOutFiles.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"



namespace CFSDat{

InputFilter::InputFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
            : BaseFilter(numWorkers,config,resMan){

  this->filtSteamType_ = INPUT_FILTER;

  CreateDummyCfsParamNode();
  PtrParamNode infoNode;

  std::string filename = params_->Get("inputFile")->GetChild()->Get("fileName")->As<std::string>();
  //Utilize the CFS way to create input files
  inFile_ = CoupledField::DefineInOutFiles::CreateSingleInputFileObject(filename,filterId_,params_->Get("inputFile")->GetChild(),infoNode);
  inFile_->InitModule();
  UInt maxDim = inFile_->GetDim();

  //now we can create the mesh
  //right now only full grids are supported
  std::string gridMode = params_->Get("gridType",ParamNode::EX)->As<std::string>();
  if(gridMode == "fullGrid"){
    ptGrid = new CF::GridCFS(maxDim,dummyXMLNode,infoNode,filterId_);
  }else if (gridMode == "simpleGrid"){
    ptGrid = new CF::GridCFS(maxDim,dummyXMLNode,infoNode,filterId_,false);
  }else{
    EXCEPTION("The chosen grid mode \'" << gridMode << "\' is not supported.")
  }

  inFile_->ReadMesh(ptGrid);
  //it would be nice not to finish the grid here
  //in order to let other filters add some entities
  //unfortunately this is not possible as we can not access anything
  //without it.
  //Another question: how can two inputs share a common grid?
  ptGrid->FinishInit();
  CreateAvailableResultInfos();
}



InputFilter::~InputFilter(){
  delete ptGrid;
}

bool InputFilter::Run(){
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();

  //=========================
  //TODO: this needs to be cleaner
  //=========================

  std::map<std::string,StdVector<Integer> > offsets;
  std::map<std::string,StdVector<boost::uuids::uuid> > offIds;
  for(; aIter != activeResults.end(); ++aIter){
    //loop over results and determine if they are offsets
    ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(*aIter);
    if(aInfo->resultName == "")
      continue;
    offsets[aInfo->resultName].Push_back(aInfo->timeStepOffset);
    offIds[aInfo->resultName].Push_back(*aIter);
  }

  //now we loop over the map and try to shift results and deactivate them
  std::map<std::string,StdVector<Integer> >::iterator mpIter = offsets.begin();
  for(; mpIter != offsets.end(); ++mpIter){
    if(mpIter->second.GetSize()==1){
      continue;
    }
    StdVector<Integer> cVec = mpIter->second;
    std::sort(cVec.Begin(),cVec.End());
    StdVector<boost::uuids::uuid>& uidVec = offIds[mpIter->first];

    //safe
    //for(UInt aV=0;aV<cVec.GetSize()-1;++aV){
    //  if(cVec[aV] != cVec[aV]-1){
    //    EXCEPTION("TIME BUFFER RESULT OFFSET TIMELINE IS NOT CONTINUOUS!");      
    //  }
    //}
    //ADDITIONAL: could be, that we compute only every other step.
    // need to check with timeline...

    for(UInt aV=0;aV<cVec.GetSize()-1;++aV){
      if(cVec[aV+1] != cVec[aV]+1)
        continue;
      CF::StdVector<UInt> eqnVec;
      //now shift the result and deactivate them
      UInt idx1 = mpIter->second.Find(cVec[aV]);
      UInt idx2 = mpIter->second.Find(cVec[aV+1]);
      //std::cout << "\tShifting result \'" << mpIter->first << "\' from step offset " <<  cVec[aV+1] << " to offset value " << cVec[aV] << std::endl;
      Vector<Double>& fullVec1 =  resultManager_->GetResultVector<Double>(uidVec[idx1],eqnVec);
      Vector<Double>& fullVec2 =  resultManager_->GetResultVector<Double>(uidVec[idx2],eqnVec);
      for(UInt aE=0;aE<fullVec1.GetSize();++aE){
        fullVec1[aE] = fullVec2[aE];
      } 
      resultManager_->DeactivateResult(uidVec[idx1]);
    }
  }
  //=========================
  // END TODO: this needs to be cleaner
  //=========================

  activeResults = resultManager_->GetActiveResults();
  aIter = activeResults.begin();
  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1)
      continue;
    ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(*aIter);
    ResultManager::ConstResPtr cRes = resultManager_->GetResultAdpter(*aIter);

    ExtendedResultInfo fileResult;
    for(UInt i=0;i<availInputResults_.GetSize();++i){
      if(availInputResults_[i].resultName == aInfo->resultName){
        fileResult = availInputResults_[i];
      }
    }
    if(fileResult.resultName != aInfo->resultName){
      EXCEPTION("Could not find requested result name in File");
    }

    if(resultManager_->GetExtInfo(*aIter)->dType == ExtendedResultInfo::COMPLEX){

    }else{

      CF::StdVector<UInt> eqnVec;
      Vector<Double>& fullVec =  resultManager_->GetResultVector<Double>(*aIter,eqnVec);

      Double reqValue = cRes->stepValue;
      CF::StdVector<Double>::iterator val= std::find_if(fileResult.timeLine->Begin(),fileResult.timeLine->End(), time_cmp(reqValue+*fileResult.timeLine->Begin(), 1E-6) );

      if(val == fileResult.timeLine->End()){
        if(reqValue+*fileResult.timeLine->Begin() < *(fileResult.timeLine->End()-1))
          std::cerr  << "ERROR: can not find a timestep for time value \'" << reqValue << "\' Either there are no more timesteps or floating point conversion errors occured" << std::endl;
        else
          std::cout << "\t\t\t WARN: Tying to access time data beyond source timeline. Results will be wrong for these steps. Take care if trying to append!" << std::endl;
        fullVec.Init();
        continue;
      }
      CF::SolutionType solType = fileResult.resultType;
      std::set<std::string>::const_iterator regIter = aInfo->regNames->begin();
      UInt idx = std::distance(fileResult.timeLine->Begin(), val);
      UInt stepNumber = (*fileResult.stepNumbers.get())[idx];
      for(; regIter != aInfo->regNames->end(); ++regIter){
        Vector<Double> resVec;
        shared_ptr<BaseResult> inResult = inFile_->GetResult(fileResult.sequenceStep,stepNumber,solType,*regIter);
        try{
          Result<Double>* myResult = dynamic_cast<Result<Double>* >(inResult.get());
          resVec =   myResult->GetVector();
        }catch(...){
          EXCEPTION("Cannot cast to desired vector type. Are you trying to load real data into a harmonic computation?");
        }
        RegionIdType regId = aInfo->ptGrid->GetRegion().Parse(*regIter);
        eqnVec.Clear(true);
        cRes->mapping->GetRegionEquations(eqnVec,regId);
        for(UInt aEq = 0; aEq<eqnVec.GetSize();++aEq){
//TODO REGION_PROBLEM: if more than one target-region is specified in the xml-scheme, e.g. fluid and one boundary, then
          //the result vector contains exactly ONE more entry than without the boundary, independent of the number of
          //additional boundary regions
          fullVec[eqnVec[aEq]] = resVec[aEq];
        }
      }
    }
  }

  return true;
}

ResultIdList InputFilter::SetUpstreamResults(){
  //nothing to be done here


  ResultIdList newList;
  return newList;
}

void InputFilter::AdaptFilterResults(){
  //here we need to set every information about our filter
  for(UInt aRes = 0; aRes < filterResIds.GetSize(); aRes++){
    ResultManager::ConstInfoPtr curResInfo = resultManager_->GetExtInfo(filterResIds[aRes]);
    ExtendedResultInfo fileResult;
    for(UInt i=0;i<availInputResults_.GetSize();++i){
      if(availInputResults_[i].resultName == curResInfo->resultName){
        fileResult = availInputResults_[i];
      }
    }
    if(fileResult.resultName != curResInfo->resultName){
      EXCEPTION("Could not find requested result name in File");
    }
    //fill the information
    //do not override the timeline if it is already filled...
    if(curResInfo->timeLine->GetSize() != 0){
      //apparently somebody downstream has already something defined
      //lets just check if we are compatible
      CF::StdVector<Double>::iterator timeIter = curResInfo->timeLine->Begin();
      Double startTime = *fileResult.timeLine->Begin();
      bool allOK = true;
      for(;timeIter != curResInfo->timeLine->End();++timeIter){
        allOK &= std::find_if(fileResult.timeLine->Begin(),fileResult.timeLine->End(), time_cmp(startTime+*timeIter, 1E-8) ) != fileResult.timeLine->End();
      }
      if(!allOK)
        EXCEPTION("The input filter cannot provide every timestep which is requested. Check the definition of input results")
    }else{
      //remove offset per default
      StdVector<Double>& curT = (*fileResult.timeLine.get());
      StdVector<Double> tfValuesOffset(curT.GetSize());
      for(UInt aTime =0;aTime<(*fileResult.timeLine.get()).GetSize();aTime++){
        tfValuesOffset[aTime] = curT[aTime] - curT[0];
      }
      resultManager_->SetTimeLine(filterResIds[aRes],tfValuesOffset);
    }
    resultManager_->SetDType(filterResIds[aRes],fileResult.dType);
    resultManager_->SetRegionNames(filterResIds[aRes],(*fileResult.regNames.get()));
    resultManager_->SetCFormat(filterResIds[aRes],fileResult.complexFormat);
    resultManager_->SetEntryType(filterResIds[aRes],fileResult.entryType);
    resultManager_->SetUnit(filterResIds[aRes],fileResult.unit);
    resultManager_->SetDofNames(filterResIds[aRes],fileResult.dofNames);
    resultManager_->SetDefOn(filterResIds[aRes],fileResult.definedOn);
    resultManager_->SetGrid(filterResIds[aRes],ptGrid);
    resultManager_->SetMeshResult(filterResIds[aRes],true);
    resultManager_->SetValid(filterResIds[aRes]);
  }
}


void InputFilter::CreateAvailableResultInfos(){

  //we traverse the input and generate our results on the fly

  std::map<UInt, BasePDE::AnalysisType> analysis;
  std::map<UInt, UInt> numSteps;
  inFile_->GetNumMultiSequenceSteps(analysis,numSteps);

  //now we obtain a List of all Results for each sequence Step
  std::map<UInt, BasePDE::AnalysisType>::iterator anaIter = analysis.begin();
  for(;anaIter!=analysis.end();++anaIter){
    std::cout << "\t\t-> Input Filter with id \"" << this->filterId_ << "\" detected sequence step #" << anaIter->first << " with " << numSteps[anaIter->first] << " steps" << std::endl;
    CF::StdVector<boost::shared_ptr<ResultInfo> > infos;
    inFile_->GetResultTypes(anaIter->first,infos,false);
    std::map<UInt, Double> steps;
    CF::StdVector<boost::shared_ptr<EntityList> > entList;
    for(UInt aRes=0;aRes < infos.GetSize(); ++aRes){
      str1::shared_ptr<ResultInfo> curRes = infos[aRes];
      inFile_->GetStepValues(anaIter->first,curRes,steps,false);
      inFile_->GetResultEntities(anaIter->first,curRes,entList,false);

      //fill Info
      ExtendedResultInfo newRes;
      newRes.ImportResultInfo(curRes);
      newRes.sequenceStep = anaIter->first;

      if(anaIter->second == BasePDE::HARMONIC)
        newRes.dType = ExtendedResultInfo::COMPLEX;
      else
        newRes.dType = ExtendedResultInfo::DOUBLE;

      //import time values
      inFile_->GetStepValues(anaIter->first,curRes,steps,false);
      std::map<UInt,Double>::iterator iter = steps.begin();
      newRes.timeLine->Reserve(steps.size());
      newRes.stepNumbers->Reserve(steps.size());
      for(;iter!=steps.end();++iter){
        newRes.timeLine->Push_back(iter->second);
        newRes.stepNumbers->Push_back(iter->first);
      }
      for(UInt aList=0;aList<entList.GetSize();aList++){
        std::string rName = entList[aList]->GetName();
        newRes.regNames->insert(rName);
        //we do not fill the entity numbers right now
        //as this can be huge. lets do that on demand
      }
      availInputResults_.Push_back(newRes);
      filtResNames.insert(newRes.resultName);
    }
  }

}

void InputFilter::CreateDummyCfsParamNode(){

  PtrParamNode meshInputNode = params_->Get("inputFile")->GetChild();
  dummyXMLNode.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));
  CoupledField::PtrParamNode iNode = dummyXMLNode->Get("fileFormats",ParamNode::APPEND)->Get("input",ParamNode::INSERT);
  iNode->AddChildNode(meshInputNode);

  //create domain node
  //UInt dim = inFile_->GetDim();
  CoupledField::PtrParamNode dNode = dummyXMLNode->Get("domain",ParamNode::APPEND);
  //dNode->Get("geometryType",ParamNode::APPEND)->SetValue(itoa(dim));

}

}
