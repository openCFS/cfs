// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ResultManager.cc
 *       \brief    <Description>
 *
 *       \date     Nov 26, 2015
 *       \author   ahueppe
 */
//================================================================================================

#include "ResultManager.hh"
#include "MatVec/Vector.hh"
#include <set>
#include <algorithm>

namespace CFSDat{

ResultManager::ResultManager(){
  isFinalized_ = false;
}

ResultManager::~ResultManager(){

}

uuids::uuid ResultManager::AddResult(std::string name, uuids::uuid filterTag, CF::StdVector<Integer> timeCacheSteps){
  uuids::uuid generated = boost::uuids::random_generator()();
  //if we do not have an offset vector we proceed with the standard
  //operation, if we have one given, we add a master result...
  if(timeCacheSteps.GetSize() == 0){
    //we add a new value for our map
    std::pair<InfoPtr,ResPtr> newInfo;
    newInfo.first.reset(new ExtendedResultInfo);
    //leave the second part as null for the moment, this will be done in finalize
    //just set the name
    newInfo.first->resultName = name;
    newInfo.first->generatorId = filterTag;
    resultMap_[generated] = newInfo;
  }else{
    //master result case...
    UInt numOffsets  = timeCacheSteps.GetSize();
    for(UInt aStep =0;aStep<numOffsets;++aStep){
      //generate new slave result
      uuids::uuid cSlave = boost::uuids::random_generator()();
      std::pair<InfoPtr,ResPtr> newInfo;
      newInfo.first.reset(new ExtendedResultInfo);
      newInfo.first->resultName = name;
      newInfo.first->generatorId = filterTag;
      newInfo.first->timeStepOffset = timeCacheSteps[aStep];
      newInfo.first->timeCacheMasterId = generated;
      resultMap_[cSlave] = newInfo;
      masterResMap_[generated][timeCacheSteps[aStep]] = cSlave;
    }
    masterInfos_[generated].reset(new ExtendedResultInfo);
    masterInfos_[generated]->isTimeCacheMaster = true;
    masterInfos_[generated]->resultName = name;
    masterInfos_[generated]->generatorId = filterTag;
    //now we already ensure that all shared pointers are identical
    std::map<Integer, uuids::uuid>::iterator resIter = masterResMap_[generated].begin();
    for(;resIter != masterResMap_[generated].end();++resIter){
      uuids::uuid cNum = resIter->second;
      resultMap_[cNum].first->timeLine = masterInfos_[generated]->timeLine;
      resultMap_[cNum].first->stepNumbers = masterInfos_[generated]->stepNumbers;
      resultMap_[cNum].first->regNames = masterInfos_[generated]->regNames;
      resultMap_[cNum].first->entityNumbers = masterInfos_[generated]->entityNumbers;
      resultMap_[cNum].first->eqnNumbers = masterInfos_[generated]->eqnNumbers;
    }
  }
  return generated;
}

bool ResultManager::IsResultVecUpToDate(uuids::uuid requestedId, Integer offsetStep){
  checkFinalized();
  uuids::uuid remappedId = remapMasterResults(requestedId,offsetStep);
  return resultMap_[remappedId].second->isUpToDate;
}
template<typename T>
CF::Vector<T>&  ResultManager::GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers, Integer timeStepOffset){
  //performance tests have reveals that all syntactic sugar
  //e.g. by providing a view on the vector or anything have not advantage
  //in runtime over just returning an index array which can be used by the caller
  //in addition, the memory requirements for an additional integer array are quite low
  //quite annoying nonetheless, if the requested result is defined on all entites,
  // a vector from 0...N is returned... bullshit
  checkFinalized();
  uuids::uuid remappedId = remapMasterResults(requestedId,timeStepOffset);
  ResultAdaptor<T> * resAdapt = dynamic_cast<ResultAdaptor<T> * >(resultMap_[remappedId].second.get());
  eqnNumbers = (*resultMap_[remappedId].first->eqnNumbers.get());
  return resAdapt->resultVector;
}

//CF::Vector<Complex>& ResultManager::GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers, Integer timeStepOffset){
//  checkFinalized();
//  uuids::uuid remappedId = remapMasterResults(requestedId,timeStepOffset);
//
//  ResultAdaptor<Complex> * resAdapt = dynamic_cast<ResultAdaptor<Complex> * >(resultMap_[remappedId].second.get());
//  //retVec = dynamic_cast<CF::Vector<Complex> & >(resultMap_[remappedId].second->GetSingleVector());
//  eqnNumbers = (*resultMap_[remappedId].first->eqnNumbers.get());
//  return resAdapt->resultVector;
//}

CF::StdVector<shared_ptr<BaseResult> >& ResultManager::GetBaseResultVector(uuids::uuid requestedId, Integer timeStepOffset){
  checkFinalized();
  uuids::uuid remappedId = remapMasterResults(requestedId,timeStepOffset);
  return resultMap_[remappedId].second->baseResultVector;
}

ResultManager::ConstResPtr ResultManager::GetResultAdpter(uuids::uuid requestedId, Integer timeStepOffset){
  checkFinalized();
  uuids::uuid remappedId = remapMasterResults(requestedId,timeStepOffset);
  return resultMap_[remappedId].second;
}

Double ResultManager::GetStepValue(uuids::uuid requestedId, Integer offsetStep){
  checkFinalized();
  uuids::uuid remappedId = remapMasterResults(requestedId,offsetStep);
  return resultMap_[remappedId].second->stepValue;
}

void ResultManager::SetTimeValue(uuids::uuid requestedId, Double value){
  checkFinalized();

  if(checkResultDef(requestedId)){
    std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[requestedId].begin();
    for(;assocIter != masterResMap_[requestedId].end();++assocIter){
      if(resultMap_[assocIter->second].first->timeCacheMasterId == requestedId){
        StdVector<Double>& tsteps = (*resultMap_[assocIter->second].first->timeLine.get());
        //assume equidistant otherwise we would need to search now
        Double delta = tsteps[1] - tsteps[0];
        resultMap_[assocIter->second].second->stepValue = value+(resultMap_[assocIter->second].first->timeStepOffset*delta);
        resultMap_[assocIter->second].second->GetSingleVector()->Init();
        resultMap_[assocIter->second].second->isUpToDate = false;
      }
    }
  }else{
    resultMap_[requestedId].second->stepValue = value;
    resultMap_[requestedId].second->isUpToDate = false;
  }
}

//  //this is avery annoying thing. since different input files have
//  //different accuracies of storing the time, we need to check for the
//  //correct step in the timeline
//
//  UInt tIdx = 0;
//
//  //check for master result
//  if(masterResMap_.find(requestedId) != masterResMap_.end()){
//    //master has no results but his children. so we take the offset
//    //here we could also establish, that we copy some values to a avoid reread
//    std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[requestedId].begin();
//    for(;assocIter != masterResMap_[requestedId].end();++assocIter){
//      StdVector<Double>& tsteps = (*resultMap_[assocIter->second].first->timeLine.get());
//      StdVector<Double>::iterator position;
//      position = std::find_if(tsteps.Begin(),tsteps.End(), time_cmp(value, 1E-8) );
//      tIdx = std::distance(tsteps.Begin(), position);
//      if(tIdx == 0 && assocIter->first < 0 ){
//        //ok here we just set the value to a negative value and deactivate the result
//        resultMap_[assocIter->second].second->stepValue = -1.0;
//        resultMap_[assocIter->second].second->GetSingleVector()->Init();
//        resultMap_[assocIter->second].second->isUpToDate = true;
//        activeIds_.erase(assocIter->second);
//      }else{
//        //in case of time caching results, it happens that the value in the future
//        //is beyond the maximum values of the the timeline. this may not be wrong e.g. if
//        // more results are available upstream then the user wanted. in this case, we just
//        //extrapolate values. The corresponding input filter will hopefully deal with the case of invalid data...
//        UInt newIdx = tIdx + assocIter->first;
//        UInt newStep = 0;
//        UInt newVal = 0;
//        if(newIdx >= tsteps.GetSize()){
//          //choose the last and the second last element and extrapolate linearly
//          UInt last = tsteps.GetSize()-1;
//          Double delta = tsteps[last] - tsteps[last-1];
//          UInt diff = abs(newIdx - last);
//          resultMap_[assocIter->second].second->stepValue = tsteps[last] + diff*delta;
//        }else{
//          resultMap_[assocIter->second].second->stepValue = tsteps[newIdx];
//        }
//        resultMap_[assocIter->second].second->isUpToDate = false;
//      }
//    }
//  }else{
//    StdVector<Double>& tsteps = (*resultMap_[requestedId].first->timeLine.get());
//    StdVector<Double>::iterator position;
//    position = std::find_if(tsteps.Begin(),tsteps.End(), time_cmp(value, 1E-8) );
//    tIdx = std::distance(tsteps.Begin(), position);
//    resultMap_[requestedId].second->stepValue = *position;
//    resultMap_[requestedId].second->stepNumber = (*resultMap_[requestedId].first->stepNumbers.get())[tIdx];
//    resultMap_[requestedId].second->isUpToDate = false;
//  }
//}

ResultManager::ConstInfoPtr ResultManager::GetExtInfo(uuids::uuid resId){
  //this check is tedious can we think of anything else?
  if(resultMap_.find(resId) == resultMap_.end()){
    //try to look within the master results
    if(masterInfos_.find(resId) != masterInfos_.end()){
      return masterInfos_[resId];
    }
    EXCEPTION("Requested Id which is not generated! Make sure every call "
                  "to ResultManager provides valid resultIds.");
  }

  return resultMap_[resId].first;
}

//has anybody an idea how to avoid this piece of code?
//and do the same without templates?
#define RESULT_MANAGER_COPY_FIELD(resId, fName , toCopy)                                     \
    bool isMaster = checkResultDef(resId);                                                   \
    if(isMaster){                                                                            \
      std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[resId].begin();    \
      for(;assocIter != masterResMap_[resId].end();++assocIter){                             \
        resultMap_[assocIter->second].first->fName = toCopy;                                 \
      }                                                                                      \
      masterInfos_[resId]->fName = toCopy;                                                   \
    }else{                                                                                   \
      InfoPtr curInfo = resultMap_[resId].first;                                             \
      curInfo->fName = toCopy;                                                               \
      if(curInfo->timeCacheMasterId!=uuids::nil_uuid() )                                     \
        masterInfos_[curInfo->timeCacheMasterId]->fName = toCopy;                             \
    }                                                                                        \

#define RESULT_MANAGER_COPY_SH_PTR_FIELD(resId, fName , toCopy)                              \
    bool isMaster = checkResultDef(resId);                                                   \
    if(isMaster){                                                                            \
      std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[resId].begin();    \
      for(;assocIter != masterResMap_[resId].end();++assocIter){                             \
        (*resultMap_[assocIter->second].first->fName.get()) = toCopy;                        \
      }                                                                                      \
      (*masterInfos_[resId]->fName.get()) = toCopy;                                          \
    }else{                                                                                   \
      (*resultMap_[resId].first->fName.get()) = toCopy;                                      \
    }                                                                                        \

#define RESULT_MANAGER_COPY_SH_PTR_SET(resId, fName , toCopy)                                \
    bool isMaster = checkResultDef(resId);                                                   \
    if(isMaster){                                                                            \
      std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[resId].begin();    \
      for(;assocIter != masterResMap_[resId].end();++assocIter){                             \
        resultMap_[assocIter->second].first->fName->clear();                                 \
        resultMap_[assocIter->second].first->fName->insert(toCopy.begin(),toCopy.end());     \
      }                                                                                      \
      masterInfos_[resId]->fName->clear();                                                   \
      masterInfos_[resId]->fName->insert(toCopy.begin(),toCopy.end());                       \
    }else{                                                                                   \
      resultMap_[resId].first->fName->clear();                                               \
      resultMap_[resId].first->fName->insert(toCopy.begin(),toCopy.end());                   \
    }                                                                                        \

#define RESULT_MANAGER_OBTAIN_FIELD(resId, fName)                  \
    bool isMaster = checkResultDef(resId);                         \
    if(isMaster){                                                  \
      return masterInfos_[resId]->fName;                           \
    }else{                                                         \
      return resultMap_[resId].first->fName;                       \
    }                                                              \


#define RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId, fName)           \
    bool isMaster = checkResultDef(resId);                         \
    if(isMaster){                                                  \
      return (*masterInfos_[resId]->fName.get());                  \
    }else{                                                         \
      return (*resultMap_[resId].first->fName.get());              \
    }                                                              \

void ResultManager::CopyResultData(uuids::uuid srcId, uuids::uuid trgId){
  checkNotFinalized();
  bool isMasterTrg = checkResultDef(trgId);
  bool isMasterSrc = checkResultDef(srcId);
  InfoPtr srcInfo = (isMasterSrc)? masterInfos_[srcId] : resultMap_[srcId].first;

  if(isMasterTrg){
    std::map<Integer, uuids::uuid>::iterator assocIter =  masterResMap_[trgId].begin();
    for(;assocIter !=  masterResMap_[trgId].end();++assocIter){
      resultMap_[assocIter->second].first->dType = srcInfo->dType;
      resultMap_[assocIter->second].first->regNames = srcInfo->regNames;
      resultMap_[assocIter->second].first->complexFormat = srcInfo->complexFormat;
      resultMap_[assocIter->second].first->definedOn = srcInfo->definedOn;
      resultMap_[assocIter->second].first->dofNames = srcInfo->dofNames;
      resultMap_[assocIter->second].first->entryType = srcInfo->entryType;
      resultMap_[assocIter->second].first->ptGrid = srcInfo->ptGrid;
    }
    (*masterInfos_[trgId]->entityNumbers.get()) = (*srcInfo->entityNumbers.get());
    (*masterInfos_[trgId]->eqnNumbers.get()) = (*srcInfo->eqnNumbers.get());
    (*masterInfos_[trgId]->timeLine.get()) = (*srcInfo->timeLine.get());
    (*masterInfos_[trgId]->stepNumbers.get()) = (*srcInfo->stepNumbers.get());
  }else{
    resultMap_[trgId].first->dType = srcInfo->dType;
    resultMap_[trgId].first->regNames = srcInfo->regNames;
    resultMap_[trgId].first->complexFormat = srcInfo->complexFormat;
    resultMap_[trgId].first->definedOn = srcInfo->definedOn;
    resultMap_[trgId].first->dofNames = srcInfo->dofNames;
    resultMap_[trgId].first->entryType = srcInfo->entryType;
    resultMap_[trgId].first->ptGrid = srcInfo->ptGrid;

    (*resultMap_[trgId].first->entityNumbers.get()) = (*srcInfo->entityNumbers.get());
    (*resultMap_[trgId].first->eqnNumbers.get()) = (*srcInfo->eqnNumbers.get());
    (*resultMap_[trgId].first->timeLine.get()) = (*srcInfo->timeLine.get());
    (*resultMap_[trgId].first->stepNumbers.get()) = (*srcInfo->stepNumbers.get());
  }
}

CF::StdVector<Double> ResultManager::GetTimeLine(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,timeLine)
}

void ResultManager::SetTimeLine(uuids::uuid resId, CF::StdVector<Double> tVec){
  checkNotFinalized();
  RESULT_MANAGER_COPY_SH_PTR_FIELD(resId,timeLine,tVec)
}

CF::StdVector<UInt> ResultManager::GetStepNumbers(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,stepNumbers)
}

void ResultManager::SetStepNumbers(uuids::uuid resId, CF::StdVector<UInt> tVec){
  checkNotFinalized();
  RESULT_MANAGER_COPY_SH_PTR_FIELD(resId,stepNumbers,tVec)
}

ExtendedResultInfo::ResDType ResultManager::GetDType(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,dType)
}

void ResultManager::SetDType(uuids::uuid resId, ExtendedResultInfo::ResDType dType){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,dType,dType)
}

ExtendedResultInfo::EntityUnknownType ResultManager::GetDefOn(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,definedOn)
}

void ResultManager::SetDefOn(uuids::uuid resId, ExtendedResultInfo::EntityUnknownType eType){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,definedOn,eType)
}


std::set<std::string>  ResultManager::GetRegionNames(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,regNames)
}

void ResultManager::SetRegionNames(uuids::uuid resId, std::set<std::string>  regions){
  checkNotFinalized();
  RESULT_MANAGER_COPY_SH_PTR_SET(resId,regNames,regions)
}

void ResultManager::SetMeshResult(uuids::uuid resId, bool isMesh){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,isMeshResult,isMesh)
}

CF::StdVector<UInt> ResultManager::GetEntityNumbers(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,entityNumbers)
}

void ResultManager::SetEntityNumbers(uuids::uuid resId, CF::StdVector<UInt> numbers){
  checkNotFinalized();
  RESULT_MANAGER_COPY_SH_PTR_FIELD(resId,entityNumbers,numbers)
}

StdVector<UInt> ResultManager::GetEqnNumbers(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,eqnNumbers)
}

CF::ComplexFormat ResultManager::GetCFormat(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,complexFormat)
}

void ResultManager::SetCFormat(uuids::uuid resId, CF::ComplexFormat cType){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,complexFormat,cType)
}

ExtendedResultInfo::EntryType ResultManager::GetEntryType(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,entryType)
}

void ResultManager::SetEntryType(uuids::uuid resId, ExtendedResultInfo::EntryType rType){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,entryType,rType)
}

std::string ResultManager::GetUnit(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,unit)
}

void ResultManager::SetUnit(uuids::uuid resId, std::string newUnit){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,unit,newUnit)
}

std::string ResultManager::GetResultName(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,resultName)
}

CF::StdVector<std::string> ResultManager::GetDofNames(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_FIELD(resId,dofNames)
}

void ResultManager::SetDofNames(uuids::uuid resId, CF::StdVector<std::string> newDofNames){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,dofNames,newDofNames)
}

void ResultManager::SetGrid(uuids::uuid resId, Grid* ptGrid){
  checkNotFinalized();
  RESULT_MANAGER_COPY_FIELD(resId,ptGrid,ptGrid)
}

void ResultManager::SetValid(uuids::uuid resId){
  //we set the result as valid anyway
  //in case of master results we check if every associated result is valid
  //and check the master only if this is the case
  bool isMaster = checkResultDef(resId);
  if(isMaster){
    WARN("trying to check a master result as valid. this is not expected and indicates a bug. Please send a report.");
  }else{
    resultMap_[resId].first->isValid = true;
    uuids::uuid masterId = resultMap_[resId].first->timeCacheMasterId;
    //now we loop over associated results and check for validity
    if(masterId != uuids::nil_uuid()){
      std::map<Integer, uuids::uuid>::iterator uIter = masterResMap_[masterId].begin();
      bool allValid = true;
      for(;uIter!=masterResMap_[masterId].end();++uIter){
        allValid &= resultMap_[uIter->second].first->isValid;
      }
      masterInfos_[masterId]->isValid = allValid;
    }
  }
}

void ResultManager::Finalize(){
  //ok, how do we do that:
  // 1. this will be a costly algorithm in which we try to identify the following
  //     a. Which results are completely the same and combine those
  //     b. From the remaining differences we try to identify if they share
  //         I. Same EntryType and same EntityType and same GridId -> oneEqnMap
  //        II. Condition I is true and same EntityList -> commonEqnVector
  //       III. Same timeline -> Delete all timelines except one and modify pointers
  //        IV. ALL the same except EntityList -> combineEqnVec and create common resultVector
  // 2. Which datastrcutures are necessary?
  //    1. We copy all resultIds into a set
  //    2. We iterate over all results if result is identical to current result, its removed from
  //       temporary set and pushed into a corresponding vector
  //       next result in remaining set and so on -> vector of vectors
  //    3. For the remainder we will check for 1b as follows
  //       loop over head entries in each vector and compare COnd 1 with all others
  //       etc. all indices which fulfill these conditions will pushed into corresponding
  //       vectors. after that, we will create eqnMaps as required
  //    4. Loop over all results marked with commonEqVector, obtain eqnVec from Map and
  //       adjust pointers
  //    5. Loop over all results marked with commonTimeLine and set pointers
  //    6. Loop over all results marked with common ResultVec, combine the entitylists and create
  //       a big result vector according to them, set unique equationvectors
  //    7. Loop over all remaining results and create unique structures
  //    8. It turns out, that shared pointers are not very well suited here due to the possibility of
  //       cyclic references. e.g. result1 gets destroyed but holds a shared_ptr to result2

  // last question: what about writing of results
  //       pos1. it gets converted if a result is requested as an output result
  //       pos2. Output results are marked on creation and cannot be combined in a simple manner
  //       pos3. Base results are cached w.r.t. the entitylists etc.

  std::set<InfoPtr> allInfos;
  std::set<Grid*> allGrids;
  UuidMap::iterator uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    if(uuidIter->second.first->isValid){
      allInfos.insert(uuidIter->second.first);
      allGrids.insert(uuidIter->second.first->ptGrid);
//      std::cerr << uuidIter->second.first->resultName << " ";
//      std::cerr << uuidIter->second.first->ptGrid << std::endl;
    }
  }

  //check if we have results which are consistent with caching results and add those to master results
  MasterUuidMap::iterator maIter = masterResMap_.begin();
  for(;maIter!=masterResMap_.end();++maIter){
    uuidIter = resultMap_.begin();
    for(;uuidIter != resultMap_.end();++uuidIter){
      if(uuidIter->second.first->timeStepOffset == 0){
        if(checkEqualStageStrict(uuidIter->second.first,masterInfos_[maIter->first]) ){
          if(maIter->second.find(0) != maIter->second.end()){
            uuidIter->second.first = masterInfos_[maIter->second[0]];
          }else{
            maIter->second[0] = uuidIter->first;
            //modify shared pointer arrays
            uuidIter->second.first->timeLine = masterInfos_[maIter->first]->timeLine;
            uuidIter->second.first->stepNumbers = masterInfos_[maIter->first]->stepNumbers;
            uuidIter->second.first->regNames = masterInfos_[maIter->first]->regNames;
            uuidIter->second.first->entityNumbers = masterInfos_[maIter->first]->entityNumbers;
            uuidIter->second.first->eqnNumbers = masterInfos_[maIter->first]->eqnNumbers;
          }
        }
      }
    }
  }

  //check consistency of time cache master results
  maIter = masterResMap_.begin();
  for(;maIter!=masterResMap_.end();++maIter){
    std::map<Integer, uuids::uuid>::iterator sIter = maIter->second.begin();
    if(!checkEqualStageStrict(masterInfos_[maIter->first],resultMap_[sIter->second].first)){
      EXCEPTION("Detected inconsistency in time cache result " << masterInfos_[maIter->first]->resultName << ". This indicates a bug.")
    }
  }

  //now it it time to create results for each unique result (compression not available yet)
  uuidIter = resultMap_.begin();
  std::set<uuids::uuid> toRemove;
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    if(!cInfo->isValid){
      WARN("Result " <<  cInfo->resultName << " is marked invalid after initialization phase. Going to remove it from pool!");
      toRemove.insert(uuidIter->first);
      continue;
    }

    if(cInfo->dType == ExtendedResultInfo::COMPLEX){
      uuidIter->second.second = str1::shared_ptr<GenericResultAdapter>(new ResultAdaptor<CF::Complex>());
    }else if(cInfo->dType == ExtendedResultInfo::DOUBLE){
      uuidIter->second.second = str1::shared_ptr<GenericResultAdapter>(new ResultAdaptor<CF::Double>());
    }else{
      EXCEPTION("Only Complex and Double results supported yet")
    }
  }
  //remove invalid results
  std::set<uuids::uuid>::iterator aiter = toRemove.begin();
  for(;aiter != toRemove.end();++aiter){
    resultMap_.erase(*aiter);
    activeIds_.erase(*aiter);
  }
  //get unique grid pointers
//
//  std::cerr << allInfos.size() << " unique resultInfos are to cover. ";
//  std::cerr << allGrids.size() << " unique grids are to cover. ";



  //start with the first result and obtain a vector of equations
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    ResPtr cRes = uuidIter->second.second;
    if(!cRes->isUpToDate){
      CreateEqnMapping(cInfo,cRes);
      CreateResultVector(cInfo,cRes);
      cRes->isUpToDate = true;
      UuidMap::iterator otherIter = resultMap_.begin();
      //now we loop over all results and try to apply same pointers as possible
      for(;otherIter != resultMap_.end();++otherIter){
        InfoPtr coInfo = otherIter->second.first;
        ResPtr coRes = otherIter->second.second;
        if(!coRes->isUpToDate){
          //try to copy shared pointer
          if(checkEqualEqnMap(uuidIter->second, otherIter->second)){
            coRes->mapping = cRes->mapping;
            if(checkEqualEqnVec(uuidIter->second, otherIter->second)){
              coInfo->eqnNumbers = cInfo->eqnNumbers;
            }
            CreateResultVector(coInfo,coRes);
            coRes->isUpToDate = true;
          }
        }
      }
    }
  }

  //create output results for the remainder
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    ResPtr cRes = uuidIter->second.second;
    if(!cInfo->isOutput)
      continue;

    //no offset results here...
    StdVector< str1::shared_ptr<CF::BaseResult> >& resVec = cRes->baseResultVector;

    if(cInfo->isMeshResult){
      //this loop body is adapted for mesh results. no matter what we have in entitynumbers from
      // upstream, we write full regions
      std::set<std::string>::iterator regIter = cInfo->regNames->begin();
      for(;regIter != cInfo->regNames->end();++regIter){
        UInt idx = resVec.GetSize();

        if(cInfo->dType == ExtendedResultInfo::COMPLEX){
          resVec.Push_back(str1::shared_ptr< BaseResult >(new Result<CF::Complex>() ) );
        }else{
          resVec.Push_back(str1::shared_ptr< BaseResult >(new Result<CF::Double>() ) );
        }

        str1::shared_ptr<EntityList> aList;
        if(cInfo->definedOn == ExtendedResultInfo::NODE){
          aList = cInfo->ptGrid->GetEntityList(EntityList::NODE_LIST,*regIter);
        }else if(cInfo->definedOn == ExtendedResultInfo::ELEMENT){
          aList = cInfo->ptGrid->GetEntityList(EntityList::ELEM_LIST,*regIter);
        }
        resVec[idx]->SetEntityList(aList);
        resVec[idx]->GetSingleVector()->Resize(aList->GetSize());
        resVec[idx]->SetResultInfo(cInfo->GetResultInfo());
      }
    }else{
      UInt idx = resVec.GetSize();
      if(cInfo->dType == ExtendedResultInfo::COMPLEX){
        resVec.Push_back(str1::shared_ptr< BaseResult >(new Result<CF::Complex>() ) );
      }else{
        resVec.Push_back(str1::shared_ptr< BaseResult >(new Result<CF::Double>() ) );
      }
      resVec[idx]->SetResultInfo(cInfo->GetResultInfo());
    }
  }

  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    uuidIter->second.second->isUpToDate = false;
  }


//  uuidIter = resultMap_.begin();
//  for(;uuidIter != resultMap_.end();++uuidIter){
//    print_ExtInfoFields((*uuidIter->second.first.get()));
//    std::cout <<uuidIter->second.second->mapping << std::endl;
//    std::cout << "--------------------------------------------------" << std::endl;
//  }
  this->isFinalized_ = true;

}
void ResultManager::CreateEqnMapping(InfoPtr cInfo, ResPtr cRes){
  if(cInfo->entryType == CF::ResultInfo::SCALAR){
    cRes->mapping.reset(new EqnMapSimple(cInfo->definedOn, cInfo->ptGrid));
  }else{
    cRes->mapping.reset(new EqnMapSimple(cInfo->definedOn, cInfo->ptGrid, cInfo->dofNames.GetSize()));
  }
  std::set<std::string>::iterator regIter = cInfo->regNames->begin();
  for(;regIter!=cInfo->regNames->end();++regIter){
    cRes->mapping->AddRegion(cInfo->ptGrid->GetRegion().Parse(*regIter));
  }

  cRes->mapping->Finalize();

  if(cInfo->entityNumbers->GetSize() != 0){
    cRes->mapping->GetSubsetEquations((*cInfo->eqnNumbers.get()),(*cInfo->entityNumbers.get()));
  }
}

void ResultManager::CreateResultVector(InfoPtr cInfo, ResPtr cRes){
  //first we create a Result
  //obtain initial equation vector
  if(cInfo->entityNumbers->GetSize() == 0){
    //we just resize the solution vector to the global number of equations
    cRes->GetSingleVector()->Resize(cRes->mapping->GetNumEquations());
  }else{
    //apparently we only need a subset of equations, so we try to cope with that
    cRes->GetSingleVector()->Resize(cInfo->eqnNumbers->GetSize());
  }
  cRes->GetSingleVector()->Init();

}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template CF::Vector<Double>&  ResultManager::GetResultVector(uuids::uuid, CF::StdVector<UInt>&, Integer);
  template CF::Vector<Complex>&  ResultManager::GetResultVector(uuids::uuid, CF::StdVector<UInt>&, Integer);
#endif
}

