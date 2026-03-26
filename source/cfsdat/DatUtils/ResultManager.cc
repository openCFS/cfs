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
// #include <algorithm>
#include <cmath>

namespace CFSDat{

ResultManager::ResultManager(){
  isFinalized_ = false;
}

ResultManager::~ResultManager(){

}

uuids::uuid ResultManager::AddResult(std::string name, uuids::uuid filterTag){
  return AddResult(name, filterTag, 0, 0);
}

uuids::uuid ResultManager::AddResult(std::string name, uuids::uuid filterTag,
       Integer minStepOffset, Integer maxStepOffset){
  uuids::uuid generated = boost::uuids::random_generator()();
  //if we do not have an offset vector we proceed with the standard
  //operation, if we have one given, we add a master result...
  
  //we add a new value for our map
  std::pair<InfoPtr,ResPtr> newInfo;
  newInfo.first.reset(new ExtendedResultInfo);
  //leave the second part as null for the moment, this will be done in finalize
  //just set the name
  newInfo.first->resultName = name;
  newInfo.first->generatorId = filterTag;
  newInfo.first->minStepOffset = minStepOffset;
  newInfo.first->maxStepOffset = maxStepOffset;
  newInfo.first->masterId = generated;
  
  resultMap_[generated] = newInfo;
  return generated;
}

uuids::uuid ResultManager::GetMasterResult(uuids::uuid resultId) {
  return resultMap_[resultId].first->masterId;
}


void ResultManager::CombineResults(uuids::uuid masterId, uuids::uuid slaveId) {
  InfoPtr masterInfo = resultMap_[masterId].first;
  ResPtr masterResult = resultMap_[masterId].second;
  
  slaveId = GetMasterResult(slaveId);
  InfoPtr slaveInfo = resultMap_[slaveId].first;
  
  if (slaveInfo->slaveIds.size() > 0) {
    std::set<uuids::uuid>& slaveIds = slaveInfo->slaveIds;
    std::set<uuids::uuid>::iterator slavItr = slaveIds.begin();
    for (; slavItr != slaveIds.end(); slavItr++) {
      uuids::uuid iSlaveId = *slavItr;
      masterInfo->slaveIds.insert(iSlaveId);
      resultMap_[iSlaveId] = resultMap_[masterId];
    }
  }
  if (slaveInfo->isOutput) {
    masterInfo->isOutput = true;
  }
  masterInfo->minStepOffset = std::min(masterInfo->minStepOffset,slaveInfo->minStepOffset);
  masterInfo->maxStepOffset = std::max(masterInfo->maxStepOffset,slaveInfo->maxStepOffset);
  masterInfo->slaveIds.insert(slaveId);
  resultMap_[slaveId] = resultMap_[masterId];
}

bool ResultManager::IsResultVecUpToDate(uuids::uuid requestedId){
  checkFinalized();
  return resultMap_[requestedId].second->IsUpToDate();
}

void ResultManager::SetResultVecUpToDate(uuids::uuid requestedId, bool upToDate) {
  checkFinalized();
  resultMap_[requestedId].second->SetUpToDate(upToDate);
}


template<typename T>
CF::Vector<T>&  ResultManager::GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers){
  //performance tests have reveals that all syntactic sugar
  //e.g. by providing a view on the vector or anything have not advantage
  //in runtime over just returning an index array which can be used by the caller
  //in addition, the memory requirements for an additional integer array are quite low
  //quite annoying nonetheless, if the requested result is defined on all entites,
  // a vector from 0...N is returned... bullshit
  checkFinalized();
  ResultCache<T> * resCache = dynamic_cast<ResultCache<T> * >(resultMap_[requestedId].second.get());
  eqnNumbers = (*resultMap_[requestedId].first->eqnNumbers.get());
  //return resAdapt->resultVector;
  return resCache->GetResultVector();
}

template<typename T>
CF::Vector<T>&  ResultManager::GetResultVector(uuids::uuid requestedId){
  //performance tests have reveals that all syntactic sugar
  //e.g. by providing a view on the vector or anything have not advantage
  //in runtime over just returning an index array which can be used by the caller
  //in addition, the memory requirements for an additional integer array are quite low
  //quite annoying nonetheless, if the requested result is defined on all entites,
  // a vector from 0...N is returned... bullshit
  checkFinalized();
  ResultCache<T> * resCache = dynamic_cast<ResultCache<T> * >(resultMap_[requestedId].second.get());
  //return resAdapt->resultVector;
  return resCache->GetResultVector();
}

//CF::Vector<Complex>& ResultManager::GetResultVector(uuids::uuid requestedId, CF::StdVector<UInt>& eqnNumbers){
//  checkFinalized();
//  ResultAdaptor<Complex> * resAdapt = dynamic_cast<ResultAdaptor<Complex> * >(resultMap_[remappedId].second.get());
//  //retVec = dynamic_cast<CF::Vector<Complex> & >(resultMap_[requestedId].second->GetSingleVector());
//  eqnNumbers = (*resultMap_[requestedId].first->eqnNumbers.get());
//  return resAdapt->resultVector;
//}

shared_ptr<EqnMapSimple> ResultManager::GetEqnMap(uuids::uuid requestedId) {
  checkFinalized();
  return resultMap_[requestedId].second->mapping;
}

CF::StdVector<shared_ptr<BaseResult> >& ResultManager::GetBaseResultVector(uuids::uuid requestedId){
  checkFinalized();
  return resultMap_[requestedId].second->baseResultVector;
}

UInt ResultManager::GetStepIndex(uuids::uuid requestedId){
  checkFinalized();
  return resultMap_[requestedId].second->GetStepIndex();
}

void ResultManager::SetStepIndex(uuids::uuid requestedId, Integer stepIndex){
  checkFinalized();
  const CF::StdVector<Double>& timeLine = *resultMap_[requestedId].first->timeLine.get();
  UInt setStepIndex = std::min<UInt>(std::max<Integer>(0,stepIndex)
                                     ,timeLine.GetSize() - 1);

  resultMap_[requestedId].second->SetStepIndex(setStepIndex);
}

Double ResultManager::GetStepValue(uuids::uuid requestedId){
  checkFinalized();
  const CF::StdVector<Double>& timeLine = *resultMap_[requestedId].first->timeLine.get();
  const UInt index = resultMap_[requestedId].second->GetStepIndex();
  return timeLine[index];
}

void ResultManager::SetStepValue(uuids::uuid requestedId, Double value){
  checkFinalized();
  const CF::StdVector<Double>& timeLine = *resultMap_[requestedId].first->timeLine.get();
  const UInt size = timeLine.GetSize();
  Double minDist = std::abs(value - timeLine[0]);
  UInt minIndex = 0;
  for (UInt i = 1; i < size; i++) {
    Double dist = std::abs(value - timeLine[i]);
    if (dist < minDist) {
      minDist = dist;
      minIndex = i;
    }
  }
  resultMap_[requestedId].second->SetStepIndex(minIndex);
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
    EXCEPTION("Requested Id which is not generated! Make sure every call "
                  "to ResultManager provides valid resultIds.");
  }

  return resultMap_[resId].first;
}

//has anybody an idea how to avoid this piece of code?
//and do the same without templates?
#define RESULT_MANAGER_COPY_FIELD(resId, fName , toCopy)                                   \
    InfoPtr curInfo = resultMap_[resId].first;                                             \
    curInfo->fName = toCopy;                                                               \

#define RESULT_MANAGER_COPY_SH_PTR_FIELD(resId, fName , toCopy)                            \
    (*resultMap_[resId].first->fName.get()) = toCopy;                                      \

#define RESULT_MANAGER_COPY_SH_PTR_SET(resId, fName , toCopy)                              \
    resultMap_[resId].first->fName->clear();                                               \
    resultMap_[resId].first->fName->insert(toCopy.begin(),toCopy.end());                   \

#define RESULT_MANAGER_OBTAIN_FIELD(resId, fName)                \
    return resultMap_[resId].first->fName;                       \


#define RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId, fName)         \
    return (*resultMap_[resId].first->fName.get());              \

void ResultManager::CopyResultData(uuids::uuid srcId, uuids::uuid trgId){
  checkNotFinalized();
  InfoPtr srcInfo = resultMap_[srcId].first;
  InfoPtr trgInfo = resultMap_[trgId].first;
  
  resultMap_[trgId].first->dType = srcInfo->dType;
  resultMap_[trgId].first->complexFormat = srcInfo->complexFormat;
  resultMap_[trgId].first->definedOn = srcInfo->definedOn;
  resultMap_[trgId].first->dofNames = srcInfo->dofNames;
  resultMap_[trgId].first->entryType = srcInfo->entryType;
  resultMap_[trgId].first->ptGrid = srcInfo->ptGrid;
  resultMap_[trgId].first->isMeshResult = srcInfo->isMeshResult;
  
  (*resultMap_[trgId].first->regNames.get())  = (*srcInfo->regNames.get());
  (*resultMap_[trgId].first->entityNumbers.get()) = (*srcInfo->entityNumbers.get());
  (*resultMap_[trgId].first->eqnNumbers.get()) = (*srcInfo->eqnNumbers.get());
  (*resultMap_[trgId].first->timeLine.get()) = (*srcInfo->timeLine.get());
  (*resultMap_[trgId].first->stepNumbers.get()) = (*srcInfo->stepNumbers.get());
}

CF::StdVector<Double> ResultManager::GetTimeLine(uuids::uuid resId){
  RESULT_MANAGER_OBTAIN_SH_PTR_FIELD(resId,timeLine)
}

bool ResultManager::IsStatic(uuids::uuid resId) {
  return resultMap_[resId].first->isStatic;
}

void ResultManager::SetStatic(uuids::uuid resId, bool isStatic) {
  resultMap_[resId].first->isStatic = isStatic;
}
/* Using IsStatic() due to merge from tuWien
bool ResultManager::IsConstant(uuids::uuid resId) {
  return resultMap_[resId].first->timeLine->GetSize() <= 1;
}
*/
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
  resultMap_[resId].first->isValid = true;
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
  // 2. Which data-structures are necessary?
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
  
  /**

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
**/
//  uuidIter = resultMap_.begin();
//  for(;uuidIter != resultMap_.end();++uuidIter){
//    print_ExtInfoFields((*uuidIter->second.first.get()));
//    //std::cout <<uuidIter->second.second->mapping << std::endl;
//    std::cout << "--------------------------------------------------" << std::endl;
//  }

  // create the result cache for each uniqe result (master result) if valid
  UuidMap::iterator uuidIter = resultMap_.begin();
  std::set<uuids::uuid> toRemove;
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    if(!cInfo->isValid){
      WARN("Result " <<  cInfo->resultName << " is marked invalid after initialization phase. Going to remove it from pool!");
      toRemove.insert(uuidIter->first);
      continue;
    }
    if (cInfo->masterId == uuidIter->first) {
      UInt cacheSteps = 1 + cInfo->maxStepOffset - cInfo->minStepOffset;
      if(cInfo->dType == ExtendedResultInfo::COMPLEX){
        uuidIter->second.second = shared_ptr<GenericResultCache>(new ResultCache<CF::Complex>(cacheSteps, cInfo->isStatic));
      } else if(cInfo->dType == ExtendedResultInfo::DOUBLE){
        uuidIter->second.second = shared_ptr<GenericResultCache>(new ResultCache<CF::Double>(cacheSteps, cInfo->isStatic));
      }else{
        EXCEPTION("Only Complex and Double results supported yet. Result: " << cInfo->resultName)
      }
    }
  }
  //remove invalid results
  std::set<uuids::uuid>::iterator aiter = toRemove.begin();
  for(;aiter != toRemove.end();++aiter){
    resultMap_.erase(*aiter);
    activeIds_.erase(*aiter);
  }
  uuidIter = resultMap_.begin();
  // setting results for slave result ids
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    if (cInfo->masterId != uuidIter->first) {
      uuidIter->second.second = resultMap_[cInfo->masterId].second;
    }
  }
  //get unique grid pointers
//
//  std::cerr << allInfos.size() << " unique resultInfos are to cover. ";
//  std::cerr << allGrids.size() << " unique grids are to cover. ";



  //start with the first result and obtain a vector of equations
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    ResPtr cRes = uuidIter->second.second;
    cRes->SetUpToDate(false);
  }
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    ResPtr cRes = uuidIter->second.second;
    if (cInfo->masterId == uuidIter->first && !cRes->IsUpToDate()) {
    //if(!cRes->IsUpToDate()){
    //if (true) {
      CreateEqnMapping(cInfo,cRes);
      SetResultVectorSize(cInfo,cRes);
      cRes->SetUpToDate(true);
      UuidMap::iterator otherIter = resultMap_.begin();
      //now we loop over all results and try to apply same pointers as possible
      for(;otherIter != resultMap_.end();++otherIter){
        InfoPtr coInfo = otherIter->second.first;
        ResPtr coRes = otherIter->second.second;
        if (coInfo->masterId == uuidIter->first && !coRes->IsUpToDate()) {
        //if (false) {
        //if(!coRes->IsUpToDate()){
          //try to copy shared pointer
          if(checkEqualEqnMap(uuidIter->second, otherIter->second)){
            coRes->mapping = cRes->mapping;
            if(checkEqualEqnVec(uuidIter->second, otherIter->second)){
              coInfo->eqnNumbers = cInfo->eqnNumbers;
            }
            SetResultVectorSize(coInfo,coRes);
            coRes->SetUpToDate(true);
          }
        }
      }
    }
  }
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    ResPtr cRes = uuidIter->second.second;
    cRes->SetUpToDate(false);
  }

  //create output results for the remainder
  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    InfoPtr cInfo = uuidIter->second.first;
    ResPtr cRes = uuidIter->second.second;
    if (cInfo->masterId != uuidIter->first || !cInfo->isOutput) {
      continue;
    }

    //no offset results here...
    StdVector< shared_ptr<CF::BaseResult> >& resVec = cRes->baseResultVector;

    if(cInfo->isMeshResult){
      //this loop body is adapted for mesh results. no matter what we have in entitynumbers from
      // upstream, we write full regions
      std::set<std::string>::iterator regIter = cInfo->regNames->begin();
      for(;regIter != cInfo->regNames->end();++regIter){
        UInt idx = resVec.GetSize();

        if(cInfo->dType == ExtendedResultInfo::COMPLEX){
          resVec.Push_back(shared_ptr< BaseResult >(new Result<CF::Complex>() ) );
        }else{
          resVec.Push_back(shared_ptr< BaseResult >(new Result<CF::Double>() ) );
        }

        shared_ptr<EntityList> aList;
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
        resVec.Push_back(shared_ptr< BaseResult >(new Result<CF::Complex>() ) );
      }else{
        resVec.Push_back(shared_ptr< BaseResult >(new Result<CF::Double>() ) );
      }
      resVec[idx]->SetResultInfo(cInfo->GetResultInfo());
    }
  }

  uuidIter = resultMap_.begin();
  for(;uuidIter != resultMap_.end();++uuidIter){
    uuidIter->second.second->SetUpToDate(false);
  }
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

void ResultManager::SetResultVectorSize(InfoPtr cInfo, ResPtr cRes){
  //first we create a Result
  //obtain initial equation vector
  if(cInfo->entityNumbers->GetSize() == 0){
    //we just resize the solution vector to the global number of equations
    cRes->SetVectorSize(cRes->mapping->GetNumEquations());
  }else{
    //apparently we only need a subset of equations, so we try to cope with that
    cRes->SetVectorSize(cInfo->eqnNumbers->GetSize());
  }
}

template CF::Vector<Double>&  ResultManager::GetResultVector(uuids::uuid, CF::StdVector<UInt>&);
template CF::Vector<Complex>&  ResultManager::GetResultVector(uuids::uuid, CF::StdVector<UInt>&);
template CF::Vector<Double>&  ResultManager::GetResultVector(uuids::uuid);
template CF::Vector<Complex>&  ResultManager::GetResultVector(uuids::uuid);
}

