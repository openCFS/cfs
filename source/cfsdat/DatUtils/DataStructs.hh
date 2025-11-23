// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DataStructs.hh
 *       \brief    <Description>
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef DATASTRUCTS_HH_
#define DATASTRUCTS_HH_

#include "General/Environment.hh"
#include "Domain/Results/BaseResults.hh"
#include "cfsdat/DatUtils/Defines.hh"
#include "cfsdat/DatUtils/ResultCache.hh"
#include "EqnNumberingSimple.hh"

#include <boost/fusion/include/sequence.hpp>
#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/vector.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/bind/bind.hpp>
#include <boost/fusion/adapted/struct/adapt_struct.hpp>
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include <typeinfo>
#if !defined(WIN32)
  #include <cxxabi.h>
#endif

using namespace boost::placeholders;

namespace CFSDat{




struct ExtendedResultInfo : public CF::ResultInfo, enable_shared_from_this<ExtendedResultInfo>{

  friend class ResultManager;

  typedef enum {UNDEF_DTYPE, INT, DOUBLE, COMPLEX } ResDType;

  ExtendedResultInfo(){
    generatorId = uuids::nil_uuid();
    dType = UNDEF_DTYPE;
    isOutput = false;
    sequenceStep = 0;
    timeLine.reset(new CF::StdVector<Double>);
    regNames.reset(new std::set<std::string>);
    entityNumbers.reset(new CF::StdVector<UInt>);
    eqnNumbers.reset(new CF::StdVector<UInt>);
    stepNumbers.reset(new CF::StdVector<UInt>);
    isValid = false;
    isStatic = false;
    ptGrid = NULL;
    isMeshResult = false;
    minStepOffset = 0;
    maxStepOffset = 0;
    masterId = uuids::nil_uuid();
  }

  str1::shared_ptr<CF::ResultInfo> GetResultInfo(){
    return str1::shared_ptr<CF::ResultInfo>(this);
  }

  uuids::uuid generatorId;

  UInt sequenceStep;

  str1::shared_ptr< CF::StdVector<Double> > timeLine;

  str1::shared_ptr< CF::StdVector<UInt> > stepNumbers;

  str1::shared_ptr< std::set<std::string> > regNames;

  str1::shared_ptr< CF::StdVector<UInt> > entityNumbers;

  str1::shared_ptr< CF::StdVector<UInt> > eqnNumbers;

  ResDType dType;

  bool isOutput;

  bool isMeshResult;

  bool isValid;

  Grid* ptGrid;
  
  Integer minStepOffset;
  
  Integer maxStepOffset;
  
  uuids::uuid masterId;
  
  std::set<uuids::uuid> slaveIds;

  bool operator==(const ExtendedResultInfo& other){
    bool isEqual = true;
    isEqual &= other.dType == this->dType;
    isEqual &= other.regNames == this->regNames;
    isEqual &= (*other.entityNumbers.get()) == (*this->entityNumbers.get());
    isEqual &= (*other.eqnNumbers.get()) == (*this->eqnNumbers.get());
    isEqual &= (*other.timeLine.get()) == (*this->timeLine.get());
    isEqual &= (*other.stepNumbers.get()) == (*this->stepNumbers.get());
    isEqual &= other.complexFormat == this->complexFormat;
    isEqual &= other.definedOn == this->definedOn;
    isEqual &= other.isStatic == this->isStatic;
    isEqual &= other.dofNames == this->dofNames;
    isEqual &= other.entryType == this->entryType;
    isEqual &= other.resultName == this->resultName;
    isEqual &= other.unit == this->unit;
    return isEqual;
  }

  void operator=(str1::shared_ptr<ExtendedResultInfo> other){
    this->dType               = other->dType;
    this->regNames            = other->regNames;
    this->entityNumbers       = other->entityNumbers;
    this->eqnNumbers          = other->eqnNumbers;
    this->timeLine            = other->timeLine;
    this->stepNumbers         = other->stepNumbers;
    this->complexFormat       = other->complexFormat;
    this->definedOn           = other->definedOn;
    this->isStatic            = other->isStatic;
    this->dofNames            = other->dofNames;
    this->entryType           = other->entryType;
    this->resultName          = other->resultName;
    this->unit                = other->unit;
    this->minStepOffset       = other->minStepOffset;
    this->maxStepOffset       = other->maxStepOffset;
    this->masterId            = other->masterId;
    this->slaveIds            = other->slaveIds;
  }

  void ImportResultInfo(str1::shared_ptr<ResultInfo> info){
    this->definedOn           = info->definedOn;
    this->isStatic            = info->isStatic;
    this->dofNames            = info->dofNames;
    this->entryType           = info->entryType;
    this->resultName          = info->resultName;
    this->resultType          = info->resultType;
    this->unit                = info->unit;
    this->complexFormat       = info->complexFormat;
  }



private:

  void CopyInfoFrom(str1::shared_ptr<ExtendedResultInfo> other){
    this->dType               = other->dType;
    this->regNames            = other->regNames;
    this->entityNumbers       = other->entityNumbers;
    this->eqnNumbers          = other->eqnNumbers;
    this->timeLine            = other->timeLine;
    this->stepNumbers         = other->stepNumbers;
    this->complexFormat       = other->complexFormat;
    this->definedOn           = other->definedOn;
    this->isStatic            = other->isStatic;
    this->dofNames            = other->dofNames;
    this->entryType           = other->entryType;
    this->resultName          = other->resultName;
    this->unit                = other->unit;
    this->minStepOffset       = other->minStepOffset;
    this->maxStepOffset       = other->maxStepOffset;
    this->masterId            = other->masterId;
    this->slaveIds            = other->slaveIds;
  }

};

}

BOOST_FUSION_ADAPT_STRUCT(
    CFSDat::ExtendedResultInfo,
    (std::string, resultName)
    (CFSDat::ExtendedResultInfo::ResDType, dType)
    (bool, isOutput)
    (CFSDat::str1::shared_ptr< CoupledField::StdVector<Double> >,  timeLine)
    (CFSDat::str1::shared_ptr< CoupledField::StdVector<UInt> >,  stepNumbers)
    (CFSDat::str1::shared_ptr< std::set<std::string> >,  regNames)
    (CFSDat::str1::shared_ptr< CoupledField::StdVector<UInt> >,  entityNumbers)
    (CFSDat::str1::shared_ptr< CoupledField::StdVector<UInt> >,  eqnNumbers)
    (CFSDat::UInt, sequenceStep )
    (CoupledField::ResultInfo::EntryType , entryType)
    (CoupledField::StdVector<std::string>, dofNames)
    (std::string, unit)
    (CoupledField::ResultInfo::EntityUnknownType, definedOn)
    (bool, isValid)
    (bool, isStatic)
    (Grid* , ptGrid)
    (bool, isMeshResult)
)

namespace CFSDat{

struct printConstField
{
    template <class Index, typename T>
    void operator()(Index,const T& t) const {
      Integer status;
      char *realname = abi::__cxa_demangle(typeid(boost::fusion::at<Index>(t)).name() , 0, 0, &status);
      std::string name = boost::fusion::extension::struct_member_name<T ,Index::value>::call();
      std::cout << realname << "  " << name  << " : " << boost::lexical_cast<std::string>( boost::fusion::at<Index>(t) );
      std::cout << std::endl;
    }
};


template<class C>
void print_ConstExtInfoFields(const C & c)
{
    typedef boost::mpl::range_c<int,0, boost::fusion::result_of::size<C>::type::value> range;
    boost::mpl::for_each<range>(boost::bind<void>(printConstField(), _1, boost::ref(c)));
}

struct printField
{
    template <class Index, typename T>
    void operator()(Index,const T& t) const {
      Integer status;
      char *realname = abi::__cxa_demangle(typeid(boost::fusion::at<Index>(t)).name() , 0, 0, &status);
      std::string name = boost::fusion::extension::struct_member_name<T ,Index::value>::call();
      std::cout << realname << "  " << name  << " : " << boost::lexical_cast<std::string>( boost::fusion::at<Index>(t) );
      std::cout << std::endl;
    }
};


template<class C>
void print_ExtInfoFields(const C & c)
{
    std::cout << std::endl << "--------------------------------------------------" << std::endl;
    std::cout << "Printing Result info Struct for " << c.resultName << std::endl;
    typedef boost::mpl::range_c<int,0, boost::fusion::result_of::size<C>::type::value> range;
    boost::mpl::for_each<range>(boost::bind<void>(printField(), _1, boost::ref(c)));
}

}
//(std::string, boost::lexical_cast<std::string>(generatorId))

//class ResultAdaptor{
//
//public:
//  enum {UNDEF_DTYPE, INT, REAL, COMPLEX } ResDType;
//
//  enum {UNDEF_ENTITIES, NODE, ELEMENT, REGION, SURFREGION } ResDefOn;
//
//  enum {UNDEF_TIMEVALS, ALL_TIMES, SOME_TIMES } DefTime;
//
//  ResultAdaptor()
//    : stepVal_(0.0){
//    solType_ = CF::NO_SOLUTION_TYPE;
//    solName_ = "";
//    isUpToDate_ = false;
//  }
//
//  void SetStepValue(Double timeFreq){
//    if(timeFreq != stepVal_){
//      stepVal_ = timeFreq;
//      isUpToDate_ = false;
//    }
//  }
//
//  void SetStepVector(Double timeFreq){
//    if(timeFreq != stepVal_){
//      stepVal_ = timeFreq;
//      isUpToDate_ = false;
//    }
//  }
//
//  shared_ptr< CF::BaseResult > GetResult(){
//    //this is not very nice because we do not know if and how the result is updated
//    return currentResult_;
//  }
//
//  bool IsUpToDate(){
//    return isUpToDate_;
//  }
//
//  UInt GetStepVal(){
//    return stepVal_;
//  }
//
//  //can we make his a const pointer?
//  boost::shared_ptr<EqnMapSimple> GetEqnMap(){
//    return mapping_;
//  }
//
//  void SetEqnMapping(shared_ptr<EqnMapSimple> map){
//    mapping_ = map;
//  }
//
//private:
//
//  //if possible, store CFS solution type
//  CF::SolutionType solType_;
//
//  shared_ptr< CF::BaseResult > currentResult_;
//
//  Double stepVal_;
//
//  boost::shared_ptr<EqnMapSimple> mapping_;
//
//  bool isUpToDate_;
//};
//
////typedef CF::StdVector<boost::shared_ptr<DataStruct> > ResultList;
//
//
//struct InOutInfo{
//
//  InOutInfo(){
//    resDType = ARB_DTYPE;
//    definedOn = ARB_ENTITIES;
//    solType = CoupledField::NO_SOLUTION_TYPE;
//    numPreviousSteps = 0;
//    numFutureSteps = 0;
//    allTimes = true;
//  }
//
//  enum {ARB_DTYPE, INT, REAL, COMPLEX } resDType;
//
//
//  enum {ARB_ENTITIES, NODE, ELEMENT, REGION, SURFREGION } definedOn;
//
//  ///solution type as a shortcut
//  CF::SolutionType solType;
//
//  UInt               numPreviousSteps;   //do we need information from the past?
//  UInt               numFutureSteps;      //do we need information from the future?
//
//  bool allTimes;
//
//  CF::StdVector<Double> timeValues;
//
//  CF::StdVector<UInt> timeSteps;
//
//  inline bool operator == (const InOutInfo &b) const {
//    bool chk = this->definedOn == b.definedOn;
//    chk = chk && this->resDType == b.resDType;
//    chk = chk && this->numPreviousSteps == b.numPreviousSteps;
//    chk = chk && this->numFutureSteps == b.numFutureSteps;
//
//    if(this->allTimes == b.allTimes){
//      return chk;
//    }else{
//      if(this->timeSteps.GetSize()!=b.timeSteps.GetSize())
//        return false;
//      if(this->timeValues.GetSize()!=b.timeValues.GetSize())
//        return false;
//
//      for(CoupledField::UInt i=0;i<timeValues.GetSize();i++){
//        chk = chk && this->timeValues[i]==b.timeValues[i];
//        chk = chk && this->timeSteps[i]==b.timeSteps[i];
//      }
//      return chk;
//    }
//  }
//
//};
//
////maybe a switch to a std::set as underlying datastructure would be more appropriate
//class InOutInfoList{
//public:
//  InOutInfoList(){
//
//  }
//  ~InOutInfoList(){
//
//  }
//
//  ///returns a subset of the vector containing a certain solType
//  InOutInfoList Get(CoupledField::SolutionType solType){
//    InOutInfoList newList;
//    for(UInt i=0;i<infos_.GetSize();++i){
//      if(infos_[i].solType == solType){
//        newList.Add(infos_[i]);
//      }
//    }
//    return newList;
//  }
//
//  void Add(InOutInfo newInfo){
//    if(infos_.Find(newInfo) == -1)
//      infos_.Push_back(newInfo);
//    return;
//  }
//
//  void Remove(InOutInfo info){
//    infos_.Erase(infos_.Find(info));
//  }
//
//  UInt GetSize(){
//    return infos_.GetSize();
//  }
//
//private:
//  CF::StdVector<InOutInfo> infos_;
//};


#endif /* DATASTRUCTS_HH_ */
