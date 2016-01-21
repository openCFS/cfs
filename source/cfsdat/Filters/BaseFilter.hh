// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseFilter.hh
 *       \brief    Basic interface for CFSDat filters
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef BASEFILTER_HH_
#define BASEFILTER_HH_

#include <map>
#include <set>
#include <cfsdat/Utils/DataStructs.hh>
#include <boost/shared_ptr.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "cfsdat/Utils/DataStructs.hh"
#include "cfsdat/Utils/ResultManager.hh"
#include "cfsdat/Utils/Defines.hh"
#include "Utils/Timer.hh"

namespace CFSDat{

class BaseFilter{

public:

  typedef enum  {FIFO_FILTER, ACCUMULATIVE_FILTER, OUTPUT_FILTER, INPUT_FILTER, NO_STREAM } StreamType;


  static void Connect(str1::shared_ptr<BaseFilter> source, str1::shared_ptr<BaseFilter> sink){
    sink->AddInput(source);
    source->AddOutput(sink);
  }

  static FilterPtr Generate(PtrParamNode filtNode,PtrResultManager resMana);


  BaseFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~BaseFilter(){

  }

  virtual bool Run() = 0;

  virtual void FinishInit(){
    CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
    for(; srcIter != sources_.End() ; srcIter++){
      // should we check here anything for success?
      (*srcIter)->FinishInit();
    }
  }

  void InitResults();

  bool IsOutput(){
    return filtSteamType_ == OUTPUT_FILTER;
  }

  std::string GetId(){
    return filterId_;
  }

  std::set<std::string> GetInputIds(){
    return inputIds_;
  }

protected:

  virtual void AddInput(str1::shared_ptr<BaseFilter> filt){
    this->sources_.Push_back(filt);
  }

  virtual void AddOutput(str1::shared_ptr<BaseFilter> filt){
    this->sinks_.Push_back(filt);
  }

  virtual void ExtractFilterResults(){

    filterResIds.Clear();

    std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
    std::set<uuids::uuid>::iterator aIter = activeResults.begin();
    for(; aIter != activeResults.end(); ++aIter){
      ResultManager::ConstInfoPtr aInfo = resultManager_->GetExtInfo(*aIter);
      //print_ConstExtInfoFields((*aInfo.get()));

      if(filtResNames.find(aInfo->resultName) != filtResNames.end()){
        filterResIds.Push_back(*aIter);
        resultManager_->DeactivateResult(*aIter);
      }
    }
  }

  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;


  virtual void FillGlobalStepsFromFilter(){
    //not  here but in input filter
  }


  ///Small difference time value callback
  struct time_cmp {
    time_cmp(double v, double d) : val(v), delta(d) { }
      inline bool operator()(const double &x) const {
          return abs(x-val) < delta;
      }
  private:
      double val, delta;
  };


  ///Starts a Timer inside the filter
  inline boost::uuids::uuid StartTime(){
    boost::uuids::uuid tId = boost::uuids::random_generator()();
    this->theTimers[tId].reset(new Timer);
    this->theTimers[tId]->Start();
    return tId;
  }

  ///Stops given Timer and returns the result
  inline std::string StopTime(boost::uuids::uuid tId){
    this->theTimers[tId]->Stop();
    std::stringstream elapsed;
    const int walltime((int) theTimers[tId]->GetWallTime());
    if(walltime > 120) {
      const int wallmin((int) (walltime / 60.0));
      if(wallmin > 60){
        elapsed << wallmin/60 << "h";
      }else{
        elapsed << wallmin << "min";
      }
    }else{
      elapsed << walltime << "s";
    }
    //remove timer from map
    this->theTimers[tId].reset();
    theTimers.erase(tId);

    return elapsed.str();
  }

  //=================================================================
  // Result Storage
  //=================================================================
  //@{ \name Storage for result Ids and strings

  /// set of all result names which can be provided by the filter
  std::set<std::string> filtResNames;

  /// set of all result names which are required from upstream filters
  std::set<std::string> upResNames;

  /// list of all result ids which are required for the filter from the upstream filters
  /// w.r.t. the entries in filterResIds
  ResultIdList upResIds;

  /// list of incomming, active result ids which match the entries in filtResNames
  ResultIdList filterResIds;

  //@}

  //currently as normal member but it will change to static...
  //! map storing the global step numbers as provided by the user
  //! mostly used by output filters to create requested time map
  //! some filters may ignore this or define their own
  std::map<UInt,Double> globalStepValueMap_;

  StreamType filtSteamType_;

  CF::StdVector< str1::shared_ptr<BaseFilter> > sources_;

  CF::StdVector< str1::shared_ptr<BaseFilter> > sinks_;


  UInt numWorkers_;

  CoupledField::PtrParamNode params_;

  std::string filterId_;

  std::set<std::string> inputIds_;

  boost::uuids::uuid filterTag_;

  str1::shared_ptr<ResultManager> resultManager_;

  /// timer for performance usage
  std::map<boost::uuids::uuid , CF::shared_ptr<CF::Timer> > theTimers;


private:

};

}



#endif /* BASEFILTER_HH_ */
