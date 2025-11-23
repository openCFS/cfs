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
#include <cfsdat/DatUtils/DataStructs.hh>
#include <boost/shared_ptr.hpp>
// https://github.com/boostorg/uuid/issues/91 applies for the gitlabrunners with cfsdat
#define BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "cfsdat/DatUtils/DataStructs.hh"
#include "cfsdat/DatUtils/ResultManager.hh"
#include "cfsdat/DatUtils/Defines.hh"
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

  //! Run the Filter to calculate actually used results.
  //! This function may be overwritten by implementing filter,
  //! alternatively, the standard implementation of BaseFilter 
  //! can be used which calls UpdateResults() if needed.
  //! The standard implementation extracts the results of the
  //! specific filter which need recomputation, 
  //! calls the UpdateResults() function, 
  //! sets the computed results UpToDate,
  //! and deactivates all upstream results
  //! \return success of the computation
  virtual bool Run();
  
  void FinishInit();
  
  void InitResults();
  
  bool IsOutput(){
    return filtStreamType_ == OUTPUT_FILTER;
  }

  std::string GetId(){
    return filterId_;
  }

  const std::set<std::string>& GetInputIds(){
    return inputIds_;
  }

protected:

  virtual void AddInput(str1::shared_ptr<BaseFilter> filt){
    this->sources_.Push_back(filt);
  }

  virtual void AddOutput(str1::shared_ptr<BaseFilter> filt){
    this->sinks_.Push_back(filt);
  }

  virtual void PrepareCalculation();

  void InitResultsUpstream();
  
  void InitResultsDownstream();

  //=================================================================
  // Helper Functions for Upstream Result Initialization
  //=================================================================
  
  //! Extracts the filter result ids, to be computed by this filter 
  //! given in filtResNames.
  //! \return Nothing
  virtual void ExtractFilterResults();

  virtual ResultIdList SetUpstreamResults()=0;
  
  bool CheckNeeded();
  
  //! Sets the default results needed from upstream filters according
  //! to the names given in upResNames, while assuming that only one
  //! result is computed by this filter.
  //! Fills the map denoted as upResNameIds.
  //! \return List of upstream result ids
  ResultIdList SetDefaultUpstreamResults();
  
  //! Registers a result needed from upstream filters,
  //! assuming that no time step offsets are requried
  //! \param (in) name name of the result to register
  //! \param (in) downStreamResultId id of result needed from downstream filter
  //! \return success of the computation
  uuids::uuid RegisterUpstreamResult(std::string name, uuids::uuid downStreamResultId);
                         
  //! Registers a result needed from upstream filters,
  //! assuming that no time step offsets are requried
  //! \param (in) name name of the result to register
  //! \param (in) downStreamResultIds ids of result needed from downstream filter
  //! \return success of the computation
  uuids::uuid RegisterUpstreamResult(std::string name, std::set<uuids::uuid> downStreamResultIds);
  
  //! Registers a result needed from upstream filters.
  //! \param (in) name name of the result to register
  //! \param (in) minOffset minimum offset step required relative to the downstream result to be computed
  //! \param (in) maxOffset maximum offset step required relative to the downstream result to be computed
  //! \param (in) downStreamResultId id of result needed from downstream filter
  //! \return success of the computation
  uuids::uuid RegisterUpstreamResult(std::string name, Integer minOffset, 
                         Integer maxOffset, uuids::uuid downStreamResultId);

  //! Registers a result needed from upstream filters.
  //! \param (in) name name of the result to register
  //! \param (in) minOffset minimum offset step required relative to the downstream result to be computed
  //! \param (in) maxOffset maximum offset step required relative to the downstream result to be computed
  //! \param (in) downStreamResultIds ids of results needed from downstream filter
  //! \return success of the computation
  uuids::uuid RegisterUpstreamResult(std::string name, Integer minOffset, 
                         Integer maxOffset, std::set<uuids::uuid> downStreamResultIds);
  
  //! Get the maximum step offset from downstream results
  //! \param (in) downStreamResultIds ids of results needed from downstream filter
  //! \return maximum step offset from downstream results
  virtual Integer GetDownStreamMaxStepOffset(std::set<uuids::uuid> downStreamResultIds);
  
  virtual void CopyTimeLineUpstream(uuids::uuid upStreamId, uuids::uuid downStreamId);
                         
  //=================================================================
  // Helper Functions for Downstream Result Initialization
  //=================================================================
  
  virtual void AdaptFilterResults()=0;

  //=================================================================
  // Helper Functions for Result Calculation in Run Function
  //=================================================================
  
  //! This function is called by the standard implementation 
  //! of the Run() function. It calculated the results of the 
  //! given ids.
  //! \param (in) upResults set of results to be updated
  //! \return success of the computation
  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  //! Extracts the active/needed results which require recomputation
  //! \return set of results which need recomputation
  virtual std::set<uuids::uuid> ExtractObsoleteResults();
  
  void VerboseSum(uuids::uuid verbResId);
  
  
  //! Returns the result vector of a result computed by this filter
  //! \param (in) resId id of the result
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetOwnResultVector(uuids::uuid resId, CF::StdVector<UInt>& eqnNumbers) {
    return resultManager_->GetResultVector<T>(resId, eqnNumbers);
  }
  
  //! Returns the result vector of a result computed by this filter
  //! \param (in) resId id of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetOwnResultVector(uuids::uuid resId) {
    return resultManager_->GetResultVector<T>(resId);
  }
  
  //! Prepares an upstream result by invoking the run functions of
  //! upstream filters, given that the result is not up to date
  //! \return Nothing
  void PrepareUpstreamResult(uuids::uuid resId);
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId, CF::StdVector<UInt>& eqnNumbers) {
    PrepareUpstreamResult(resId);
    return resultManager_->GetResultVector<T>(resId, eqnNumbers);
  }

  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName, CF::StdVector<UInt>& eqnNumbers) {
    return GetUpstreamResultVector<T>(upResNameIds[resName], eqnNumbers);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \param (in) stepValue step value required
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId, Double stepValue, 
                                    CF::StdVector<UInt>& eqnNumbers) {
    resultManager_->SetStepValue(resId,stepValue);
    return GetUpstreamResultVector<T>(resId, eqnNumbers);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \param (in) stepValue step value required
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName, Double stepValue,
                                    CF::StdVector<UInt>& eqnNumbers) {
    return GetUpstreamResultVector<T>(upResNameIds[resName], stepValue, eqnNumbers);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \param (in) stepIndex step index in the timeline
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId, Integer stepIndex,
                                    CF::StdVector<UInt>& eqnNumbers) {
    resultManager_->SetStepIndex(resId,stepIndex);
    return GetUpstreamResultVector<T>(resId, eqnNumbers);
  }

  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \param (in) stepIndex step index in the timeline
  //! \param (out) eqnNumbers equation number of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName, Integer stepIndex, 
                                    CF::StdVector<UInt>& eqnNumbers) {
    return GetUpstreamResultVector<T>(upResNameIds[resName], stepIndex, eqnNumbers);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \param (in) stepValue step value required
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId, Double stepValue) {
    resultManager_->SetStepValue(resId,stepValue);
    return GetUpstreamResultVector<T>(resId);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \param (in) stepValue step value required
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName, Double stepValue) {
    return GetUpstreamResultVector<T>(upResNameIds[resName], stepValue);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \param (in) stepIndex step index in the timeline
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId, Integer stepIndex) {
    resultManager_->SetStepIndex(resId,stepIndex);
    return GetUpstreamResultVector<T>(resId);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \param (in) stepIndex step index in the timeline
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName, Integer stepIndex) {
    return GetUpstreamResultVector<T>(upResNameIds[resName], stepIndex);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resId id of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(uuids::uuid resId) {
    PrepareUpstreamResult(resId);
    return resultManager_->GetResultVector<T>(resId);
  }
  
  //! Function for receiving a result vector computed by upstream filters.
  //! Therefore, this function invokes Run() of the usptream filters
  //! \param (in) resName name of the result
  //! \return a reference to the result vector
  template<typename T>
  Vector<T>& GetUpstreamResultVector(std::string resName) {
    return GetUpstreamResultVector<T>(upResNameIds[resName]);
  }
  
  //! Deactivates all results needed from upstream filters
  //! \return Nothing
  void DeactivateUpstreamResults();

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

  /// list of incomming, active result ids which match the entries in filtResNames
  ResultIdList filterResIds;
  
  /// map from filter result names to filter result Ids
  /// filled by ExtractFilterResults
  std::map<std::string, uuids::uuid> filterResNameIds;
  
  /// set of all result names which are required from upstream filters
  std::set<std::string> upResNames;

  /// map from upstream result names to upstream Ids
  /// filled by RegisterUpstreamResult
  std::map<std::string, uuids::uuid> upResNameIds;
  
  /// list of all result ids which are required for the filter from the upstream filters
  /// w.r.t. the entries in filterResIds
  ResultIdList upResIds;

  

  //@}

  //currently as normal member but it will change to static...
  //! map storing the global step numbers as provided by the user
  //! mostly used by output filters to create requested time map
  //! some filters may ignore this or define their own
  std::map<UInt,Double> globalStepValueMap_;

  StreamType filtStreamType_;

  /// Pointers to upstream filters
  CF::StdVector< str1::shared_ptr<BaseFilter> > sources_;

  /// Pointers to downstream filters
  CF::StdVector< str1::shared_ptr<BaseFilter> > sinks_;
  
  /// counter for InitResults function, should not exceed the number of sources
  UInt initSourceResults_;

  /// counter for InitResults function, should not exceed the number of sinks
  UInt initSinkResults_;
  
  /// counter for FinishInit function, should not exceed the number of sinks
  UInt finishInitSinkResults_;

  UInt numWorkers_;

  CoupledField::PtrParamNode params_;

  std::string filterId_;

  std::set<std::string> inputIds_;

  boost::uuids::uuid filterTag_;

  str1::shared_ptr<ResultManager> resultManager_;

  /// timer for performance usage
  std::map<boost::uuids::uuid , CF::shared_ptr<CF::Timer> > theTimers;

  // Start time to remove offset of data
  Double startTime_ = 0.0;

private:

};

}



#endif /* BASEFILTER_HH_ */
