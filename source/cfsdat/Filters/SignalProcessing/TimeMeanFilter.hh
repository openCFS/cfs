// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TimeMeanFilter.hh
 *       \brief    <Description>
 *
 *       \date     Apr 1, 2017
 *       \author   sschoder
 */
//================================================================================================

#ifndef TIMEMEANFILTER_HH_
#define TIMEMEANFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{
class TimeMeanFilter : public BaseFilter{
public:
  struct tMeanInfo{
    CoupledField::SolutionType d1Type;
    CoupledField::SolutionType d2Type;
  };

  typedef std::map<CoupledField::SolutionType, tMeanInfo> TimeMeanMap;
  static TimeMeanMap tMeanMap_;

  //!  Constructor.
  TimeMeanFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
    : BaseFilter(numWorkers,config, resMan){
    baseQuantity_ = CoupledField::NO_SOLUTION_TYPE;
    isFirstRun_ = true;
    filtStreamType_ = FIFO_FILTER;
  }

  virtual ~TimeMeanFilter(){

  }


protected:

  //! Setup of Results of the filter
  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

  CoupledField::SolutionType baseQuantity_;


  //maps filterResultIds to timestep values
  std::map<uuids::uuid, Double> timeSteps_;


private:
  //! standard map for results all tag in xml
  static TimeMeanMap init_mean_map() {
    TimeMeanMap myMap;
    myMap[CoupledField::FLUIDMECH_PRESSURE].d1Type = CoupledField::FLUIDMECH_PRESSURE_DERIV_1;
    myMap[CoupledField::FLUIDMECH_PRESSURE].d2Type = CoupledField::FLUIDMECH_PRESSURE_DERIV_2;
    myMap[CoupledField::FLUIDMECH_VELOCITY].d1Type = CoupledField::FLUIDMECH_VELOCITY_DERIV_1;
    myMap[CoupledField::FLUIDMECH_VELOCITY].d2Type = CoupledField::FLUIDMECH_VELOCITY_DERIV_2;
    return myMap;
  }

  bool isFirstRun_;
};

class TimeMeanFilter1 : public TimeMeanFilter{

public:


  TimeMeanFilter1(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~TimeMeanFilter1(){

  }

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

  struct filtRes  {};
  struct upRes    {};

  typedef boost::bimaps::bimap
  <
    boost::bimaps::tagged< std::string , filtRes >,
    boost::bimaps::tagged< std::string , upRes   >
  > InOutBiMap;

 //stores association of upstream result name to downstream name
  InOutBiMap inOutNames_;

 //stores for each filter result the dependend upstream ids
 std::map<std::string, uuids::uuid> filtResToUpResIds_;

// Vector<Double>& returnVec;
 Vector<Double> tmp_;
 Integer N_;

};
}
#endif /* DERIVFILTERD1_HH_ */
