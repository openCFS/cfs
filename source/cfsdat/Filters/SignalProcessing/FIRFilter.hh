// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FIRFilter.hh
 *       \brief    <Description>
 *
 *       \date     Mar 12, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef FIRFILTER_HH_
#define FIRFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{
class FIRFilter : public BaseFilter{
public:

  //!  Constructor.
  FIRFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~FIRFilter(){
  };
  
  bool UpdateResults(std::set<uuids::uuid>& upResults);

protected:

  bool upsampling_;
  bool downsampling_;
  UInt resampleFactor_;
  
  StdVector<Double> coefficients_;
  StdVector<Integer> coefficientOffsets_;
  
  Integer minOffset_;
  Integer maxOffset_;
  
  virtual void SetDefaultCoefficients();
           
  virtual void ApplyFIRFilter(Vector<Double>& returnVec, uuids::uuid upRes, Integer actStepIndex, 
          StdVector<Double>& coefficients, StdVector<Integer>& coefficientOffsets);

  //! Setup of Results of the filter
  virtual ResultIdList SetUpstreamResults();

  //! Get the maximum step offset from downstream results
  //! \param (in) downStreamResultIds ids of results needed from downstream filter
  //! \return maximum step offset from downstream results
  virtual Integer GetDownStreamMaxStepOffset(std::set<uuids::uuid> downStreamResultIds);
  
  virtual void CopyTimeLineUpstream(uuids::uuid upStreamId, uuids::uuid downStreamId);
  
  void UpsampleTimeLine(StdVector<Double>& inputTimeLine, StdVector<Double>& upSampleTimeLine);
  
  void DownsampleTimeLine(StdVector<Double>& inputTimeLine, StdVector<Double>& downSampleTimeLine);
  
  virtual void AdaptFilterResults();

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

private:

};

}

#endif /* FIRFILTER_HH_ */
