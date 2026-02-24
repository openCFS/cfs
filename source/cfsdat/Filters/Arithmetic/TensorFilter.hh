// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BinOpFilter.hh
 *       \brief    <Description>
 *
 *       \date     Apr 13, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_TENSOR_HH_
#define SOURCE_CFSDAT_FILTERS_TENSOR_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{

class TensorFilter : public BaseFilter{
  
public:

  TensorFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~TensorFilter();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();
  
  virtual void PrepareCalculation();
  
  StdVector<std::string> resNames;
  StdVector<uuids::uuid> resIds;

  std::string outName;
  uuids::uuid outId;
  
  UInt size;

};

}

#endif /* SOURCE_CFSDAT_FILTERS_TENSOR_HH_ */
