// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Cell2NodeInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Apr 18, 2016
 *       \author   sschoder
 */
//================================================================================================

#pragma once

#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"

namespace CFSDat{

//! Class for choosing a very simple but fast interpolation approach
//! Upon initialization we determine cell NearesNeighbours
//! additionally we set the local coordinates accoring to source
//! during traversal, we just apply those loads
class Cell2NodeInterpolator : public MeshFilter{

public:

  Cell2NodeInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~Cell2NodeInterpolator();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();



private:


  std::vector<QuantityStruct> interpolData_;
  StdVector<UInt> nodeNeighbours_;

  //! Global Factor for scaling the result
  Double globalFactor_;

};

}
