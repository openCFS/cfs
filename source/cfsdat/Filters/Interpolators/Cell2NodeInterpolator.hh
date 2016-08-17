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

//#ifndef SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CELL2NODEINTERPOLATOR_HH_
//#define SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CELL2NODEINTERPOLATOR_HH_

#pragma once

#include "MeshBasedInterpolator.hh"
#include "DataInOut/SimInput.hh"

namespace CFSDat{

//! Class for choosing a very simple but fast interpolation approach
//! Upon initialization we determine cell NearesNeighbours
//! additionally we set the local coordinates accoring to source
//! during traversal, we just apply those loads
class Cell2NodeInterpolator : public MeshBasedInterpolator{

  struct InpolationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    UInt srcEqn;
    Double volume;

    InpolationStruct() : tENum(0),srcEqn(0),volume(.0){
      localCoords.Resize(3);
    }

    bool operator < (const InpolationStruct& str) const
    {
        return (srcEqn < str.srcEqn);
    }
  };

public:

  Cell2NodeInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~Cell2NodeInterpolator();

  virtual bool Run();

protected:

  virtual void PrepareInterpolation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:


  std::vector<InpolationStruct> interpolData_;
  StdVector<UInt> nodeNeighbours_;

};

}
//#endif /* SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CELL2NODEINTERPOLATOR_HH_ */
