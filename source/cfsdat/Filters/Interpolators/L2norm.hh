// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     L2norm.hh
 *       \brief    <Description>
 *
 *       \date     Feb 13, 2017
 *       \author   sschoder
 */
//================================================================================================

#pragma once


#include "DataInOut/SimInput.hh"
#include <cfsdat/DatUtils/Point.hh>
#include <Filters/MeshFilter.hh>


namespace CFSDat{


class L2norm : public MeshFilter{

public:
  struct InpolationStruct{
    CF::StdVector<UInt> tNNum; //target node numbers
    CF::StdVector<UInt> srcEqn; //source equations (nodes of target-element)
    UInt trgElemNum;
    Double volume;
    InpolationStruct() : trgElemNum(0),volume(.0){
      tNNum.Resize(1);
      srcEqn.Resize(0);
    }
  };

  L2norm(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~L2norm();


protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

    //! Vector containing target-informations
    std::vector<InpolationStruct> interpolData_;

};

}
