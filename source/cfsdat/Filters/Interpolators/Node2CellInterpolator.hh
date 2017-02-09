// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Node2CellInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Jan 3, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "DataInOut/SimInput.hh"
#include <cfsdat/Utils/Point.hh>
#include <Filters/MeshFilter.hh>


namespace CFSDat{


class Node2CellInterpolator : public MeshFilter{

public:
  struct InpolationStruct{
    CF::StdVector<UInt> tNNum; //target node numbers
    CF::StdVector<UInt> srcEqn; //source equations (nodes of target-element)
    UInt trgElemNum;
    InpolationStruct() : trgElemNum(0){
      tNNum.Resize(1);
      srcEqn.Resize(0);
    }
  };

  Node2CellInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~Node2CellInterpolator();

  virtual bool Run();



protected:

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

    //! Vector containing target-informations
    std::vector<InpolationStruct> interpolData_;

};

}
