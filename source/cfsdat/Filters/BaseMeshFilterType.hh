// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseMeshFilterType.hh
 *       \brief    <Description>
 *
 *       \date     Jan 11, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "cfsdat/Filters/BaseFilter.hh"


namespace CFSDat{

//! Class for choosing the type of (mesh-based) filter, which shall be created
//! according to the received parameter node pointer
class BaseMeshFilterType : public BaseFilter {

public:

  //! Choose the filter according to ptrNode
  static FilterPtr Generate(PtrParamNode ptrNode, PtrResultManager resMana);

  BaseMeshFilterType(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
    :BaseFilter(numWorkers,config,resMan){

  }

  virtual ~BaseMeshFilterType(){
  }

protected:

  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

};

}
