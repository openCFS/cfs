// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     PostLighthillSource.hh
 *       \brief    <Description>
 *
 *       \date     Feb 16, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include "cfsdat/Filters/BaseFilter.hh"

namespace CFSDat{

//! Simple filter to extract the physically scalar acoustic-source-term from the
//! 2D storage format. We just have to extract the x-component, because the source
//! term filter writes on the x-component (y is zero)

class PostLighthillSource : public BaseFilter{

public:
  PostLighthillSource(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~PostLighthillSource();


protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

  int component_;

};



}
