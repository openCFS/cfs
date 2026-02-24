// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BinOpFilter.hh
 *       \brief    <Description>
 *
 *       \date     Mai 13, 2018
 *       \author   mtautz
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_
#define SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include "BinOpFunctions.hh"
#include <boost/bimap.hpp>

namespace CFSDat{

class BinOpFilter : public BaseFilter{
  
public:

  BinOpFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~BinOpFilter();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();
  
  virtual void PrepareCalculation();
  
  std::string opType;
  std::string multType;

  std::string resAName;
  uuids::uuid resAId;
  bool isConstantA;
  Vector<Double> constantA;

  std::string resBName;
  uuids::uuid resBId;
  bool isConstantB;
  Vector<Double> constantB;

  std::string outName;
  uuids::uuid outId;
  
  BinOpFunctions::BinOpFctStruct<Double>::BinOpFctPtr applyFctPtr;
  UInt size;

};

}

#endif /* SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_ */
