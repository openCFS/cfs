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

#ifndef SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_
#define SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{

class BinOpFilter : public BaseFilter{
public:

  static FilterPtr GenerateOperator(PtrParamNode interpolNode, PtrResultManager resMana);

  BinOpFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
    :BaseFilter(numWorkers,config,resMan){

  }


  virtual ~BinOpFilter(){

  }

  virtual bool Run() = 0;

protected:

  struct  PlusOp{
    static inline Double Apply(Double a1,Double a2){
      return a1+a2;
    }
  };

  struct  MultOp{
    static inline Double Apply(Double a1,Double a2){
      return a1*a2;
    }
  };

  struct  DivOp{
    static inline Double Apply(Double a1,Double a2){
      return a1/a2;
    }
  };

  struct  MinusOp{
    static inline Double Apply(Double a1,Double a2){
      return a1-a2;
    }
  };


  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

  std::string res1Name;
  uuids::uuid res1Id;

  std::string res2Name;
  uuids::uuid res2Id;

  std::string outName;
  uuids::uuid outId;

};

template<class Operator>
class GenericBinOpFilter : public BinOpFilter{

public:

  GenericBinOpFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~GenericBinOpFilter(){

  }

  virtual bool Run();

protected:

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

};

}

#endif /* SOURCE_CFSDAT_FILTERS_ARITHMETIC_BINOPFILTER_HH_ */
