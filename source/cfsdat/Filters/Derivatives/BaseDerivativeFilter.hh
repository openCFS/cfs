// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseDerivativeFilger.hh
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "cfsdat/Filters/BaseFilter.hh"


namespace CFSDat{

class BaseDerivativeFilter : public BaseFilter {
public:
  static FilterPtr GenerateSpatialDerivative(PtrParamNode derivNode, PtrResultManager resMana);

protected:
  ///Basic associative structure for mapping weighted values between source and target results
  struct AssocStruct{
    //Needs to be as lightweight as possible as it will be created millions of times
    //designed to perform differentiation in the following steps:

    ///weight value
    Double weight;

    ///target entity number
    UInt tENum;

    ///source entity number
    UInt sENum;

    ///Compare two instances according to the source entity number
    bool operator < (const AssocStruct& str) const
    {
        return (sENum < str.sENum);
    }
  };

public:

  BaseDerivativeFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
    :BaseFilter(numWorkers,config,resMan){

  }


  virtual ~BaseDerivativeFilter(){

  }

  virtual bool Run() = 0;

protected:

  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

  CF::StdVector<AssocStruct> derivativeData;


};

}
