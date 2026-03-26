// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     VecOpFilter.hh
 *       \brief    <Description>
 *
 *       \date     Oct 19, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>

namespace CFSDat{

/* This serves as a base class for vector operations (Vector x Vector, Vector * Vector)
 * It can be extended to other operation, e.g. tensor operations
 */



class VecOpFilter : public BaseFilter{
public:

  static FilterPtr GenerateVectorOperator(PtrParamNode node, PtrResultManager resMana);

  VecOpFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
    :BaseFilter(numWorkers,config,resMan){

  }


  virtual ~VecOpFilter(){

  }

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
class GenericVecOpFilter : public VecOpFilter{

public:

  GenericVecOpFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~GenericVecOpFilter(){

  }

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

};

}
