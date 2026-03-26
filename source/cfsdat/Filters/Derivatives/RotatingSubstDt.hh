// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     RotatingSubstDt.hh
 *       \brief    <Description>
 *
 *       \date     Nov 1, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_DERIVATIVES_ROTATINGSUBSTDT_HH_
#define SOURCE_CFSDAT_FILTERS_DERIVATIVES_ROTATINGSUBSTDT_HH_

#include "cfsdat/Filters/BaseFilter.hh"

namespace CFSDat{



//!This filter is basically a shortcut in computing the convective time derivative
//! with rotating domains. Right now everything is implemented here,
//! TODO: andi try to make it more ++, less C
class RotatingSubstDt : public BaseFilter{



public:
  RotatingSubstDt(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
    : BaseFilter(numWorkers,config, resMan){

    this->filtStreamType_ = FIFO_FILTER;

    if(params_->Has("rotatingDomain"))
      rpm_ = params_->Get("rotatingDomain")->Get("rpm")->As<Double>();

    // outputresultname
    std::string outResult = params_->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(outResult);
    
    // inputresultnames
    timeName_ = params_->Get("pressure")->Get("resultName")->As<std::string>();
    upResNames.insert(timeName_);
    
    gradName_ = params_->Get("presGrad")->Get("resultName")->As<std::string>();
    upResNames.insert(gradName_);
    
    hasMeanFlow_ = params_->Has("meanFlow");
    if(hasMeanFlow_){
      // meanflow is third
      meanFlowName_ = params_->Get("meanFlow")->Get("resultName")->As<std::string>();
      upResNames.insert(meanFlowName_);
    }
  }

  virtual ~RotatingSubstDt(){

  }

  virtual void PrepareCalculation();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:

  template<unsigned int D>
  void ExtractCylinderVelocities(PtrParamNode cylNode);

  Double rpm_;

  StdVector<UInt> rotEnts_;

  StdVector< Vector<Double> > rotField_;

  std::string timeName_;
  uuids::uuid timeId_;

  std::string gradName_;
  uuids::uuid gradId_;

  std::string meanFlowName_;
  uuids::uuid meanFlowId_;

  UInt gradDim_;

  Double dt_;

  bool hasMeanFlow_;
  
  bool hasRotation_;

};



}

#endif /* SOURCE_CFSDAT_FILTERS_DERIVATIVES_ROTATINGSUBSTDT_HH_ */
