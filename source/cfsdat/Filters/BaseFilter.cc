// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseFilter.cc
 *       \brief    <Description>
 *
 *       \date     Aug 28, 2015
 *       \author   ahueppe
 */
//================================================================================================



#include "BaseFilter.hh"
#include "Filters/Arithmetic/BinOpFilter.hh"
#include "Filters/Arithmetic/VecOpFilter.hh"
#include "Filters/Input/InputFilter.hh"
#include "Filters/Output/OutputFilter.hh"
#include "Filters/Interpolators/BaseInterpolationFilter.hh"
#include "Filters/Derivatives/BaseDerivativeFilter.hh"
#include "Filters/Derivatives/RotatingSubstDt.hh"
#include "Filters/Derivatives/TimeDerivFilter.hh"
#include <boost/tokenizer.hpp>

namespace CFSDat{

FilterPtr BaseFilter::Generate(PtrParamNode filtNode,PtrResultManager resMana){
  FilterPtr newPtr;
if(filtNode->GetName() == "meshInput"){

  newPtr = FilterPtr(new CFSDat::InputFilter(0,filtNode,resMana));
}else if(filtNode->GetName() == "meshOutput"){
  newPtr = FilterPtr(new CFSDat::OutputFilter(0,filtNode,resMana));
}else if(filtNode->GetName() == "timeDeriv1"){
  newPtr = FilterPtr(new CFSDat::TimeDerivFilterD1(0,filtNode,resMana));
}else if(filtNode->GetName() == "substantialDeriv"){
  newPtr = FilterPtr(new CFSDat::RotatingSubstDt(0,filtNode,resMana));
}else if(filtNode->GetName() == "interpolation"){
  newPtr = BaseInterpolationFilter::GenerateInterpolator(filtNode,resMana);
}else if(filtNode->GetName() == "binaryOperation"){
  newPtr = BinOpFilter::GenerateOperator(filtNode,resMana);
}else if(filtNode->GetName() == "differentiation"){
  newPtr = BaseDerivativeFilter::GenerateSpatialDerivative(filtNode,resMana);
}else if (filtNode->GetName() == "vectorOperation"){
  newPtr = VecOpFilter::GenerateVectorOperator(filtNode,resMana);
}
return newPtr;
}


BaseFilter::BaseFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
            : numWorkers_(numWorkers),
              params_(config){

  std::cout << "\t---> Creating Filter with ID \"" << config->Get("id", CF::ParamNode::EX)->As<std::string>() << "\"" << std::endl;

  filtSteamType_ = BaseFilter::NO_STREAM;
  filterTag_ = boost::uuids::random_generator()();
  filterId_ = config->Get("id", CF::ParamNode::EX)->As<std::string>();

  std::string inIds = config->Get("inputFilterIds", CF::ParamNode::EX)->As<std::string>();

  typedef  boost::char_separator<char>  c_sep;
  typedef  boost::tokenizer< c_sep > tok;

  c_sep sep(" ,");
  tok varTok(inIds, sep);

  for(tok::iterator beg=varTok.begin(); beg!=varTok.end();++beg){
    inputIds_.insert(beg->c_str());
  }

  resultManager_ = resMan;

  //determine global step value map
  PtrParamNode stepNode = params_->GetParent()->Get("stepValueDefinition");
  if(stepNode->Has("startStop")){
    PtrParamNode stNode = stepNode->Get("startStop");
    UInt start    = stNode->Get("startStep")->Get("value")->As<UInt>();
    UInt numSteps = stNode->Get("numSteps")->Get("value")->As<UInt>();
    Double delta    = stNode->Get("delta")->Get("value")->As<Double>();
    //first step is always 1...
    for(UInt i=0;i<numSteps;++i){
      globalStepValueMap_[start+i] = (start+i)*delta;
    }
  }
}

void BaseFilter::InitResults(){
  std::cout << "\t---> Initializing Filter id " << this->filterId_ << std::endl;
  //search for results provided by filter itself
  //deactivate everything which should not go up
  ExtractFilterResults();

  //obtain results of upstream data
  //and give the implementing filter the possiblity to modify
  upResIds = this->SetUpstreamResults();

  //activate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->ActivateResult(upResIds[aRes]);
  }

  //loop over sources
  for(UInt aSrc=0;aSrc< sources_.GetSize(); aSrc++){
    sources_[aSrc]->InitResults();
  }

  //adapt the filter results according to the upstream data
  AdaptFilterResults();

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  for(UInt aRes=0;aRes<filterResIds.GetSize();aRes++){
    resultManager_->ActivateResult(filterResIds[aRes]);
  }
}

}
