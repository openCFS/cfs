// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseDerivativeFilter.cc
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================


#include <Filters/Derivatives/CurlDifferentiator.hh>
#include "BaseDerivativeFilter.hh"
#include "Filters/Derivatives/GradientDifferentiator.hh"
#include "Filters/Derivatives/DivergenceDifferentiator.hh"
#include "Filters/Derivatives/Lighthill.hh"

namespace CFSDat{


FilterPtr BaseDerivativeFilter::GenerateSpatialDerivative(PtrParamNode derivNode, PtrResultManager resMana){
  FilterPtr newFilter;
 if(derivNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Gradient"){
   newFilter = FilterPtr(new CFSDat::GradientDifferentiator(0,derivNode,resMana));
 }else if(derivNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Divergence"){
   newFilter = FilterPtr(new CFSDat::DivergenceDifferentiator(0,derivNode,resMana));
 }else if(derivNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Curl"){
   newFilter = FilterPtr(new CFSDat::CurlDifferentiator(0,derivNode,resMana));
 }else if(derivNode->Get("type")->As<std::string>() == "AeroacousticSource_Lighthill" || "AeroacousticSource_LambVector"){
   newFilter = FilterPtr(new CFSDat::Lighthill(0,derivNode,resMana));
 }
 return newFilter;
}

}
