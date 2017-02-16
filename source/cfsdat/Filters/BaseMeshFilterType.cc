// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseMeshFilterType.cc
 *       \brief    <Description>
 *
 *       \date     Jan 11, 2017
 *       \author   kroppert
 */
//================================================================================================


#include <Filters/BaseMeshFilterType.hh>

// Interpolators
#include "Filters/Interpolators/CentroidInterpolator.hh"
#include "Filters/Interpolators/NearestNeighbourInterpolator.hh"
#include "Filters/Interpolators/Cell2NodeInterpolator.hh"
#include "Filters/Interpolators/FEBasedInterpolator.hh"
#include "Filters/Interpolators/GridIntersectionFilter.hh"
#include "Filters/Interpolators/RBFInterpolator.hh"
#include "Filters/Interpolators/Node2CellInterpolator.hh"

// Differentiators
#include "Filters/Derivatives/CurlDifferentiator.hh"
#include "Filters/Derivatives/GradientDifferentiator.hh"
#include "Filters/Derivatives/DivergenceDifferentiator.hh"

// Acoustic Source Terms
#include "Filters/Derivatives/Lighthill.hh"



namespace CFSDat{


FilterPtr BaseMeshFilterType::Generate(PtrParamNode ptrNode, PtrResultManager resMana){
  FilterPtr newFilter;
 if(ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_Conservative_CellCentroid"){
   newFilter = FilterPtr(new CFSDat::CentroidInterpolator(0,ptrNode,resMana));
 }
 else if(ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_Conservative_CutCell"){
   newFilter = FilterPtr(new CFSDat::GridIntersectionFilter(0,ptrNode,resMana));
 }
 else if (ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_NearestNeighbour"){
   newFilter = FilterPtr(new CFSDat::NearestNeighbourInterpolator(0,ptrNode,resMana));
 }
 else if (ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_Cell2Node"){
   newFilter = FilterPtr(new CFSDat::Cell2NodeInterpolator(0,ptrNode,resMana));
 }
 else if (ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_RBF"){
   newFilter = FilterPtr(new CFSDat::RBFInterpolator(0,ptrNode,resMana));
 }
 else if (ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_Node2Cell"){
   newFilter = FilterPtr(new CFSDat::Node2CellInterpolator(0,ptrNode,resMana));
 }
 else if(ptrNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Gradient"){
   newFilter = FilterPtr(new CFSDat::GradientDifferentiator(0,ptrNode,resMana));
 }
 else if(ptrNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Divergence"){
   newFilter = FilterPtr(new CFSDat::DivergenceDifferentiator(0,ptrNode,resMana));
 }
 else if(ptrNode->Get("type")->As<std::string>() == "SpaceDifferentiation_Curl"){
   newFilter = FilterPtr(new CFSDat::CurlDifferentiator(0,ptrNode,resMana));
 }
 else if(ptrNode->Get("type")->As<std::string>() == "AeroacousticSource_LighthillSourceTerm" ||
         ptrNode->Get("type")->As<std::string>() == "AeroacousticSource_LighthillSourceVector" ||
         ptrNode->Get("type")->As<std::string>() == "AeroacousticSource_LambVector"){
   newFilter = FilterPtr(new CFSDat::Lighthill(0,ptrNode,resMana));
 }
 else if (ptrNode->Get("type")->As<std::string>() == "FieldInterpolation_FEBased"){
    newFilter = FilterPtr(new CFSDat::FEBasedInterpolator(0,ptrNode,resMana));
  }
 return newFilter;
}


}
