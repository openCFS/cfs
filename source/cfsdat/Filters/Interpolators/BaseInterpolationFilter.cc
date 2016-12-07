// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CaseInterpolationFilter.cc
 *       \brief    <Description>
 *
 *       \date     Dec 2, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "BaseInterpolationFilter.hh"
#include "Filters/Interpolators/CentroidInterpolator.hh"
#include "Filters/Interpolators/NearestNeighbourInterpolator.hh"
#include "Filters/Interpolators/Cell2NodeInterpolator.hh"
#include "Filters/Interpolators/GridIntersectionFilter.hh"
#include "Filters/Interpolators/RBFInterpolator.hh"


namespace CFSDat{


FilterPtr BaseInterpolationFilter::GenerateInterpolator(PtrParamNode interpolNode, PtrResultManager resMana){
  FilterPtr newFilter;
 if(interpolNode->Get("type")->As<std::string>() == "FieldInterpolation_Conservative_CellCentroid"){
   newFilter = FilterPtr(new CFSDat::CentroidInterpolator(0,interpolNode,resMana));
 } else if(interpolNode->Get("type")->As<std::string>() == "FieldInterpolation_Conservative_CutCell"){
   newFilter = FilterPtr(new CFSDat::GridIntersectionFilter(0,interpolNode,resMana));
 } else if (interpolNode->Get("type")->As<std::string>() == "FieldInterpolation_NearestNeighbour"){
   newFilter = FilterPtr(new CFSDat::NearestNeighbourInterpolator(0,interpolNode,resMana));
 } else if (interpolNode->Get("type")->As<std::string>() == "FieldInterpolation_Cell2Node"){
   newFilter = FilterPtr(new CFSDat::Cell2NodeInterpolator(0,interpolNode,resMana));
 } else if (interpolNode->Get("type")->As<std::string>() == "FieldInterpolation_RBF"){
   newFilter = FilterPtr(new CFSDat::RBFInterpolator(0,interpolNode,resMana));
 }
 return newFilter;
}


}
