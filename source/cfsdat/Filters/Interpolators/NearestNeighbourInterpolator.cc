// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearestNeighbourInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================



#include "AbstractInterpolator.hh"
#include "NearestNeighbourInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "cfsdat/Utils/KNNIndexSearch.hh"
//#include "cfsdat/Utils/KNNSearch.hh"
#include <algorithm>
#include <vector>

#include "DataInOut/Logging/LogConfigurator.hh"

// check if Intel MKL is available
#include <def_use_blas.hh>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef USE_MKL
  #include <mkl.h>
#endif

namespace CFSDat{

  // declare logging stream
  DEFINE_LOG(nearestNeighbourInterpolator, "NearestNeighbourInterpolator")

// constructor
NearestNeighbourInterpolator::NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan):AbstractInterpolator(numWorkers,config,resMan) {
  interpolatorName_ = "NearestNeighbourInterpolator";
#ifndef USE_CGAL
    EXCEPTION(interpolatorName_ << " needs to be compiled with USE_CGAL=ON!");
#endif
  // get xml info
  p_ = config->Get("IntSchemeNN")->Get("interpolationExponent")->As<UInt>();
  numNeighbors_ = config->Get("IntSchemeNN")->Get("numNeighbours")->As<UInt>();
  globalFactor_ = config->Get("IntSchemeNN")->Get("globalFactor")->As<Double>();
  maxSearchDist_ = config->Get("IntSchemeNN")->Get("maxDistance")->As<Double>();
  if ((maxSearchDist_ < 0.0) && (maxSearchDist_ != -1.0))
    WARN("Ignoring negative maxDistance for NearestNeighbourInterpolator.");
  useElemAsTarget_ = false;
  if(params_->Has("useElemAsTarget"))
    useElemAsTarget_ = params_->Get("useElemAsTarget")->As<bool>();
}

// destructor
NearestNeighbourInterpolator::~NearestNeighbourInterpolator(){}

void NearestNeighbourInterpolator::FillMatrix(StdVector<CF::UInt>& globSrcEntity, StdVector<CF::UInt>& globTrgEntity) {
  UInt maxNumSrcEntities = globSrcEntity.GetSize();
  UInt maxNumTrgEntities = globTrgEntity.GetSize();
  numNeighbors_ = std::min(numNeighbors_, maxNumSrcEntities);
  
  std::cout << "\t\t Creating nearest neighbour search tree " << std::endl;
  KNNIndexSearch indexSearch;
  // get the coordinates
  for(CF::UInt srcEnt = 0; srcEnt < maxNumSrcEntities; srcEnt++) {
    CF::UInt globEntityNumber = globSrcEntity[srcEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      if (useElemAsSource_) {
        inGrid_->GetElemCentroid(pCoord,globEntityNumber,true);
      } else {
        inGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      }
      indexSearch.AddPoint(pCoord,srcEnt);
    }
  }
  indexSearch.BuildTree();
  
  std::cout << "\t\t Creating interpolation matrix " << std::endl;
  matrix_.inCount = maxNumSrcEntities;
  matrix_.outCount = maxNumTrgEntities;
  
  StdVector<CF::UInt>& outInIndex = matrix_.outInIndex;
  StdVector<CF::UInt>& outIn = matrix_.outIn;
  StdVector<CF::Double>& outInFactor = matrix_.outInFactor;
  
  // creating target -> source indices
  outInIndex.Resize(maxNumTrgEntities + 1);
  CF::UInt index = 0;
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    outInIndex[trgEnt] = index;
    if (globTrgEntity[trgEnt] != UnusedEntityNumber) {
      index += numNeighbors_;
    }
  }
  outInIndex[maxNumTrgEntities] = index;
  
  // creating target -> source factors and filling source indices
  outIn.Resize(index);
  outInFactor.Resize(index);
  #pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      if(useElemAsTarget_){
        trgGrid_->GetElemCentroid(pCoord, globEntityNumber,true);
      } else {
        trgGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      }
      StdVector<UInt> nearestIndices;
      StdVector<Double> distances;
      indexSearch.QueryPoint(pCoord, nearestIndices, distances, numNeighbors_);
      const CF::UInt ITrgSrcIndex = outInIndex[trgEnt];
      if (numNeighbors_ > 1) {
        CF::UInt* srcIndices = &outIn[ITrgSrcIndex];
        CF::Double* srcFactors = &outInFactor[ITrgSrcIndex];
        CF::UInt i;
        CF::Double dmax = 0.0; // set to zero so that all distances will be >= max at init
        CF::Double dmin = 1e10; // set to large number so that all distances will be < dmin at init
        for (i = 0; i < numNeighbors_; i++) {
          srcIndices[i] = nearestIndices[i];
          srcFactors[i] = distances[i];
          dmax = std::max(distances[i],dmax);
          dmin = std::min(distances[i],dmin);
        }
        // if we have defined a maximum search distance and there is no closest neighbor within this distance, set the solution to 0.0
        if ((maxSearchDist_ >= 0) && (dmin > maxSearchDist_)) {
          for (i = 0; i < numNeighbors_; ++i) {
            srcFactors[i] = 0.0;
          }
        } else {
          // Apply Shepard interpolation cf. Numerical Recipes 3rd ed. p. 143ff.
          // or http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes \
          // /Inverse_Distance_Weighted/Shepards_Method.htm
          const CF::Double R = 1.01 * dmax;
          CF::Double weights = 0.0;
          // The point which is farthest away, should at least have a non-zero
          // weight of 0.01. If we would choose R = dmax, it would not contribute
          // at all.
          for (i = 0; i < numNeighbors_; i++) {
            Double w;
            if(srcFactors[i] == 0){
              w = 1.0;
            } else {
              w = std::pow((R-srcFactors[i])/(R*srcFactors[i]), p_);
            }
            srcFactors[i] = w;
            weights += w;
          }
          for (i = 0; i < numNeighbors_; i++) {
            srcFactors[i] /= weights;
            LOG_TRACE(nearestNeighbourInterpolator) << "\t\t\t scrFactor " << srcFactors[i] << std::endl;

            srcFactors[i] = srcFactors[i] * globalFactor_;
            LOG_TRACE(nearestNeighbourInterpolator) << "\t\t\t scrFactor mult " << srcFactors[i] << std::endl;
          }
        }
      } else {
        // here, only one neighbor is included.
        outIn[ITrgSrcIndex] = nearestIndices[0];
        outInFactor[ITrgSrcIndex] = 1.0;
      }
    }
  }
}
}
