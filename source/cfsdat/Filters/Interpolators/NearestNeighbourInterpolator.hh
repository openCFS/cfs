// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearestNeighbourInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once //instead of the #ifndef #def ...


#include "DataInOut/SimInput.hh"
#include <Filters/MeshFilter.hh>
#include "Utils/InterpolationMatrix.hh"
#include "AbstractInterpolator.hh"



namespace CFSDat{

//! Class for calculating a nearest neighbour interpolation using CGAL
//! for neighbour search

class NearestNeighbourInterpolator : public AbstractInterpolator {
  

public:

  NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~NearestNeighbourInterpolator();

protected:
  
  virtual void FillMatrix(StdVector<CF::UInt>& globSrcEntity, StdVector<CF::UInt>& globTrgEntity);

private:


  //! Number of neighbor points to include in interpolation.
  UInt numNeighbors_;
  //! Maximum distance between source and target positions to include into the interpolation operation
  Double maxSearchDist_;

#ifdef USE_CGAL // clang otherwise complains about private field not used
  //! Global Factor for scaling the result
  Double globalFactor_;
#endif

  //! Exponent for calculation of interpolation weight function.
  UInt p_;
  
};

}
