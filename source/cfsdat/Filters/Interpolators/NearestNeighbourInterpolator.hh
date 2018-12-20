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

  //! Global Factor for scaling the result
  Double globalFactor_;

  //! Exponent for calculation of interpolation weight function.
  UInt p_;

  //! number of euqations per entity
  UInt numEquPerEnt_;
  
  //! Entity map used for source values
  str1::shared_ptr<EqnMapSimple> scrMap_;

  //! Entity map used for target values
  str1::shared_ptr<EqnMapSimple> trgMap_;
  
  //! Boolean, which stores the binary input from xml-file, if we
  //! perform a mesh-check. If mesh-check is set in xml-> mCheck = true
  bool mCheck_;

  //! if true, then element values are the interpolation target
  bool useElemAsTarget_;

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
  Grid* inGrid_;
  
  //! index in the static matrices vector to use
  UInt matrixIndex_;
  
  //! contains pointers to every interpolator which created a matrix
  static CF::StdVector<NearestNeighbourInterpolator*> interpolators_;

  //! contains the matrices createb by the Interpolators from interpolators_
  static CF::StdVector<Matrix> matrices_;

};

}
