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
#include <cfsdat/Utils/Point.hh>
#include <Filters/MeshFilter.hh>

#include <limits>



namespace CFSDat{

//! Class for calculating a nearest neighbour interpolation using CGAL
//! for neighbour search

class NearestNeighbourInterpolator : public MeshFilter{
  
  //! Entity number for unused indices in in- and output results
  static const CF::UInt UnusedEntityNumber = UINT_MAX;
  
  //! struct containing an interpolation matrix, which may be applied to scalars and vector
  struct Matrix {
    CF::UInt numTargets;
    StdVector<CF::UInt> targetSourceIndex;
    StdVector<CF::UInt> targetSource;
    StdVector<CF::Double> targetSourceFactor;
  };
  
public:

  NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~NearestNeighbourInterpolator();

  virtual bool Run();

protected:
  
  CF::UInt CountUsedEntities(StdVector<CF::UInt>& entities);
  
  void GetUsedMappedEntities(str1::shared_ptr<EqnMapSimple>& map, StdVector<CF::UInt>& entities, std::set<std::string>& regions, Grid* grid);

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  bool CreatesEqualMatrix(NearestNeighbourInterpolator* otherInterpolator);

private:

  //! Method which performs comparison of original results and twice interpolated ones
  //! in order to check the difference between the initial source-values and the
  //! two-times-interpolated ones. With these two results we can compute the
  //! correlation between the two results and throw an exception if the interpolation
  //! is not good enough, so the target-mesh has to be refined or the interpolation-
  //! exponent has to be adapted

  //! \param (in) origVec original input vector of filter (not interpolated)
  //! \param (in) newVec twice interpolated vector Src->Trg->Src
  void CheckFilterResults(Vector<Double>& origVec,
                          Vector<Double>& newVec);


  //! Number of neighbor points to include in interpolation.
  UInt numNeighbors_;

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
