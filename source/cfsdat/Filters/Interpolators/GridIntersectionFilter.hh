// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GridIntersectionFilter.hh
 *       \brief    Class for providing cut cell based interpolation
 *                 capabilities
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_INTERPOLATORS_GRIDINTERSECTIONFILTER_HH_
#define SOURCE_CFSDAT_FILTERS_INTERPOLATORS_GRIDINTERSECTIONFILTER_HH_

#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/CoordFormat.hh"
#include "Domain/Mesh/MeshUtils/Intersection/IntersectAlgos/ElemIntersect.hh"

namespace CFSDat{

//! This filter provides a conservative interpolation based on
//! volume Grid intersection

//! The basic algorithms is as follows
//!  - call the volume intersection class to intersect the given grids
//!  - obtain intersection volumes and centroids as a result
//!  - based on these pairs, create a sparse matrix
//!  - each pair can be converted to a weight and is thereby summed to the matrix entries
//!  - during runtime only a matrix vector multiplication is performed, all costly objects are cleared
class GridIntersectionFilter : public MeshFilter{

public:

  GridIntersectionFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~GridIntersectionFilter();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:
  void CreateCRS(const StdVector<ElemIntersect::VolCenterInfo> & infos);

  void FillInterpolationMatrix(const StdVector<ElemIntersect::VolCenterInfo> & infos);

  CRS_Matrix<Double>* InterpolationMatrix;

  //! Global Factor for scaling the result
  Double globalFactor_;

};

}

#endif /* SOURCE_CFSDAT_FILTERS_INTERPOLATORS_GRIDINTERSECTIONFILTER_HH_ */
