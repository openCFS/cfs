// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CentroidInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Nov 14, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_
#define SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_


#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"
#include "MatVec/CRS_Matrix.hh"


namespace CFSDat{

//! Class for choosing a very simple but fast interpolation approach
//! Upon initialization we determine cell centroids and volume from input
//! additionally we set the local coordinates accoring to source
//! during traversal, we just apply those loads
class CentroidInterpolator : public MeshFilter{

public:

  CentroidInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~CentroidInterpolator();

protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

private:
//  void CreateCRS(const Grid* inGrid, const Grid* trgGrid);

//  void FillInterpolationMatrix(const StdVector<ElemIntersect::VolCenterInfo> & infos);

  CRS_Matrix<Double>* InterpolationMatrix;


//  std::vector<QuantityStruct> interpolData_;

  //! Global Factor for scaling the result
  Double globalFactor_;

  bool checkSum_;
};

}
#endif /* SOURCE_CFSDAT_FILTERS_INTERPOLATORS_CENTROIDINTERPOLATOR_HH_ */
