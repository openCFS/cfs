// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearesNeighbourInterpolator.hh
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



namespace CFSDat{

//! Class for calculating a nearest neighbour interpolation using CGAL
//! for neighbour search

class NearestNeighbourInterpolator : public MeshFilter{

public:

  NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~NearestNeighbourInterpolator();

  virtual bool Run();

protected:

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  //! Read scattered data
  void ReadScatteredData(CF::StdVector< CF::Vector<Double> > elemCentroids, CF::StdVector< CF::Vector<Double> > scatteredData);


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


  //! Unfortunately we need the coordinates in such a "hardcoded" way, because of the CGAL-nearest neighbour search
  //! Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  //! Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;

  //! Number of neighbor points to include in interpolation.
  UInt numNeighbors_;

  //! Exponent for calculation of interpolation weight function.
  Double p_;

  //! Boolean, which stores the binary input from xml-file, if we
  //! perform a mesh-check. If mesh-check is set in xml-> mCheck = true
  bool mCheck_;

  //! primary vector of struct for interpolation from Src->Trg
  //! forward interpolation
  std::vector<QuantityStruct> interpolDataTrg_;

  //! vector for the mesh-check, for interpolation from Trg->Src
  //! reversed forward interpolation
  std::vector<QuantityStruct> interpolDataSrc_;

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
  Grid* inGrid_;

};

}
