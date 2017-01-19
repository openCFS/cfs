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

  //! Struct stores the LOCAL coordinates of all elements in the defined regions
  //! and the element nr. for the specific elements
  struct InterpolationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    InterpolationStruct() : tENum(0){
      localCoords.Resize(3);
    }
  };


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

  //! Method, containing the actual interpolation algorithm

  //! \param returnVec (out) vector containing the interpolated results
  //! \param vec (in) vector containing the interpolated result for every
  //!                 target-node, can be 1D for scalar values or 2D/3D
  //!                 for vector values, no return-parameter, only needed
  //!                 for intern-calculation
  //! \param downMap (in) pointer to equation numbers of target-mesh
  //! \param srcCoords (in) coordinate matrix of source-mesh
  //! \param trgCoords (in) coordinate matrix of target-mesh
  //! \param interpolationData (in) information of elements of target-regions
  //! \param grid (in) pointer to the target-grid
  void Interpolation(Vector<Double>& returnVec,
                     const CF::StdVector< CF::Vector<Double> >&  scatteredData,
                     CF::Vector<Double>& vec,
                     const str1::shared_ptr<EqnMapSimple>& downMap,
                     const CF::StdVector< CF::Vector<double> >& srcCoords,
                     const CF::StdVector< CF::Vector<double> >& trgCoords,
                     std::vector<InterpolationStruct>& interpolationData,
                     Grid* grid);

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


  //! Here we bring the input-results and the input-coordinate into the
  //! appropriate form a nearest neighbour search with CGAL.
  //! The problem is that the input-results, provided by result-manager, is
  //! only a vector containing the results in the form (u1,v1,w1,u2,v2,w2,...)
  //! and for CGAL we need something like ([u1,v1,w1],[u2,v2,w2],...)

  //! \param scatteredData (out) Vector containing the correct ordered input-data
  //! \param vec (in) vector containing the interpolated result for every
  //!                 target-node, can be 1D for scalar values or 2D/3D
  //!                 for vector values, no return-parameter, only needed
  //!                 for intern-calculation
  //! \param inVec (in) vector containing the input-results, provided by the
  //!                   result-manager
  //! \coords (in) locations, where the values of inVec are defined
  void FillScatteredDataVec(CF::StdVector< CF::Vector<Double> >& scatteredData,
                            CF::Vector<Double>& vec,
                            const Vector<Double>& inVec,
                            const CF::StdVector< CF::Vector<double> >& coords);


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
  std::vector<InterpolationStruct> interpolDataTrg_;

  //! vector for the mesh-check, for interpolation from Trg->Src
  //! reversed forward interpolation
  std::vector<InterpolationStruct> interpolDataSrc_;

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
  Grid* inGrid_;

};

}
