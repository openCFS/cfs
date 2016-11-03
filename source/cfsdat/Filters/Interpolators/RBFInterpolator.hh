// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     RBFInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Sep 8, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include "MeshBasedInterpolator.hh"
#include "DataInOut/SimInput.hh"
#include <boost/tr1/type_traits.hpp>
#include <def_use_cgal.hh>
#include <def_use_flann.hh>
#include <cfsdat/Utils/Point.hh>


namespace CFSDat{

//! Class for calculating interpolation using radial basis functions (RBF)


class RBFInterpolator : public MeshBasedInterpolator{

  struct InpolationStruct{
    CF::Vector<Double> localCoords; //target local coordinates
    UInt tENum; //target element number


    InpolationStruct() : tENum(0){
      localCoords.Resize(3);
    }


  };

public:

  RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~RBFInterpolator();

  virtual bool Run();



protected:

  virtual void PrepareInterpolation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  // Read scattered data
  void ReadScatteredData(CF::StdVector< CF::Vector<Double> > elemCentroids,
                         CF::StdVector< CF::Vector<Double> > scatteredData);

  // build the local RBF interpolation matrix, based on the nearest neighbour source points
  // build the local interpolation-value vector and solve the system ALoc*c=vector for
  // the local RBF coefficients c
  void CalcLocRBFCoefs(CF::Matrix<Double>& coefVec,
                                        CF::Vector<Double>& globPoint,
                                        CF::StdVector< Vector<Double> >& neighbors,
                                        CF::StdVector< Double >& l2Distances,
                                        CF::StdVector< Vector<Double> >& vectors,
                                        UInt numNN,
                                        Double alpha);

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;

  //! Exponent for calculation of interpolation weight function.
  Double p_;

  //! Library used to find the k nearest neighbors of a point.
  //KNNLibary knnLib_;
  //0 for CGAL, 1 for FLANN
  UInt knnLib_;

#ifdef USE_CGAL
    boost::shared_ptr<Tree> searchTree_;

    void KNNSearch_CGAL(const Vector<Double> globPoint,
                        StdVector< Vector<Double> >& neighbors,
                        StdVector< Double >& l2Distances,
                        StdVector< Vector<Double> >& vectors,
                        UInt numNeighbors);
#endif

#ifdef USE_FLANN
    boost::shared_ptr< flann::Index<flann::L2<Double> > > index_;
    boost::shared_ptr< flann::Matrix<Double> > dataset_;

    void KNNSearch_FLANN(const Vector<Double> globPoint,
                         StdVector< Vector<Double> >& neighbors,
                         StdVector< Double >& l2Distances,
                         StdVector< Vector<Double> >& vectors,
                         StdVector< Vector<Double> > scatteredData,
                         UInt numNeighbors);
#endif


private:

  std::vector<InpolationStruct> interpolData_;
  StdVector<UInt> nodeNeighbours_;

};

}
