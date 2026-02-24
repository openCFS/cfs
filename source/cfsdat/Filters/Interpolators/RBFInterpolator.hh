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


#include "DataInOut/SimInput.hh"
//#include <cfsdat/DatUtils/Point.hh>
#include <Filters/MeshFilter.hh>
#include <def_use_cgal.hh>


namespace CFSDat{

//! Class for calculating interpolation using radial basis functions (RBF)
class RBFInterpolator : public MeshFilter{


  //! struct containing an interpolation matrix, which may be applied to scalars and vector
  struct Matrix {
    CF::UInt numTargets;
    StdVector<CF::UInt> targetSourceIndex;
    StdVector<CF::UInt> targetSource;
    StdVector<CF::Double> targetSourceFactor;
    StdVector<CF::Double> targetSourceFactor2;
    StdVector< CF::Matrix<CF::Double> > targetRBFInvMat;
  };


public:

  RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~RBFInterpolator();

protected:


  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

#ifdef USE_CGAL
  Double DistanceEUCLID(CF::Vector<Double> p1, CF::Vector<Double> p2);

  void PreparePATCH();

  void PrepareCGAL();

  //! number of equations per entity
  UInt numEquPerEnt_ = 0;

  Grid* inGrid_ = NULL;

  //! index in the static matrices vector to use
  UInt matrixIndex_ = 0;

  MathParser* mp_ = NULL;
#endif

private:


    //! Number of neighbor points to include in interpolation.
    UInt numNN_;
    UInt numNW_;

    //! Global Factor for scaling the result
    Double globalFactor_;

    //! Entity map used for source values
    shared_ptr<EqnMapSimple> scrMap_;

    //! Entity map used for target values
    shared_ptr<EqnMapSimple> trgMap_;


    //! Coordinates of input data
    CF::StdVector< CF::Vector<double> > sourceCoords_;

    //! Coordinates of target data
    CF::StdVector< CF::Vector<double> > targetCoords_;

    //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
    UInt inDim_;

    //! Exponent for calculation of interpolation weight function.
    Double p_;

    //! Vector containing target-informations
    std::vector<QuantityStruct> interpolData_;

    //! if true, then all nodes on region "wall" will have 0.0 as return-vector entry
    bool noSlip_;
    //! boundary set
    StdVector<UInt> boundary_;

    //! if true, then element values are the interpolation target
    bool useElemAsTarget_;


    //! contains pointers to every interpolator which created a matrix
    static CF::StdVector<RBFInterpolator*> interpolators_;

    //! contains the matrices createb by the Interpolators from interpolators_
    static CF::StdVector<Matrix> matrices_;

    std::string expr_;
};

}
