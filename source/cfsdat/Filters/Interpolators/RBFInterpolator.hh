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
//#include <cfsdat/Utils/Point.hh>
#include <Filters/MeshFilter.hh>


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

  RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~RBFInterpolator();

  virtual bool Run();



protected:



  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual Double DistanceECLID(CF::Vector<Double> p1, CF::Vector<Double> p2);

  virtual void PreparePATCH();

  virtual void PrepareCGAL();

  virtual void AdaptFilterResults();

private:


    //! Number of neighbor points to include in interpolation.
    UInt numNN_;
    UInt numNW_;


    //! number of euqations per entity
    UInt numEquPerEnt_;

    //! Entity map used for source values
    str1::shared_ptr<EqnMapSimple> scrMap_;

    //! Entity map used for target values
    str1::shared_ptr<EqnMapSimple> trgMap_;

    Grid* inGrid_;

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

    //! index in the static matrices vector to use
    UInt matrixIndex_;

    //! contains pointers to every interpolator which created a matrix
    static CF::StdVector<RBFInterpolator*> interpolators_;

    //! contains the matrices createb by the Interpolators from interpolators_
    static CF::StdVector<Matrix> matrices_;

    std::string expr_;

    MathParser* mp_;
};

}
