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
#include <cfsdat/Utils/Point.hh>


namespace CFSDat{

//! Class for calculating interpolation using radial basis functions (RBF)
class RBFInterpolator : public MeshBasedInterpolator{

public:
  struct InpolationStruct{
    CF::Vector<Double> localCoords; //target local coordinates
    UInt tENum; //target element number

    InpolationStruct() : tENum(0){
      localCoords.Resize(3);
    }
  };

  RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~RBFInterpolator();

  virtual bool Run();



protected:

  virtual void PrepareInterpolation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();


  // build the local RBF interpolation matrix, based on the nearest neighbour source points
  // build the local interpolation-value vector and solve the system ALoc*c=vector for
  // the local RBF coefficients c
  void CalcLocRBFCoefs(CF::Matrix<Double>& coefVec,
                       const CF::Vector<Double>& globPoint,
                       const CF::StdVector< Vector<Double> >& neighbors,
                       const CF::StdVector< Double >& l2Distances,
                       const CF::StdVector< Vector<Double> >& vectors,
                       const UInt numNN,
                       const Double alpha);

private:
    // Coordinates of input data
    CF::StdVector< CF::Vector<double> > sourceCoords_;

    // Coordinates of target data
    CF::StdVector< CF::Vector<double> > targetCoords_;

    //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
    UInt inDim_;

    //! Exponent for calculation of interpolation weight function.
    Double p_;

    //! Vector containing target-informations
    std::vector<InpolationStruct> interpolData_;

};

}
