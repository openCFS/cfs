// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DivergenceDifferentiator.hh
 *       \brief    <Description>
 *
 *       \date     Oct 11, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include <Filters/MeshFilter.hh>
#include "DataInOut/SimInput.hh"
#include <boost/tr1/type_traits.hpp>
#include <def_use_cgal.hh>
#include <def_use_flann.hh>
#include <cfsdat/Utils/Point.hh>

namespace CFSDat{



class DivergenceDifferentiator : public MeshFilter{

  struct DifferentiationStruct{
    CF::Vector<Double> localCoords;
    UInt tENum;
    UInt srcEqn;
    Double volume;

    DifferentiationStruct() : tENum(0),srcEqn(0),volume(.0){
      localCoords.Resize(3);
    }

    bool operator < (const DifferentiationStruct& str) const
    {
        return (srcEqn < str.srcEqn);
    }
  };

public:

  DivergenceDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~DivergenceDifferentiator();

  virtual bool Run();



protected:

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  // build the local RBF interpolation matrix, based on the nearest neighbour source points
  // build the local interpolation-value vector and solve the system ALoc*c=vector for
  // the local RBF coefficients c
  void CalcLocRBFDerivativeCoefs(CF::Matrix<Double>& vec,
                                 const CF::Vector<Double>& globPoint,
                                 const CF::StdVector< Vector<Double> >& neighbors,
                                 const CF::StdVector< Double >& l2Distances,
                                 const CF::StdVector< Vector<Double> >& vectors,
                                 const UInt numNN,
                                 const Double alpha);

private:

  std::vector<DifferentiationStruct> derivData_;

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;


};

}
