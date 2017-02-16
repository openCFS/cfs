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
#include <cfsdat/Utils/Point.hh>

namespace CFSDat{



class DivergenceDifferentiator : public MeshFilter{

public:

  DivergenceDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~DivergenceDifferentiator();

  virtual bool Run();



protected:

  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();


private:

  std::vector<QuantityStruct> derivData_;

  // Coordinates of input data
  CF::StdVector< CF::Vector<double> > sourceCoords_;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;

  //! Dimension of input values (0=scalar, 1=two-dim vector, 2=three-dim vector).
  UInt inDim_;


};

}
