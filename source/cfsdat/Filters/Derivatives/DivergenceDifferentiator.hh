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


namespace CFSDat{

class DivergenceDifferentiator : public MeshFilter{

  //! struct containing an interpolation matrix, which may be applied to scalars and vector
  struct Matrix {
    CF::UInt numTargets;
    StdVector< StdVector<CF::UInt> > targetSourceIndex;
    StdVector<CF::UInt> targetSource;
    StdVector< CF::Matrix<CF::Double> > targetSourceFactor;
  };


public:

  DivergenceDifferentiator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~DivergenceDifferentiator();


protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();


private:

  Grid* inGrid_;

  //! Entity map used for source values
  shared_ptr<EqnMapSimple> scrMap_;

  //! Entity map used for target values
  shared_ptr<EqnMapSimple> trgMap_;

  //! number of euqations per entity
  UInt numEquPerEnt_;

  std::vector<QuantityStruct> derivData_;

  //! Scaling of epsilon-parameter for RBF-basis function
  Double epsScal_;

  //! Scaling of beta-parameter for RBF-basis function
  Double betaScal_;

  //! Scaling of k-parameter for RBF-basis function
  Double kScal_;

  //! if true, a console output of [minimal distance, maximal distance, optimized epsilon]
  //! will be performed
  bool logEps_;

  //! index in the static matrices vector to use
  UInt matrixIndex_;

  //! contains pointers to every interpolator which created a matrix
  static CF::StdVector<DivergenceDifferentiator*> differentiators_;

  //! contains the matrices created by the Interpolators from interpolators_
  static CF::StdVector<Matrix> matrices_;


};

}
