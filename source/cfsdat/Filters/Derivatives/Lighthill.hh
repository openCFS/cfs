// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Lighthill.hh
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#pragma once


#include <Filters/Derivatives/AeroacousticBase.hh>
#include "DataInOut/SimInput.hh"


namespace CFSDat{

//! Class which provides methods to compute the Lamb- and LighthillSource-vector as well as
//! the whole LighthillSource-term
class Lighthill : public AeroacousticBase{

  //! struct containing an interpolation matrix, which may be applied to scalars and vector
  struct Matrix {
    CF::UInt numTargets;
    CF::UInt numSources;
    StdVector< StdVector<CF::UInt> > targetSourceIndexNtE;
    StdVector< StdVector<CF::UInt> > targetSourceIndexDiv;
    StdVector< StdVector<CF::UInt> > targetSourceIndexGrad;
    StdVector< StdVector<CF::UInt> > targetSourceIndexCurl;
    StdVector< StdVector<CF::UInt> > targetSourceIndexEtN; // EtN...element to node
    StdVector<CF::UInt> targetSourceNtE;
    StdVector<CF::UInt> targetSourceEtN;
    StdVector< CF::Matrix<CF::Double> > targetSourceFactorDiv;
    StdVector< CF::Matrix<CF::Double> > targetSourceFactorGrad;
    StdVector< CF::Matrix<CF::Double> > targetSourceFactorCurl;
    StdVector< CF::Matrix<CF::Double> > targetSourceNNFactor;
  };


public:

  Lighthill(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~Lighthill();


protected:

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
  
  virtual void PrepareCalculation();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  std::string resVelocityName;
  uuids::uuid resVelocityId;

  std::string resVorticityName;
  uuids::uuid resVorticityId;

  std::string resDensityName;
  uuids::uuid resDensityId;

private:

#ifdef USE_CGAL
  void LambVector(Vector<Double>& tempRetVec);

  void LighthillTensor(Vector<Double>& tempRetVec);

  void LighthillSourceVector(Vector<Double>& tempRetVec);

  void LighthillSourceTerm(Vector<Double>& tempRetVec, bool isTensorForm);

  Grid* inGrid_;

  //! number of euqations per entity
  UInt numEquPerEnt_;

  //! index in the static matrices vector to use
  UInt matrixIndex_;

  Integer stepIndex_;

#endif



  //! Entity map used for source values
  shared_ptr<EqnMapSimple> scrMap_;

  //! Entity map used for target values
  shared_ptr<EqnMapSimple> trgMap_;

  //! Number of neighbor points to include in differentiation.
  UInt numNeighbors_;

  std::vector<QuantityStruct> derivData_;

  //! Scaling of epsilon-parameter for RBF-basis function
  Double epsScal_;

  //! Scaling of beta-parameter of linear term of RBF-basis function
  Double betaScal_;

  //! Scaling of k-parameter of constant term of RBF-basis function
  Double kScal_;

  //! if true, a console output of [minimal distance, maximal distance, optimized epsilon]
  //! will be performed
  bool logEps_;


  //! contains pointers to every interpolator which created a matrix
  static CF::StdVector<Lighthill*> differentiators_;

  //! contains the matrices created by the Interpolators from interpolators_
  static CF::StdVector<Matrix> matrices_;

  //! Density, if not specified in xml-scheme it is automatically set to one
  Double density_ = 1.0;

  //! String if the full Lighthill or only the Lamb-vector is computed
  std::string Form_;

  //! Boolean if an extern vorticity-input is provided or if we have to compute it
  bool externVorticity_;

  //! Boolean if an extern density-input is provided when we have a variable density field (e. compressible flow)
  bool externDensity_;

  bool checkSum_;

};



}
