// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     sngrBase.hh
 *       \brief    <Description>
 *
 *       \date     Jul 24, 2017
 *       \author   r.krusche
 */
//================================================================================================

#ifndef SOURCE_CFSDAT_FILTERS_SNGR_HH_
#define SOURCE_CFSDAT_FILTERS_SNGR_HH_

#include "cfsdat/Filters/BaseFilter.hh"
#include <boost/bimap.hpp>
#include <algorithm>
#include <math.h>
#include <stdlib.h>

namespace CFSDat{

class SNGRFilter : public BaseFilter{
public:

  SNGRFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~SNGRFilter(){

  }

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);
//  virtual bool Run();

protected:

  virtual void FinishInit();

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

// utility functions
  virtual void GetTkeThreshold();
  virtual void SetRandVectors(UInt k, UInt j);
  virtual void SetVelocityAmplitude(UInt k, UInt j);

// method blocks
  virtual void PrepareMethodBailly();

  virtual void UpdateResultMethodBailly();

  virtual void SetTimelineMethodBailly();

  //virtual void PrepareMethodBillsonSNGR(); //BillsonSNGR, weil es gibt noch BillsonSAT

//  virtual void InitTimelineBillsonSNGR(UInt i);
//virtual void FinalTimelineMethodBillsonSNGR();

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter

// variables to initialise
  Grid* inGrid_;


  std::string inTKE_;
  uuids::uuid tkeId_;

  std::string inTEF_;
  uuids::uuid tefId_;

  std::string inVelocity_;
  uuids::uuid velocityId_;

  std::string inDensity_;
  uuids::uuid densityId_;

  std::string inTemp_;
  uuids::uuid temperatureId_;

  std::string outName_;
  uuids::uuid outId_;

  std::string interName_;
  uuids::uuid interId_;

  std::string incrModes_;
  uuids::uuid iModeId_;

  std::string method_;
  uuids::uuid methodId_;


  Double TKEcrit_;
  Double minTKE_;
  Double sigLength_;
  Double fL_;
  Double ft_;
  Double fa_;
  Double maxWN_;
  Double minWN_;
  Double c_mu_ = 0.09; // constand of k-epsilon model
  Double A_ = 1.452762;   // constand of SNGR
  Double deltaT_;
  Double deltaWN_; // wave number step in energy spectrum TODO see line 111 dWN_
  Double rmsOfVelFluct_;
  Double kolmogorovWN_;
  Double peakWN_;
  UInt maxFreq_;
  UInt minFreq_;
  UInt numModes_;
  UInt ensemble_;
  UInt numSteps_;
  UInt numNodes_;
  UInt tkeFAIL_ = 0;     // counter for nodes for which the deviation of reconstructed and read TKE is not small enough.
  UInt perpFAIL_ = 0;    // counter for nodes for which waveVec and dirVec are not sufficiently perpendicular to one another
//  UInt flg_=0; // counter for number of nodes the TKE-Criterion is met
  Vector<Double> idsNodesToProcess_;
  Vector<Double> waveVec_; // wave vector
  Vector<Double> dirVec_;  // direction vector of n-th mode
  Vector<Double> dWN_; // TODO, wenn incrModes=="logarithmic" ändert sich die step size mit jedem Mode, array benötigt, sonst reicht deltaWN_
//  Vector<Double> turbLengthScale_;
  Vector<Double> velAmplitude_;
  Vector<Double> omega_;
  Vector<Double> reconstTKE_;
  Vector<Double> turbReconstVelocity_;
  Vector<Double> stepValues_; // time step values of calculated time line
  Vector<Double> phase_;
  Vector<Double> kn_;  // array of all wave number values
  Vector<Double> randAngles_; // array of theta and phi for defining wave vector in space
  //INPUT from RANS solution
  Vector<Double> TKE_;
  Vector<Double> meanVelocity_;
  Vector<Double> localDensity_;
  Vector<Double> localTemp_;
  Vector<Double> TEF_;

  // constants and parameters SNGR
  //Double fL = 2.5;
  //Double ft_ = 0.0002;

};

}

#endif /* SOURCE_CFSDAT_FILTERS_SNGR_HH_ */
