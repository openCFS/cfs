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

  virtual void PrepareMethodBailly();

  virtual void UpdateResultMethodBailly();

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
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

  Double TKEcrit_;
  Double sigLength_;
  Double fL_;
  Double ft_;
  Double fa_;
  Double maxWN_;
  Double minWN_;
  UInt maxFreq_;
  UInt minFreq_;
  UInt numModes_;
  UInt ensemble_;
  
  Vector<Double> turbReconstVelocity_;

  // time step of output signal
  Double deltaT_;
  

};

}

#endif /* SOURCE_CFSDAT_FILTERS_SNGR_HH_ */
