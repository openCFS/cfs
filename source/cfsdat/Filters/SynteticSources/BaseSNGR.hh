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
//  BinOpFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
//    :BaseFilter(numWorkers,config,resMan){
//
//  }


  virtual ~SNGRFilter(){

  }

  virtual bool UpdateResults(std::set<uuids::uuid>& upResults);

protected:

  virtual ResultIdList SetUpstreamResults();

  virtual void AdaptFilterResults();

  //! for the mesh-check this mesh also needs to be stored, trgGrid_ is
  //! stored in MeshFilter
  Grid* inGrid_;


  std::string inTKE;
  uuids::uuid tkeId;

  std::string inTEF;
  uuids::uuid tefId;

  std::string inVelocity;
  uuids::uuid velocityId;

  std::string inDensity;
  uuids::uuid densityId;

  std::string inTemp;
  uuids::uuid temperatureId;

  std::string outName;
  uuids::uuid outId;

  std::string incrModes;
  uuids::uuid iModeId;

  Double TKEcrit;
  Double sigLength;
  Double fL;
  Double ft;
  Double fa;
  Double maxWN;
  Double minWN;
  UInt maxFreq;
  UInt minFreq;
  UInt numModes;
  UInt ensemble;
  
  Vector<Double> turbReconstVelocity;

  
  // time step of output signal
  Double deltaT;
  const Double pi = 3.14159265358979323846;
  

};

}

#endif /* SOURCE_CFSDAT_FILTERS_SNGR_HH_ */
