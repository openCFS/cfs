// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AbstractInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Apr 4, 2018
 *       \author   mtautz
 */
//================================================================================================



#include "NearestNeighbourInterpolator.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>


// check if Intel MKL is available
#include <def_use_blas.hh>
#ifdef USE_MKL
  #include <mkl.h>
#endif

namespace CFSDat{

AbstractInterpolator::AbstractInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  interpolatorName_ = "AbstractInterpolator";
  useElemAsSource_ = false;
  useElemAsTarget_ = false;
  if(params_->Has("useElemAsTarget")){useElemAsTarget_ = params_->Get("useElemAsTarget")->As<bool>();}
  verboseSum_ = false;

  // Global factor removed from xml scheme. Therefore, it is hard-coded. Use binary operation filter if required.
  globalFactor_= 1.0;

}

AbstractInterpolator::~AbstractInterpolator(){

}

bool AbstractInterpolator::UpdateResults(std::set<uuids::uuid>& upResults){
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);
  
  Double* out = &returnVec[0];
  Double* in = &inVec[0];
  matrix_.Interpolate(out, in, numEquPerEnt_);
  returnVec.ScalarMult(globalFactor_);
  
  if (verboseSum_) {
    VerboseSum(upResIds[0]);
    VerboseSum(filterResIds[0]);
  }
  return true;
}


void AbstractInterpolator::PrepareCalculation(){
  std::cout << "\t ---> " << interpolatorName_  << " preparing for interpolation" << std::endl;
  
  std::cout << "\t\t Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  scrMap_ = resultManager_->GetEqnMap(upRes);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  useElemAsSource_ = inInfo->definedOn == ExtendedResultInfo::ELEMENT;
  
  //const CF::UInt maxNumSrcEntities = scrMap_->GetNumEntities();
  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);

  std::cout << "\t\t Obtaining target entities " << std::endl;
  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);
  //const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);
    
  std::cout << "\t\t\t Interpolator is dealing with " << numSrcEntities <<
               " source " << (useElemAsSource_ ? "elements" : "nodes") << 
               " and "<< numTrgEntities << " target " << (useElemAsTarget_ ? "elements" : "nodes") << std::endl;
  FillMatrix(globSrcEntity, globTrgEntity);
  std::cout << "\t\t Interpolation prepared!" << std::endl;
}


ResultIdList AbstractInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void AbstractInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("NearesNeighbour interpolation required input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Interpolator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  useElemAsSource_ = resultManager_->GetDefOn(upResIds[0]) == ExtendedResultInfo::ELEMENT;
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);
  //after this filter we have nodal values on different regions
  //on a different grid
  if(useElemAsTarget_){
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  } else {
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  }
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}

void AbstractInterpolator::PutIntoMatrix(UInt out, UInt in, Double factor) {
  UInt index = matrix_.outInIndex[out];
  for (;matrix_.outIn[index] != UnusedEntityNumber; index++) {}
  matrix_.outIn[index] = in;
  matrix_.outInFactor[index] = factor;
}

}
