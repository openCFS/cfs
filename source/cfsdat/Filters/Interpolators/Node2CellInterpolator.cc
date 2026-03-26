// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Node2CellInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Jan 3, 2017
 *       \author   kroppert
 */
//================================================================================================


#include "Node2CellInterpolator.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

Node2CellInterpolator::Node2CellInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;

  globalFactor_ = 1.0; // if scaling required, it must be done by dedicated binary operation - filter

}

Node2CellInterpolator::~Node2CellInterpolator(){

}

bool Node2CellInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);
  
	if(resultManager_->GetExtInfo(filterResIds[0])->dType == ExtendedResultInfo::COMPLEX){
	  Vector<Complex>& returnVec = GetOwnResultVector<Complex>(filterResIds[0]);
	  // vector, containing the source data values
	  Vector<Complex>& inVec = GetUpstreamResultVector<Complex>(upResIds[0], stepIndex);
    Node2Cell(returnVec, filterResIds[0], inVec, interpolData_);

    returnVec.ScalarMult(globalFactor_);
	} else {
    Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
    // vector, containing the source data values
    Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

    Node2Cell(returnVec, filterResIds[0], inVec, interpolData_);

    returnVec.ScalarMult(globalFactor_);
	}

  return true;
}






void Node2CellInterpolator::PrepareCalculation(){
  std::cout << "\t ---> Node2CellInterpolator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  // input (source) grid and output (target) grid are the same!
  Grid* inGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;

  std::vector<UInt> allTrgElems;

  StdVector<UInt> curElems;

  //loop over target regions and add element numbers to vector
  std::cout << "\t\t 1/2 Reading elements in target region " << std::endl;
  std::set<std::string>::iterator tRegIter = trgRegions_.begin();
  for(;tRegIter != trgRegions_.end();++tRegIter){
    trgGrid_->GetElemNumsByName(curElems,*tRegIter);
    allTrgElems.insert(allTrgElems.end(),curElems.Begin(),curElems.End());
  }
  std::cout << "\t\t\t Interpolator is dealing with " << allTrgElems.size() <<
               " target element centroids" << std::endl;


  std::cout << "\t\t 2/2 Generating interpolation info ..." << std::endl;
  shared_ptr<EqnMapSimple> upMap = resultManager_->GetEqnMap(upRes);
  interpolData_.reserve(allTrgElems.size());
  StdVector<UInt> tempNodeNums;
  StdVector<UInt> sEqn;
  for(UInt aMatch = 0;aMatch < allTrgElems.size();++aMatch){
      QuantityStruct newStruct;
      //get nodenumbers of containing src-element
      inGrid->GetElemNodes(tempNodeNums, allTrgElems[aMatch]);
      newStruct.trgElemNum = allTrgElems[aMatch];
      newStruct.tNNum.Resize(tempNodeNums.GetSize());
      newStruct.tNNum = tempNodeNums;
      for(UInt n = 0; n < tempNodeNums.GetSize(); ++n){
        upMap->GetEquation(sEqn, tempNodeNums[n], ExtendedResultInfo::NODE);
        for(UInt d = 0; d < sEqn.GetSize(); ++d){
          newStruct.srcEqn.Push_back(sEqn[d]);
        }
      }
      interpolData_.push_back(newStruct);
  }

  allTrgElems.clear();

  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList Node2CellInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void Node2CellInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require node result input
  if(!inInfo->isMeshResult){
    EXCEPTION("Node2Cell interpolation requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Interpolator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  //after this filter we have cell values on different regions
  //on a different grid
  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}
