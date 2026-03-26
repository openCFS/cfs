// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     L2norm.hh
 *       \brief    <Description>
 *
 *       \date     Feb 13, 2017
 *       \author   sschoder
 */
//================================================================================================


#include "L2norm.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

L2norm::L2norm(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;
}

L2norm::~L2norm(){

}

bool L2norm::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);
  

  Double returnVal = 0.0;

  //current element
  UInt curE;
  CF::StdVector<UInt> eqns;
  shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(filterResIds[0]);


  // for every element in the target mesh
  for(UInt i=0;i < interpolData_.size();++i){
    InpolationStruct& aStru = interpolData_[i];
    curE = aStru.trgElemNum;

    //Sum up the contributions from all nodes of element curE
    StdVector<UInt> sEqn;
    downMap->GetEquation(eqns,curE,ExtendedResultInfo::ELEMENT);


    //UInt dim = eqns.GetSize(); //dimension of input data (scalar, 2vector, 3vector)
    for(UInt aNode =0;aNode < aStru.tNNum.GetSize(); ++aNode){
      for(UInt aDof=0;aDof < eqns.GetSize(); aDof++){
        returnVec[eqns[aDof]] += inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]]
                                       * inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]]
                                               * aStru.volume
                                               / aStru.tNNum.GetSize();
        returnVal += inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]]
                                       * inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]]
                                               * aStru.volume
                                               / aStru.tNNum.GetSize();
      }
    }
  }

  std::cout << "L2 Norm" << returnVal << std::endl;
  
  return true;
}





void L2norm::PrepareCalculation(){
  std::cout << "\t ---> L2norm preparing for interpolation" << std::endl;

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
      InpolationStruct newStruct;
      //get nodenumbers of containing src-element
      inGrid->GetElemNodes(tempNodeNums, allTrgElems[aMatch]);
      shared_ptr<ElemShapeMap> eShape = trgGrid_->GetElemShapeMap(inGrid->GetElem(aMatch+1),true);
      newStruct.volume = eShape->CalcVolume();
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

ResultIdList L2norm::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void L2norm::AdaptFilterResults(){
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
