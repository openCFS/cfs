// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Cell2NodeInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Apr 20, 2016
 *       \author   sschoder
 */
//================================================================================================


#include "Cell2NodeInterpolator.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>

namespace CFSDat{

Cell2NodeInterpolator::Cell2NodeInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;

  globalFactor_ = 1.0; // if scaling required, it must be done by dedicated binary operation - filter

}

Cell2NodeInterpolator::~Cell2NodeInterpolator(){

}

bool Cell2NodeInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);
  
  Cell2Node(returnVec,filterResIds[0], inVec, interpolData_, nodeNeighbours_);

  returnVec.ScalarMult(globalFactor_);


  return true;
}


void Cell2NodeInterpolator::PrepareCalculation(){
  //1. Get Cell points from input
  //3. Search for containing elements in trg
  //4. Store for each src cell local Coordinates, src cell Idx, trg cell idx

  std::cout << "\t ---> Cell2NodeInterpolator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  Grid* inGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;

  //lets declare some variables and estimate the memory
  std::vector<UInt> allSrcElems;
  CF::StdVector<const CF::Elem*> trgElements;
  CF::StdVector< LocPoint > locPoints;
  CF::StdVector< CF::Vector<Double> > elemCentroids;

  //loop over source regions and add element numbers to vector

  std::set<std::string>::iterator sRegIter = srcRegions_.begin();
  for(;sRegIter != srcRegions_.end();++sRegIter){
    StdVector<UInt> curElems;
    inGrid->GetElemNumsByName(curElems,*sRegIter);
    allSrcElems.insert(allSrcElems.end(),curElems.Begin(),curElems.End());
  }


  std::cout << "\t\t 1/6 Obtaining source element centroids " << std::endl;
  StdVector<shared_ptr<EntityList> > lists;
  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType aReg = trgGrid_->GetRegion().Parse(*destRegIt);
    shared_ptr<ElemList> newList(new ElemList(trgGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);
  }
  std::cout << "\t\t\t Interpolator is dealing with " << allSrcElems.size() <<
               " source element centroids" << std::endl;
  //should not be necessary to make it unique
  elemCentroids.Resize(allSrcElems.size());
  locPoints.Resize(allSrcElems.size());
  for(UInt i=0;i<allSrcElems.size();++i){
    CF::Vector<Double> cCoord;
    inGrid->GetElemCentroid(cCoord,allSrcElems[i],true);
    if(trgGrid_->GetDim() == 2){
      elemCentroids[i].Resize(2);
      elemCentroids[i][0] = cCoord[0];
      elemCentroids[i][1] = cCoord[1];
    }else{
      elemCentroids[i].Resize(3);
      elemCentroids[i][0] = cCoord[0];
      elemCentroids[i][1] = cCoord[1];
      elemCentroids[i][2] = cCoord[2];
    }
   // std::cout << elemCentroids[i].GetSize() << std::endl;
  }

  std::cout << "\t\t 2/6 Searching for containing target elements (can take a while)..." << std::endl;
  trgGrid_->GetElemsAtGlobalCoords(elemCentroids,locPoints,trgElements,
                                   lists,1e-6, 1e-3);

  std::cout << "\t\t 3/6 Generating interpolation info ..." << std::endl;
  interpolData_.reserve(trgElements.GetSize());
  UInt foundCounter = 0;
  for(UInt aMatch = 0;aMatch < trgElements.GetSize();++aMatch){
    if(trgElements[aMatch]!= NULL){
      //obtain element volume
      QuantityStruct newStruct;
      shared_ptr<ElemShapeMap> eShape = trgGrid_->GetElemShapeMap(trgElements[aMatch],true);
      newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.srcEqnSingle = allSrcElems[aMatch];
      newStruct.trgElemNum = trgElements[aMatch]->elemNum;
      interpolData_.push_back(newStruct);
      ++foundCounter;
    }
  }

  std::cout << "\t\t\t Number of interpolation pairs computed: " << foundCounter << std::endl;
  std::cout << "\t\t 4/6 Clear generated temporary data storage ..." << std::endl;
  trgElements.Clear(false);
  elemCentroids.Clear(false);
  locPoints.Clear(false);
  allSrcElems.clear();

  //for an export import step, here would be the right place

  std::cout << "\t\t 5/6 Remap data to equation numbers ..." << std::endl;

    shared_ptr<EqnMapSimple> upMap = resultManager_->GetEqnMap(upRes);
    CF::StdVector<UInt> sEqn;


  trgGrid_->SetNodesToElemsMap();
  for(UInt i=0;i<interpolData_.size();++i){
    upMap->GetEquation(sEqn,interpolData_[i].srcEqnSingle,ExtendedResultInfo::ELEMENT);
    //save, assuming a scalar type
    interpolData_[i].srcEqnSingle = sEqn[0];

    QuantityStruct& aStru = interpolData_[i];

    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);
    const CF::StdVector<UInt>& eConn = curE->connect;

    UInt neighbourElems;
    UInt node;
    StdVector<RegionIdType>  volRegions;
    trgGrid_->GetVolRegionIds(volRegions);
    nodeNeighbours_.Resize(eConn.GetSize()*interpolData_.size());
    for(UInt aNode =0;aNode < eConn.GetSize(); ++aNode){
      node = eConn[aNode];
      trgGrid_->GetNumOfElemsNextToNodes(neighbourElems,node,volRegions);
      nodeNeighbours_[eConn[aNode]] = neighbourElems;
      }

  }

  std::cout << "\t\t 6/6 Sort Data according to eqn numbers ..." << std::endl;
  std::sort(interpolData_.begin(),interpolData_.end());

  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList Cell2NodeInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void Cell2NodeInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("Cell2Node interpolation required input to be defined on mesh");
  }
  //require defined on elems
  if(inInfo->definedOn != ExtendedResultInfo::ELEMENT){
    EXCEPTION("Cell2Node interpolation can only handle element results");
  }
  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Interpolator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);
  //after this filter we have nodal values on different regions
  //on a different grid
  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}
