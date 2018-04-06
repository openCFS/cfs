// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FEBasedInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Jan 18, 2017
 *       \author   sschoder
 */
//================================================================================================


#include "FEBasedInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>

namespace CFSDat{

FEBasedInterpolator::FEBasedInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
}

FEBasedInterpolator::~FEBasedInterpolator(){

}

bool FEBasedInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

  //perform interpolation From Elem result to nodes
  CF::Vector<Double> srcShFnc;
  CF::StdVector<UInt> srcEqns;
  CF::StdVector<UInt> trgEqns;
  CF::shared_ptr<ElemShapeMap> srcShape;
  str1::shared_ptr<EqnMapSimple> srcDownMap = resultManager_->GetEqnMap(upResIds[0]);
  str1::shared_ptr<EqnMapSimple> trgDownMap = resultManager_->GetEqnMap(filterResIds[0]);

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];
  Grid* srcGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;

  returnVec.Init(0.0);
  for(UInt i=0;i < interpolData_.size();++i){

    QuantityStruct& srcStru = interpolData_[i];


    const Elem* srcE = srcGrid_->GetElem(srcStru.srcElemNum);
    srcShape = srcGrid_->GetElemShapeMap(srcE,true);

    const CF::StdVector<UInt>& srcConn = srcE->connect;

    FeH1 * srcElem = dynamic_cast<FeH1*>(srcShape->GetBaseFE());
    //we assume scalar shape functions
    srcShFnc.Resize(srcConn.GetSize());
    srcShFnc.Init();
    srcElem->GetShFnc(srcShFnc,srcStru.localCoords,srcE);

    Double curval = 0.0;
    for(UInt aNode =0;aNode < srcConn.GetSize(); ++aNode){
      srcDownMap->GetEquation(srcEqns,srcConn[aNode],ExtendedResultInfo::NODE);
      trgDownMap->GetEquation(trgEqns,srcStru.tNNum[aNode] ,ExtendedResultInfo::NODE);
      // Finite element based interpolation of a target node point on a src element
      curval  = srcShFnc[aNode];

      for(UInt aDof=0;aDof < srcEqns.GetSize(); aDof++){
        returnVec[trgEqns[aDof]] += curval * inVec[srcEqns[aDof]];
      }
    }
  }

  return true;
}

void FEBasedInterpolator::PrepareCalculation(){
  //1. Get the source values on its nodes
  //2. Search for containing elements in trg


  std::cout << "\t ---> FEBasedInterpolator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  Grid* srcGrid_   = resultManager_->GetExtInfo(upRes)->ptGrid;

  //lets declare some variables and estimate the memory for source Points
  StdVector<CF::UInt> allSrcElems;
  StdVector<CF::UInt> allSrcNodes;
  StdVector<shared_ptr<EntityList> > lists;
  //loop over source regions and add element numbers and nodes to vector
  std::set<std::string>::iterator sRegIter = srcRegions_.begin();
  for(;sRegIter != srcRegions_.end();++sRegIter){
    StdVector<CF::UInt> curElems;
    StdVector<CF::UInt> nodeList;
    srcGrid_->GetElemNumsByName(curElems,*sRegIter);
    srcGrid_->GetNodesByName(nodeList, *sRegIter);
    //insert the element numbers into the allSrcElems list
    //unfortunately we have to loop over all entries, because in StdVector,
    //only single entries can be "pushed back"
    for (UInt eIter = 0; eIter < curElems.GetSize(); ++eIter){
      allSrcElems.Push_back(curElems[eIter]);
    }
    //do the same for the nodes
    for (UInt nIter = 0; nIter < nodeList.GetSize(); ++nIter){
      allSrcNodes.Push_back(nodeList[nIter]);
    }

    RegionIdType aReg = srcGrid_->GetRegion().Parse(*sRegIter);
    shared_ptr<ElemList> newList(new ElemList(srcGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);
  }

  //loop over target regions and add element numbers to vector and also nodes
  StdVector<const CF::Elem*> allTrgElems;
  StdVector<CF::UInt> allTrgNodes;

  std::set<std::string>::iterator tRegIter = trgRegions_.begin();
  for(;tRegIter != trgRegions_.end();++tRegIter){
    StdVector<UInt> nodeList;
    trgGrid_->GetNodesByName(nodeList, *tRegIter);

    //do the same for the nodes
    for (UInt nIter = 0; nIter < nodeList.GetSize(); ++nIter){
      allTrgNodes.Push_back(nodeList[nIter]);
    }

  }
  std::cout << "\t\t\t Interpolator is dealing with " << allSrcNodes.GetSize() <<
               " source nodes and "<< allTrgNodes.GetSize() << " target nodes." << std::endl;

  std::cout << "\t\t 1/3 Obtaining target coordinates" << std::endl;

  // Coordinates of target data
  CF::StdVector< CF::Vector<double> > targetCoords_;
  // target (output) is solely defined on nodes
  targetCoords_.Resize(allTrgNodes.GetSize());
  for(UInt nIter=0; nIter < allTrgNodes.GetSize(); ++nIter){
      CF::Vector<Double> SCoord;
      trgGrid_->GetNodeCoordinate3D(SCoord, allTrgNodes[nIter]);
      if(trgGrid_->GetDim() == 2){
        targetCoords_[nIter].Resize(2);
        targetCoords_[nIter][0] = SCoord[0];
        targetCoords_[nIter][1] = SCoord[1];
      }else{
        targetCoords_[nIter].Resize(3);
        targetCoords_[nIter][0] = SCoord[0];
        targetCoords_[nIter][1] = SCoord[1];
        targetCoords_[nIter][2] = SCoord[2];
      }
    }

  CF::StdVector< LocPoint > locPoints;
  //tempElems are just dummy vectors, get deleted right after we get the local coordinates
  StdVector<const CF::Elem*> tempElems;
  //mapping of global point targetCoords_ to local locPoints
  srcGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints, tempElems, lists, 1e-6, 1e-3);

  std::cout << "\t\t 2/3 Generating interpolation info ..." << std::endl;

    interpolData_.reserve(allTrgNodes.GetSize());
    for(UInt aMatch = 0;aMatch < tempElems.GetSize();++aMatch){
      if(tempElems[aMatch]!= NULL){
        QuantityStruct newStruct;
        newStruct.localCoords = locPoints[aMatch].coord;
        newStruct.srcElemNum = tempElems[aMatch]->elemNum;
        newStruct.tNNum = allTrgNodes[aMatch];
        interpolData_.push_back(newStruct);
      }
    }

  std::cout << "\t\t 3/3 Sort Data according to eqn numbers ..." << std::endl;
  std::sort(interpolData_.begin(),interpolData_.end());

  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList FEBasedInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void FEBasedInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("FEBased interpolation required input to be defined on mesh");
  }
  //require defined on elems
  if(inInfo->definedOn != ExtendedResultInfo::NODE){
    EXCEPTION("FEBased interpolation can only handle node results - If you have Element results use Element2Node filter infront of this filter!");
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
