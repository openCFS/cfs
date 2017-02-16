// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     RBFInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Sep 8, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "RBFInterpolator.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

RBFInterpolator::RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;
  inDim_ = 0;
  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();

  noSlip_ = false;
  if(config->Has("noSlipWall")){noSlip_ = true;}

}

RBFInterpolator::~RBFInterpolator(){

}

bool RBFInterpolator::Run(){
  // we deactivate every result, except for our own
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();

  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1){
      WARN(" There are still active results when reaching the interpolation filter. This indicates an unexpected use of the pipeline.")
    }
    resultManager_->DeactivateResult(*aIter);
  }
  Double aTF = resultManager_->GetStepValue(filterResIds[0]);
  resultManager_->SetTimeValue(upResIds[0],aTF);
  // now we deactivate our own result and activate the others
  resultManager_->ActivateResult(upResIds[0]);

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }


  CF::StdVector<UInt> eqnNums;

  /// this is the vector, which will be filled with the interpolation result
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init(0.0);

  // vector, containing the source data values
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());


  // vector containing the interpolated result for every trg node
  // can be 1D for scalar values or 2D/3D for vector values
  CF::Vector<Double> vec;


  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  FillScatteredDataVec(scatteredData, vec, inVec, sourceCoords_, inDim_);

  // Perform actual interpolation algorithm
  RBFInterpolation(returnVec, scatteredData, vec, downMap, sourceCoords_, targetCoords_, interpolData_, trgGrid_, inDim_, p_, noSlip_);



  resultManager_->ActivateResult(filterResIds[0]);



  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}




void RBFInterpolator::PrepareCalculation(){
  //1. Get get the source coordinates and the values, defined on those coordinates (Source...Src)
  //2. Get the target coordinates (trg)
  //3. Store for each trg element local Coordinates, elem number, volume, ...
  //4. Small clean-up


  std::cout << "\t ---> RBFInterpolator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  Grid* inGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;


  StdVector<CF::UInt> allSrcElems;
  StdVector<CF::UInt> allSrcNodes;
  //loop over source regions and add element numbers and nodes to vector
  std::set<std::string>::iterator sRegIter = srcRegions_.begin();
  for(;sRegIter != srcRegions_.end();++sRegIter){
    StdVector<CF::UInt> curElems;
    StdVector<CF::UInt> nodeList;
    inGrid->GetElemNumsByName(curElems,*sRegIter);
    inGrid->GetNodesByName(nodeList, *sRegIter);
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
  }


  //loop over target regions and add element numbers to vector and also nodes
  StdVector<const CF::Elem*> allTrgElems;
  StdVector<CF::UInt> allTrgElemNums;
  StdVector<CF::UInt> allTrgNodes;
  StdVector<shared_ptr<EntityList> > lists;
  std::set<std::string>::iterator tRegIter = trgRegions_.begin();
  for(;tRegIter != trgRegions_.end();++tRegIter){
    StdVector<UInt> curElemNums;
    StdVector<UInt> nodeList;

    trgGrid_->GetElemNumsByName(curElemNums,*tRegIter);
    trgGrid_->GetNodesByName(nodeList, *tRegIter);
    //insert the element numbers and elements into the allTrgElems list
    //unfortunately we have to loop over all entries, because in StdVector,
    //only single entries can be "pushed back"
    for (UInt eIter = 0; eIter < curElemNums.GetSize(); ++eIter){
      allTrgElemNums.Push_back(curElemNums[eIter]);

      const CF::Elem* curElem;
      curElem = trgGrid_->GetElem(curElemNums[eIter]);
      allTrgElems.Push_back(curElem);
    }
    //do the same for the nodes
    for (UInt nIter = 0; nIter < nodeList.GetSize(); ++nIter){
      allTrgNodes.Push_back(nodeList[nIter]);
    }
    RegionIdType aReg = trgGrid_->GetRegion().Parse(*tRegIter);
    shared_ptr<ElemList> newList(new ElemList(trgGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);
  }

  std::cout << "\t\t\t Interpolator is dealing with " << allSrcElems.GetSize() <<
      " source elements and "<< allTrgElems.GetSize() << " target elements" << std::endl;


  std::cout << "\t\t 1/5 Obtaining source coordinates " << std::endl;
  //here we also find out, if the input is defined on elements (-centroids) or nodes
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    sourceCoords_.Resize(allSrcElems.GetSize());
    for(UInt i=0;i<allSrcElems.GetSize();++i){
      CF::Vector<Double> cCoord;
      inGrid->GetElemCentroid(cCoord,allSrcElems[i],true);
      if(trgGrid_->GetDim() == 2){
        sourceCoords_[i].Resize(2);
        sourceCoords_[i][0] = cCoord[0];
        sourceCoords_[i][1] = cCoord[1];
      }else{
        sourceCoords_[i].Resize(3);
        sourceCoords_[i][0] = cCoord[0];
        sourceCoords_[i][1] = cCoord[1];
        sourceCoords_[i][2] = cCoord[2];
      }
    }
  }else{
    //that's the case, if the source values are defined on src-nodes
    sourceCoords_.Resize(allSrcNodes.GetSize());
    for(UInt nIter=0; nIter < allSrcNodes.GetSize(); ++nIter){
      CF::Vector<Double> nCoord;
      //inGrid->GetNodeCoordinate3D(nCoord, srcNodes[nIter]);
      inGrid->GetNodeCoordinate3D(nCoord, allSrcNodes[nIter]);
      if(trgGrid_->GetDim() == 2){
        sourceCoords_[nIter].Resize(2);
        sourceCoords_[nIter][0] = nCoord[0];
        sourceCoords_[nIter][1] = nCoord[1];
      }else{
        sourceCoords_[nIter].Resize(3);
        sourceCoords_[nIter][0] = nCoord[0];
        sourceCoords_[nIter][1] = nCoord[1];
        sourceCoords_[nIter][2] = nCoord[2];
      }
    }
  }



  std::cout << "\t\t 2/5 Obtaining target coordinates" << std::endl;

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
  StdVector<const CF::Elem*> tempElems;
  //mapping of global point targetCoords_ to local locPoints
  trgGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints, tempElems, lists, 1e-6, 1e-3);
//  tempElems.Clear();

  std::cout << "\t\t 3/5 Generating interpolation info ..." << std::endl;

  interpolData_.reserve(allTrgElems.GetSize());
  for(UInt aMatch = 0;aMatch < tempElems.GetSize();++aMatch){
    if(tempElems[aMatch]!= NULL){
      QuantityStruct newStruct;
      newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.trgElemNum = tempElems[aMatch]->elemNum;
      interpolData_.push_back(newStruct);
    }
  }

  std::cout << "\t\t 4/5 Clear generated temporary data storage ..." << std::endl;
  allTrgElems.Clear(false);
  allSrcElems.Clear(false);
  allSrcNodes.Clear(false);
  allTrgElemNums.Clear(false);
  allTrgNodes.Clear(false);

  std::cout << "\t\t 5/5 Remap data to equation numbers ..." << std::endl;
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdapter(upRes)->mapping;


  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList RBFInterpolator::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

  //add input result to manager
  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  uuids::uuid newId = resultManager_->AddResult(inRes,this->filterTag_);

  //set the timeline of upstream data if already set
  resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(newId);

  return generated;

}

void RBFInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("RBF interpolation requires input to be defined on mesh");
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
