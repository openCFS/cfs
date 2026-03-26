// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     FEBasedInterpolator.cc
 *       \brief    This is a finte element based interpolator. It works similar as the interpolation
 * in CoefFunctionGridNodel but with extended functionality:
 * 
 *  - the input values are taken either from the nodes or from the elements (centroids)
 * 
 *  - for conservative interpolation the elements in the target grid are found in which the source
 *    points (node coordinates of element centroids) lay in, then the shape function of the element nodes are evaluated at the input point to 
 *    get the interpolation coefficients
 *     - In case of 100% interpolation, also target elements/nodes are searched for, if the sources are not located inside a target element
 *           - if connectivity is on, the neighbours of sources which are evaluated successfully are evaluated and the non-interpolated sources get targets accordingly. The neighbour search is done either based on finite elements (element->node->element) or finite volumes (element->face->element)
 *               - if the distance from source to target exceeds maxDistance only numNeighbours neighbours are taken into account to avoid too big matrices. At such places there should be no relavant data.
 *           - if connectivity is off, numNeigbhours defines how many nearest neighbours are searched for to make the interpolation conservative
 *     - if 100% interpolation and no connectivity is selectet is selected, the nearest nodes to the source points are set for the 
 * 
 *  - for non-conservative interpolation the search is done in reversed order. The target points are
 *    located in the source elements and the shape functions are evaluated from these elements.
 *     - In case of 100% interpolation sorces elements/nodes are searched for, if no source element was found for a target. This can also be done either by connectivity or nearest neighbours
 * 
 * The option exponent regulates the weightings in case of 100& option is needed to find a corresponding source or target. The higher "exponent" is , the higher is the weighting for the near replacements in contrast to the far away ones
 *      
 *
 *       \date     Jun 18, 2018
 *       \author   mtautz
 */
//================================================================================================


#include "FEBasedInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "Domain/Mesh/Grid.hh"
#include "cfsdat/DatUtils/KNNIndexSearch.hh"

#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <limits>

namespace CFSDat{

FEBasedInterpolator::FEBasedInterpolator(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                     :AbstractInterpolator(numWorkers,config,resMan){
  interpolatorName_ = "FEBasedInterpolator";
  
  nSteps_ = 6;
  iSteps_ = 0;
  
  conservative_ = false;
  hundredPercent_ = true;
  connectivity_ = true;
  useFV_ = false;
  numNeighbours_ = 1;
  maxDistance_ = 0.01;
  exponent_ = 1.0;
  if (config->Has("fetype")) {
    CF::PtrParamNode fetype = config->Get("fetype");
    if (fetype->Has("conservative")) {
      conservative_ = fetype->Get("conservative")->As<bool>();
    }
    if (fetype->Has("hundredPercent")) {
      hundredPercent_ = fetype->Get("hundredPercent")->As<bool>();
    }
    if (fetype->Has("connectivity")) {
      connectivity_ = fetype->Get("connectivity")->As<bool>();
    }
    if (fetype->Has("useFV")) {
      useFV_ = fetype->Get("useFV")->As<bool>();
    }
    if (fetype->Has("numNeighbours")) {
      numNeighbours_ = std::abs(fetype->Get("numNeighbours")->As<int>());
    }
    if (fetype->Has("maxDistance")) {
      maxDistance_ = std::abs(fetype->Get("maxDistance")->As<double>());
    }
    if (fetype->Has("exponent")) {
      exponent_ = std::abs(fetype->Get("exponent")->As<double>());
    }
  }
  if (conservative_) {
    verboseSum_ = true;
  }
  if (hundredPercent_) {
    nSteps_ += 2;
  }
}

FEBasedInterpolator::~FEBasedInterpolator(){

}

void FEBasedInterpolator::FillMatrix(StdVector<CF::UInt>& globSrcEntity, StdVector<CF::UInt>& globTrgEntity) {
  UInt maxNumSrcEntities = globSrcEntity.GetSize();
  UInt maxNumTrgEntities = globTrgEntity.GetSize();
  
  matrix_.inCount = maxNumSrcEntities;
  matrix_.outCount = maxNumTrgEntities;
  StdVector<CF::UInt>& outInIndex = matrix_.outInIndex;
  StdVector<CF::UInt>& outIn = matrix_.outIn;
  StdVector<CF::Double>& outInFactor = matrix_.outInFactor;
  outInIndex.Resize(maxNumTrgEntities + 1);
  outInIndex.Init(0);
  
  if (conservative_) {
    // For conservative interpolation we search for a target element for each source entity coordinate
    StdVector<UInt> foundTargetElements;
    StdVector< LocPoint > locCoords;
    std::vector<bool> hasFoundTargetElement;
    //StdVector<CF::UInt> foundReplacementEntity;
    StdVector< StdVector<CF::UInt>* > foundReplacementEntities;
    StdVector< StdVector<CF::Double>* > foundReplacementDistances;
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Seraching target entities for conservative interpolation" << std::endl;
    FindTargetElements(trgGrid_, trgRegions_, globTrgEntity, useElemAsTarget_, trgMap_,
                       inGrid_, srcRegions_, globSrcEntity, useElemAsSource_, scrMap_,
                       foundTargetElements, locCoords, hasFoundTargetElement,
                       foundReplacementEntities, foundReplacementDistances);
    
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Computing interpolation matrix size: ";
    // determining matrix size
    for (UInt i = 0; i < maxNumSrcEntities; i++) {
      UInt sourceEntNumber = globSrcEntity[i];
      if (sourceEntNumber != UnusedEntityNumber) {
        if (hasFoundTargetElement[i]) {
          UInt targetElemNumber = foundTargetElements[i];
          if (useElemAsTarget_) {
            outInIndex[trgMap_->GetEntityIndex(targetElemNumber)]++;
          } else {
            const Elem* elem = trgGrid_->GetElem(targetElemNumber);
            const StdVector<UInt>& elemNodes = elem->connect;
            const UInt numElemNodes = elemNodes.GetSize();
            for (UInt n = 0; n < numElemNodes; n++) {
              outInIndex[trgMap_->GetEntityIndex(elemNodes[n])]++;
            }
          }
        } else if (hundredPercent_) {
          StdVector<CF::UInt>* replacements = foundReplacementEntities[i];
          UInt size = replacements->GetSize();
          for (UInt k = 0; k < size; k++) {
            outInIndex[replacements->operator[](k)]++;
          }
        }
      }
    }
    
    // initializing matrix
    UInt index = 0; 
    for (UInt i = 0; i <= maxNumTrgEntities; i++) {
      UInt plus = outInIndex[i];
      outInIndex[i] = index;
      index += plus;
    }
    std::cout << index << " NNZ" << std::endl;
    outIn.Resize(index);
    outIn.Init(UnusedEntityNumber);
    outInFactor.Resize(index);
    outInFactor.Init(0);
    
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Filling interpolation matrix" << std::endl;
    // filling matrix
    for (UInt i = 0; i < maxNumSrcEntities; i++) {
      UInt sourceEntNumber = globSrcEntity[i];
      if (sourceEntNumber != UnusedEntityNumber) {
        if (hasFoundTargetElement[i]) {
          UInt targetElemNumber = foundTargetElements[i];
          if (useElemAsTarget_) {
            PutIntoMatrix(trgMap_->GetEntityIndex(targetElemNumber), i,  1.0);
          } else {
            const Elem* elem = trgGrid_->GetElem(targetElemNumber);
            shared_ptr<ElemShapeMap> esm = this->trgGrid_->GetElemShapeMap( elem, true );
            FeH1 * trgElemFE = dynamic_cast<FeH1*>(esm->GetBaseFE());
            Vector<Double> ShFns;
            trgElemFE->GetShFnc(ShFns,locCoords[i],elem);
            const StdVector<UInt>& elemNodes = elem->connect;
            const UInt numElemNodes = elemNodes.GetSize();
            for (UInt n = 0; n < numElemNodes; n++) {
              PutIntoMatrix(trgMap_->GetEntityIndex(elemNodes[n]), i, ShFns[n]);
            }
          }
        } else if (hundredPercent_) {
          StdVector<CF::UInt>* replacements = foundReplacementEntities[i];
          StdVector<CF::Double>* distances = foundReplacementDistances[i];
          UInt size = replacements->GetSize();
          if (size == 1) {
            PutIntoMatrix(replacements->operator[](0), i, 1.0);
          } else if (size > 0) {
            Double maxDistance = 0.0;
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              maxDistance = std::max(kDist, maxDistance);
            }
            Double sum = 0.0;
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              sum += std::pow(1.0  - kDist / maxDistance, 2.0);
            }
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              PutIntoMatrix(replacements->operator[](k), i, std::pow(1.0  - kDist / maxDistance, 2.0)/sum);
            }
          }
          delete replacements;
          delete distances;
          //PutIntoMatrix(foundReplacementEntity[i], i, 1.0);
        }
      }
    }
    /**
    // check conservativeness
    StdVector<Double> blasum(maxNumSrcEntities);
    blasum.Init(0.0);
    for (UInt i = 0; i < matrix_.outIn.GetSize(); i++) {
      blasum[outIn[i]] += outInFactor[i];
    }
    UInt maxdiffIdx = 0;
    Double maxDiff = 0.0;
    for (UInt i = 0; i < maxNumTrgEntities; i++) {
      if (globSrcEntity[i] != UnusedEntityNumber) {
        Double diff = std::abs(blasum[i] - 1.0);
        if (diff > maxDiff) {
          maxDiff = diff;
          maxdiffIdx = i;
        }
      }
    }
    std::cout << " maxdiffIdx " << maxdiffIdx << "     maxDiff " << maxDiff << std::endl;
    **/
  } else {
    // For non-conservative interpolation we search for a source element for each target entity
    StdVector<UInt> foundSourceElements;
    std::vector<bool> hasFoundSourceElement;
    //StdVector<CF::UInt> foundReplacementEntity;
    StdVector< StdVector<CF::UInt>* > foundReplacementEntities;
    StdVector< StdVector<CF::Double>* > foundReplacementDistances;
    StdVector< LocPoint > locCoords;
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Seraching source entities for interpolation" << std::endl;
    FindTargetElements(inGrid_, srcRegions_, globSrcEntity, useElemAsSource_, scrMap_,
                       trgGrid_, trgRegions_, globTrgEntity, useElemAsTarget_, trgMap_, 
                       foundSourceElements, locCoords, hasFoundSourceElement,
                       foundReplacementEntities, foundReplacementDistances);
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Computing interpolation matrix size: ";
    // determining matrix size
    for (UInt i = 0; i < maxNumTrgEntities; i++) {
      UInt targetEntNumber = globTrgEntity[i];
      if (targetEntNumber != UnusedEntityNumber) {
        if (hasFoundSourceElement[i]) {
          if (useElemAsSource_) {
            outInIndex[i]++;
          } else {
            UInt sourceElemNumber = foundSourceElements[i];
            const Elem* elem = inGrid_->GetElem(sourceElemNumber);
            outInIndex[i] += elem->connect.GetSize();
          }
        } else if (hundredPercent_) {
          StdVector<CF::UInt>* replacements = foundReplacementEntities[i];
          outInIndex[i] += replacements->GetSize();
        }
      }
    }
    
    // initializing matrix
    UInt index = 0; 
    for (UInt i = 0; i <= maxNumTrgEntities; i++) {
      UInt plus = outInIndex[i];
      outInIndex[i] = index;
      index += plus;
    }
    //std::cout << index << " NNZ" << std::endl;
    outIn.Resize(index);
    outIn.Init(UnusedEntityNumber);
    outInFactor.Resize(index);
  
    // fill matrix
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Filling interpolation matrix";
    for (UInt i = 0; i < maxNumTrgEntities; i++) {
      UInt targetEntNumber = globTrgEntity[i];
      if (targetEntNumber != UnusedEntityNumber) {
        if (hasFoundSourceElement[i]) {
          UInt sourceElemNumber = foundSourceElements[i];
          if (useElemAsSource_) {
            PutIntoMatrix(i, scrMap_->GetEntityIndex(sourceElemNumber), 1.0);
          } else {
            const Elem* elem = inGrid_->GetElem(sourceElemNumber);
            shared_ptr<ElemShapeMap> esm = this->inGrid_->GetElemShapeMap( elem, true );
            FeH1 * trgElemFE = dynamic_cast<FeH1*>(esm->GetBaseFE());
            Vector<Double> ShFns;
            trgElemFE->GetShFnc(ShFns,locCoords[i],elem);
            const StdVector<UInt>& elemNodes = elem->connect;
            const UInt numElemNodes = elemNodes.GetSize();
            for (UInt n = 0; n < numElemNodes; n++) {
              PutIntoMatrix(i, scrMap_->GetEntityIndex(elemNodes[n]), ShFns[n]);
            }
          }
        } else if (hundredPercent_) {
          StdVector<CF::UInt>* replacements = foundReplacementEntities[i];
          StdVector<CF::Double>* distances = foundReplacementDistances[i];
          UInt size = replacements->GetSize();
          if (size == 1) {
            PutIntoMatrix(i, replacements->operator[](0), 1.0);
          } else if (size > 0) {
            Double maxDistance = 0.0;
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              maxDistance = std::max(kDist, maxDistance);
            }
            Double sum = 0.0;
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              sum += std::pow(1.0  - kDist / maxDistance, 2.0);
            }
            for (UInt k = 0; k < size; k++) {
              Double kDist = distances->operator[](k);
              PutIntoMatrix(i, replacements->operator[](k), std::pow(1.0  - kDist / maxDistance, 2.0)/sum);
            }
          }
          delete replacements;
          delete distances;
        }
      }
    }
  }
}

void FEBasedInterpolator::CopyElemeVectorEntries(UInt sourceIdx, UInt targetIdx, StdVector< StdVector<UInt>* >& elements, 
                             std::vector<bool>& hasElementsVector, std::vector<bool>& newAddedElements,
                             UInt& countEls, bool& hasPutSomeNewElements) {
  if (!hasElementsVector[targetIdx]) {
    hasElementsVector[targetIdx] = true;
    elements[targetIdx] = new StdVector<UInt>();
    newAddedElements[targetIdx] = true;
    hasPutSomeNewElements = true;
    countEls++;
  }
  StdVector<UInt>* sourceElems = elements[sourceIdx];
  StdVector<UInt>* targetElems = elements[targetIdx];
  for (UInt ee = 0; ee < sourceElems->GetSize(); ee++) {
    UInt eeI = sourceElems->operator[](ee);
    bool found = false;
    for (UInt nn = 0; nn < targetElems->GetSize() && !found; nn++) {
      found = eeI == targetElems->operator[](nn);
    }
    if (!found) {
      targetElems->Push_back(eeI);
    }
  }
}

void FEBasedInterpolator::FindTargetElements(Grid* findGrid, std::set<std::string>& findRegions, 
                          StdVector<CF::UInt>& findEntitites, bool findElems, shared_ptr<EqnMapSimple> findMap,
                          Grid* searchGrid, std::set<std::string>& searchRegions, 
                          StdVector<CF::UInt>& searchEntitites, bool searchElems, shared_ptr<EqnMapSimple> searchMap,
                          StdVector<CF::UInt>& foundElements, StdVector< LocPoint >& locCoords, std::vector<bool>& hasFoundElement,
                          StdVector< StdVector<CF::UInt>* >& foundReplacementEntities, StdVector< StdVector<CF::Double>* >& foundReplacementDistances) {
  std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Getting Coordinates" << std::endl;
  // getting search coordinates
  UInt numSearchEntities = searchEntitites.GetSize();
  std::vector<bool> needSearchElem(numSearchEntities);
  Vector<double> dummy(3,-1000000.0);
  Vector<double> pCoord;
  StdVector<Vector<double> > globCoords(numSearchEntities);
  for (UInt i = 0; i < numSearchEntities; i++) {
    UInt globEntityNumber = searchEntitites[i];
    if (globEntityNumber != UnusedEntityNumber) {
      if (searchElems) {
        searchGrid->GetElemCentroid(pCoord,globEntityNumber,true);
      } else {
        searchGrid->GetNodeCoordinate3D(pCoord, globEntityNumber);
      }
      globCoords[i] = pCoord;
      needSearchElem[i] = true;
    } else {
      globCoords[i] = dummy;
      needSearchElem[i] = false;
    }
  }
  
  // local coordinates and found Elements
  StdVector< const Elem* > elems;
  
  std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Collecting Elements" << std::endl;
  // entitielsists for the target elements
  StdVector<shared_ptr<EntityList> > lists;
  std::set<std::string>::const_iterator destRegIt = findRegions.begin();
  for(; destRegIt != findRegions.end(); ++destRegIt ) {
    shared_ptr<ElemList> newList(new ElemList(findGrid));
    newList->SetRegion(findGrid->GetRegionId(*destRegIt));
    lists.Push_back(newList);
  }
  
  std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Searching coordinates in elements" << std::endl;
  // searching elements
  findGrid->GetElemsAtGlobalCoords( globCoords, locCoords, elems, lists, 1e-8, 1e-3, false);
  
  // getting element numbers
  foundElements.Resize(numSearchEntities);
  hasFoundElement.resize(numSearchEntities);
  for (UInt i = 0; i < numSearchEntities; i++) {
    if (needSearchElem[i]) {
      if (!elems[i]) {
        foundElements[i] = UnusedEntityNumber;
        hasFoundElement[i] = false;
      } else {
        foundElements[i] = elems[i]->elemNum;
        hasFoundElement[i] = true;
      }
    }
  }
  
  if (!hundredPercent_) {
    return;
  }

  foundReplacementEntities.Resize(numSearchEntities);
  foundReplacementDistances.Resize(numSearchEntities);
  if (!connectivity_) { // nearest neighbour replacement
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Creating search tree for 100 % interpolation" << std::endl;
    // creating search tree
    KNNIndexSearch indexSearch;
    UInt numFindEntities = findEntitites.GetSize();
    for(UInt i = 0; i < numFindEntities; i++) {
      UInt globEntityNumber = findEntitites[i];
      if (globEntityNumber != UnusedEntityNumber) {
        Vector<Double> pCoord;
        if (findElems) {
          findGrid->GetElemCentroid(pCoord,globEntityNumber,true);
        } else {
          findGrid->GetNodeCoordinate3D(pCoord, globEntityNumber);
        }
        indexSearch.AddPoint(pCoord,i);
      }
    }
    indexSearch.BuildTree();
    
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Localizing nearest entites" << std::endl;
    // finding replacement for 100 % interpolation
    for (UInt i = 0; i < numSearchEntities; i++) {
      if (needSearchElem[i] && !hasFoundElement[i]) {
        UInt globEntityNumber = searchEntitites[i];
        Vector<Double> pCoord;
        if (searchElems) {
          searchGrid->GetElemCentroid(pCoord,globEntityNumber,true);
        } else {
          searchGrid->GetNodeCoordinate3D(pCoord, globEntityNumber);
        }
        StdVector<CF::UInt>* replacements = new StdVector<CF::UInt>();
        StdVector<CF::Double>* distances = new StdVector<CF::Double>();
        indexSearch.QueryPoint(pCoord, *replacements, *distances, numNeighbours_);
        foundReplacementEntities[i] = replacements;
        foundReplacementDistances[i] = distances;
      }
    }
  } else { // replacement with use of mesh relations
    if (!searchElems) {
      EXCEPTION("Nonono connectivity not working for elem->* interpolation conservative or *->elem interpolaton non-conservative, please implement it");
    }
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Traversing through mesh connectivity for 100% interpolation" << std::endl;
    StdVector<int> searchRegionIds;
    std::set<std::string>::const_iterator itr = searchRegions.begin();
    for(; itr != searchRegions.end(); ++itr ) {
      searchRegionIds.Push_back(searchGrid->GetRegionId(*itr));
    }
    
    std::vector<bool> hasElementsVector(numSearchEntities,false);
    StdVector< StdVector<UInt>* > elements(numSearchEntities);
    for (UInt i = 0; i < numSearchEntities; i++) {
      if (needSearchElem[i] && hasFoundElement[i]) {
        StdVector<UInt>* iVec = new StdVector<UInt>();
        iVec->Push_back(foundElements[i]);
        elements[i] = iVec;
        hasElementsVector[i] = true;
      }
    }
    std::vector<bool> addedElements = hasElementsVector;
    std::vector<bool> finishedElements = hasElementsVector;
    std::vector<bool> needSearchGlobElem;
    if (useFV_) {
      needSearchGlobElem.resize(searchGrid->GetNumElems()+1,false);
      for (UInt i = 0; i < numSearchEntities; i++) {
        if (needSearchElem[i]) {
          needSearchGlobElem[searchEntitites[i]] = true;
        }
      }
    }
    bool hasPutSomeNewElements = true;
    UInt countEls = 0;
    while (hasPutSomeNewElements) {
      
      countEls = 0;
      hasPutSomeNewElements = false;
      std::vector<bool> newAddedElements(numSearchEntities,false);
      if (useFV_) {
        // copy from element to master element
        Grid::FiniteVolumeRepresentation& fvr = searchGrid->GetFiniteVolumeRepresentation();
        for (UInt iElem = 1; iElem <= searchGrid->GetNumElems(); iElem++) {
          UInt masterElem = fvr.masterElement[iElem];
          if ((masterElem != iElem) && needSearchGlobElem[iElem] && needSearchGlobElem[masterElem]) {
            UInt masterIdx = searchMap->GetEntityIndex(masterElem);
            UInt iIdx = searchMap->GetEntityIndex(iElem);
            if (addedElements[iIdx]) {
              CopyElemeVectorEntries(iIdx, masterIdx, elements, hasElementsVector, addedElements, countEls, hasPutSomeNewElements);
              finishedElements[masterIdx] = true;
            }
          }
        }
        // copy from master element to slave elements
        for (UInt iElem = 1; iElem <= searchGrid->GetNumElems(); iElem++) {
          UInt masterElem = fvr.masterElement[iElem];
          if ((masterElem != iElem) && needSearchGlobElem[iElem] && needSearchGlobElem[masterElem]) {
            UInt masterIdx = searchMap->GetEntityIndex(masterElem);
            UInt iIdx = searchMap->GetEntityIndex(iElem);
            if (addedElements[masterIdx]) {
              CopyElemeVectorEntries(masterIdx, iIdx, elements, hasElementsVector, addedElements, countEls, hasPutSomeNewElements);
              finishedElements[iIdx] = true;
            }
          }
        }
        // grow along faces
        UInt faceCount = fvr.faceCount;
        for (UInt iFace = 0; iFace < faceCount; iFace++) {
          UInt faceMaster = fvr.faceMasterElement[iFace];
          UInt faceSlave = fvr.faceSlaveElement[iFace];
          if ((faceMaster != faceSlave) && needSearchGlobElem[faceMaster] && needSearchGlobElem[faceSlave]) {
            UInt masterIdx = searchMap->GetEntityIndex(faceMaster);
            UInt slaveIdx = searchMap->GetEntityIndex(faceSlave);
            if (addedElements[slaveIdx]) {
              std::swap(masterIdx, slaveIdx);
              std::swap(faceMaster, faceSlave);
            }
            if (addedElements[masterIdx] && !finishedElements[slaveIdx]) {
              CopyElemeVectorEntries(masterIdx, slaveIdx, elements, hasElementsVector, newAddedElements, countEls, hasPutSomeNewElements);
            }
          }
        }
        
      } else {
        for (UInt i = 0; i < numSearchEntities; i++) {
          if (needSearchElem[i] && addedElements[i]) {
            UInt elemNum = searchEntitites[i];
            const Elem* elem = searchGrid->GetElem(elemNum);
            const StdVector<UInt>& eCon = elem->connect;
            StdVector<const Elem*> nextElemsList;
            searchGrid->GetElemsNextToNodes( nextElemsList, eCon, searchRegionIds);
            //StdVector<UInt>* iElemElems = elements[i];
            for (UInt nE = 0; nE < nextElemsList.GetSize(); nE++) {
              UInt nextElemNum = nextElemsList[nE]->elemNum;
              UInt nextElemNumMapIdx = searchMap->GetEntityIndex(nextElemNum);
              if (!finishedElements[nextElemNumMapIdx]) {
                CopyElemeVectorEntries(i, nextElemNumMapIdx, elements, hasElementsVector, newAddedElements, countEls, hasPutSomeNewElements);
              }
            }
          }
        }
      }
      //std::cout << " countEls " << countEls << std::endl;
      for (UInt iif = 0; iif < numSearchEntities; iif++) {
        finishedElements[iif] = finishedElements[iif] || newAddedElements[iif];
      }
      addedElements = newAddedElements;
    }
    searchGrid->ClearNodeToElemConnectivity();
    std::cout << "\t\t\t" << ++iSteps_ << "/" << nSteps_ << " Checking found 100% interpolations" << std::endl;
    for (UInt i = 0; i < numSearchEntities; i++) {
      if (needSearchElem[i] && !hasFoundElement[i]) {
        StdVector<UInt>* replacementElements = elements[i];
        const UInt size = replacementElements->GetSize();
        StdVector<CF::UInt>* replacements = new StdVector<CF::UInt>();
        StdVector<CF::Double>* distances = new StdVector<CF::Double>();
        foundReplacementEntities[i] = replacements;
        foundReplacementDistances[i] = distances;
        Vector<double> pCoord = globCoords[i];
        Double minDistance = std::numeric_limits<Double>::max();
        if (findElems) { // elements are to be found
          replacements->Reserve(size);
          distances->Reserve(size);
          for (UInt j = 0; j < size; j++) {
            replacements->Push_back(findMap->GetEntityIndex(replacementElements->operator[](j)));
            Vector<double> elemCentroid;
            findGrid->GetElemCentroid(elemCentroid, replacementElements->operator[](j), true);
            elemCentroid -= pCoord;
            Double distance = elemCentroid.NormL2();
            distances->Push_back(distance);
            minDistance = std::min(minDistance, distance);
          }
        } else { // nodes are to be found
          std::set<UInt> vertSet;
          for (UInt j = 0; j < size; j++) {
            const Elem* elem = findGrid->GetElem(replacementElements->operator[](j));
            const StdVector<UInt>& eCon = elem->connect;
            for (UInt k = 0; k < eCon.GetSize(); k++) {
              vertSet.insert(eCon[k]);
            }
          }
          UInt vertSize = vertSet.size();
          replacements->Reserve(vertSize);
          distances->Reserve(vertSize);
          for (std::set<UInt>::iterator itr = vertSet.begin(); itr != vertSet.end(); itr++) {
            UInt vertNum = *itr;
            Vector<double> vertCoord;
            findGrid->GetNodeCoordinate3D(vertCoord, vertNum);
            vertCoord -= pCoord;
            replacements->Push_back(findMap->GetEntityIndex(vertNum));
            Double distance = vertCoord.NormL2();
            distances->Push_back(distance);
            minDistance = std::min(minDistance, distance);
          }
        }
        // caring for maximum distance
        if (minDistance > maxDistance_ && distances->GetSize() > numNeighbours_) {
          StdVector< std::pair<Double, UInt> > mypairs;
          for (UInt r = 0; r < distances->GetSize(); r++) {
            mypairs.Push_back(std::pair<Double, UInt>(distances->operator[](r),replacements->operator[](r)));
          }
          std::sort(mypairs.begin(),mypairs.end());
          UInt iNeighbours = numNeighbours_;
          for (; iNeighbours < distances->GetSize() && mypairs[iNeighbours].first == mypairs[iNeighbours - 1].first; iNeighbours++) {}
          distances->Clear(false);
          distances->Resize(iNeighbours);
          replacements->Clear(false);
          replacements->Resize(iNeighbours);
          for (UInt r = 0; r < iNeighbours; r++) {
            distances->operator[](r) = mypairs[r].first;
            replacements->operator[](r) = mypairs[r].second;
          }
        }
      }
    }
    for (UInt i = 0; i < numSearchEntities; i++) {
      if (hasElementsVector[i]) {
        delete elements[i];
      }
    }
  }
  
}

}
