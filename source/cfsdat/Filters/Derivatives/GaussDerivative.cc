// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GaussDerivative.cc
 *       \brief    <Description>
 *
 *       \date     Apr 23, 2018
 *       \author   mtautz
 */
//================================================================================================

#include "GaussDerivative.hh"

namespace CFSDat{
  
StdVector<GaussDerivative*> GaussDerivative::allDerivatives_;

GaussDerivative::GaussDerivative(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
: BaseFilter(numWorkers,config, resMan){

  // This filter is a filter of the type First In First Out (FIFO)
  this->filtStreamType_ = FIFO_FILTER;
  this->resName = config->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  this->outName = config->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outName);
  upResNames.insert(resName);

}

GaussDerivative::~GaussDerivative(){
}

void GaussDerivative::PrepareCalculation() {
  // checking if needed matrices were already created by another filter
  for (UInt i = 0; i < allDerivatives_.GetSize(); i++) {
    if ((allDerivatives_[i]->grid_ == grid_) && (operator==(*(allDerivatives_[i]->regNames_),*regNames_))) {
      usedDerivative_ = allDerivatives_[i];
      // there exists another Gauss derivative filter using the same matrices
      return;
    }
  }

  Grid::FiniteVolumeRepresentation& fvr = grid_->GetFiniteVolumeRepresentation();
  if (!fvr.isSet) {
    EXCEPTION("No finite volume representation provided on grid for Gauss derivative computation");
  }
  
  // getting some basic information
  const UInt numFaces = fvr.faceCount;
  shared_ptr<EqnMapSimple> scrMap = resultManager_->GetEqnMap(resId);
  
  // check used entities
  std::vector<bool> isFaceInterior(numFaces,false);
  std::vector<bool> isFaceBoundary(numFaces,false);
  std::vector<bool> isFaceSlaveBoundary(numFaces,false);
  StdVector<UInt> faceMasterIdx(numFaces);
  faceMasterIdx.Init(0);
  StdVector<UInt> faceSlaveIdx(numFaces);
  faceSlaveIdx.Init(0);

  // comparing face information to used entities in map and changing accordingly
  for (UInt iFace = 0; iFace < numFaces; iFace++) {
    UInt masterElem = fvr.faceMasterElement[iFace];
    bool isUsedMasterElem = scrMap->IsEntityUsed(masterElem);
    UInt slaveElem = fvr.faceSlaveElement[iFace];
    bool isUsedSlaveElem = scrMap->IsEntityUsed(slaveElem);
    if (fvr.isFaceBoundary[iFace]) {
      if (isUsedMasterElem) {
        isFaceBoundary[iFace] = true;
      }
    } else {
      if (isUsedMasterElem) {
        if (isUsedSlaveElem) {
          isFaceInterior[iFace] = true;
        } else {
          isFaceBoundary[iFace] = true;
        }
      } else {
        if (isUsedSlaveElem) {
          isFaceBoundary[iFace] = true;
          isFaceSlaveBoundary[iFace] = true;
        }
      }
    }
    if (isUsedMasterElem) {
      faceMasterIdx[iFace] = scrMap->GetEntityIndex(masterElem);
    }
    if (isUsedSlaveElem) {
      faceSlaveIdx[iFace] = scrMap->GetEntityIndex(slaveElem);
    }
  }
  
  //const UInt usedElems = revMap.GetSize();
  // creating master elemnts according to mapping
  //StdVector<UInt> revMap = resultManager_->GetEntityNumbers(resId);
  StdVector<UInt> revMap;
  scrMap->GetReverseEntityMap(revMap);
  const UInt usedElems = revMap.GetSize();
  StdVector<UInt> masterElemIdx(usedElems);
  for (UInt i = 0; i < usedElems; i++) {
    UInt elemNum = revMap[i];
    UInt masterElemNum = fvr.masterElement[elemNum];
    if ((elemNum != masterElemNum) && scrMap->IsEntityUsed(masterElemNum)) {
      masterElemIdx[i] = scrMap->GetEntityIndex(masterElemNum);
    } else {
      masterElemIdx[i] = i;
    }
  }
  
  CalcFeVolume(usedElems, masterElemIdx);
  
  std::vector<bool> dummyFBV(numFaces, false);
  CreateGaussComputation(fvr, interiorMatrix_, usedElems, 
      isFaceInterior, dummyFBV, dummyFBV, 
      masterElemIdx, faceMasterIdx, faceSlaveIdx,
      true, hasFvVolume_, fvVolume_);
  
  CreateGaussComputation(fvr, boundaryMatrix_, usedElems, 
      dummyFBV, isFaceBoundary, isFaceSlaveBoundary, 
      masterElemIdx, faceMasterIdx, faceSlaveIdx,
      false, hasFvVolume_, fvVolume_);
  
  for (UInt i = 0; i < usedElems; i++) {
    UInt elemNum = revMap[i];
    UInt masterElemNum = fvr.masterElement[elemNum];
    if ((elemNum != masterElemNum) && scrMap->IsEntityUsed(masterElemNum)) {
      UInt masIdx = masterElemIdx[i];
      fvVolume_[i] = fvVolume_[masIdx];
      hasFvVolume_[i] = hasFvVolume_[masIdx];
      feVolume_[i] = feVolume_[masIdx];
      hasFeVolume_[i] = hasFeVolume_[masIdx];
    }
  }
  
  
  usedDerivative_ = this;
  allDerivatives_.Push_back(this);
}

void GaussDerivative::CalcFeVolume(UInt usedElems, StdVector<UInt>& masterElemIdx) {
  // creating master elemnts according to mapping
  //StdVector<UInt> revMap = resultManager_->GetEntityNumbers(resId);
  StdVector<UInt> revMap;
  shared_ptr<EqnMapSimple> scrMap = resultManager_->GetEqnMap(resId);
  scrMap->GetReverseEntityMap(revMap);
  hasFeVolume_.clear();
  hasFeVolume_.resize(usedElems, false);
  feVolume_.Resize(usedElems);
  feVolume_.Init(0.0);
  for (UInt i = 0; i < usedElems; i++) {
    UInt elemNum = revMap[i];
    const Elem* elem = grid_->GetElem(elemNum);
    shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap(elem);
    Elem::ShapeType st = Elem::GetShapeType(elem->type);
    ElemShape& eShape = Elem::GetShape(st);
    if (eShape.dim == 3) {
      UInt masterElem = masterElemIdx[i];
      feVolume_[masterElem] += esm->CalcVolume();
      hasFeVolume_[masterElem] = true;
    }
  }
}


void GaussDerivative::CreateGaussComputation(Grid::FiniteVolumeRepresentation& fvr, 
      DerivativeMatrix&  derivMatrix, UInt usedElems, 
      std::vector<bool>& isFaceInterior, std::vector<bool>& isFaceBoundary, std::vector<bool>& isFaceSlaveBoundary, 
      StdVector<UInt>& masterElemIdx, StdVector<UInt>& faceMasterIdx, StdVector<UInt>& faceSlaveIdx, 
      bool calcFvVolume, std::vector<bool>& hasFvVolume, StdVector<Double>& fvVolume) {
  const UInt numFaces = fvr.faceCount;

  // setting size of matrix
  derivMatrix.inCount = usedElems;
  derivMatrix.outCount = usedElems;
  if (calcFvVolume) {
    hasFvVolume.clear();
    hasFvVolume.resize(usedElems, false);
    fvVolume.Resize(usedElems);
    fvVolume.Init(0.0);
  }
  
  // computing start indices for matrix rows
  StdVector<UInt>& outInIndex = derivMatrix.outInIndex; 
  outInIndex.Resize(usedElems + 1);
  outInIndex.Init(0);
  std::vector<bool> useSelf(usedElems + 1, false);
  std::vector<bool> useFace(numFaces, false);
  for (UInt iFace = 0; iFace < numFaces; iFace++) {
    if (isFaceInterior[iFace]) {
      outInIndex[faceMasterIdx[iFace]]++;
      outInIndex[faceSlaveIdx[iFace]]++;
      useSelf[faceSlaveIdx[iFace]] = true;
      useSelf[faceMasterIdx[iFace]] = true;
      useFace[iFace] = true;
    } else if (isFaceBoundary[iFace]) {
      useFace[iFace] = true;
      if (isFaceSlaveBoundary[iFace]) {
        useSelf[faceSlaveIdx[iFace]] = true;
      } else {
        useSelf[faceMasterIdx[iFace]] = true;
      }
    }
  }
  for (UInt iElem = 0; iElem < usedElems; iElem++) {
    if (useSelf[iElem]) {
      outInIndex[iElem]++;
    }
  }
  
  // caring for cells splitted into tetrahedrons
  for (UInt iElem = 0; iElem < usedElems; iElem++) {
    if (masterElemIdx[iElem] != iElem) {
      outInIndex[iElem] = outInIndex[masterElemIdx[iElem]];
    }
  }
  
  // creating indices in factor array
  UInt factorIdx = 0;
  for (UInt iElem = 0; iElem <= usedElems; iElem++) {
    UInt iFactors = outInIndex[iElem];
    outInIndex[iElem] = factorIdx;
    factorIdx += iFactors;
  }
  
  StdVector<Double> calcedVolume(usedElems);
  calcedVolume.Init(0.0);
  
  StdVector<UInt>& faceNode = fvr.facePoint;
  StdVector<UInt>& faceNodeIndex = fvr.facePointIndex;
  Grid* grid = fvr.grid;
  StdVector<UInt>& outIn = derivMatrix.outIn; 
  outIn.Resize(factorIdx);
  const UInt dummyIdx = usedElems * 2;
  outIn.Init(dummyIdx);
  StdVector<Double>& outInVector = derivMatrix.outInVector; 
  outInVector.Resize(factorIdx * 3);
  outInVector.Init(0.0);
  for (UInt iFace = 0; iFace < numFaces; iFace++) {
    Vector<Double> finalFaceVector(3, 0.0);
      
    // summing up all face normals by division of polygon into triangles
    UInt nodeIdx = faceNodeIndex[iFace];
    UInt endNodeIdx = faceNodeIndex[iFace + 1];
    UInt numFaceNodes = endNodeIdx - nodeIdx;
    Vector<Double> baseNodeCoord;
    grid->GetNodeCoordinate(baseNodeCoord, faceNode[nodeIdx]);
    Vector<Double> center = baseNodeCoord;
    nodeIdx++;
    
    Vector<Double> coordB;
    grid->GetNodeCoordinate(coordB, faceNode[nodeIdx]);
    center += coordB;
    nodeIdx++;
    Vector<Double> vectorB;
    vectorB = baseNodeCoord - coordB;
    
    Vector<Double> coordC;
    Vector<Double> vectorC;
    Vector<Double> vectorPlus;
    for(; nodeIdx < endNodeIdx; nodeIdx++) {
      grid->GetNodeCoordinate(coordC, faceNode[nodeIdx]);
      center += coordC;
      vectorC = baseNodeCoord - coordC;
      vectorB.CrossProduct(vectorC, vectorPlus);
      finalFaceVector += vectorPlus;
      vectorB = vectorC;
    }
    Double volumeContribution = (0.5 / (Double)numFaceNodes) * (center * finalFaceVector) / 3.0;
    
    // division by 2 to obtain final face normal
    if (isFaceInterior[iFace]) {
      // additional division due to internal face
      finalFaceVector *= 0.25;
    } else {
      finalFaceVector *= 0.5;
    }
    
    UInt masterIdx = faceMasterIdx[iFace];
    UInt slaveIdx = faceSlaveIdx[iFace];
    bool useMaster = true;
    bool useSlave = true;
    if (isFaceBoundary[iFace]) {
      if (isFaceSlaveBoundary[iFace]) {
        useMaster = false;
      } else {
        useSlave = false;
      }
    }
    if (calcFvVolume) {
      if (useMaster) {
        hasFvVolume[masterIdx] = true;
        fvVolume[masterIdx] += volumeContribution;
      }
      if (useSlave) {
        hasFvVolume[slaveIdx] = true;
        fvVolume[slaveIdx] -= volumeContribution;
      }
    }
      
    if (useFace[iFace]) {
      if (useMaster) {
        UInt startIdx = outInIndex[masterIdx];
        UInt writeVecIdx = startIdx * 3;
        outIn[startIdx] = masterIdx;
        outInVector[writeVecIdx] += finalFaceVector[0];
        outInVector[writeVecIdx+1] += finalFaceVector[1];
        outInVector[writeVecIdx+2] += finalFaceVector[2];
        if (useSlave) {
          startIdx++;
          while (outIn[startIdx] != dummyIdx) {
            startIdx++;
          }
          writeVecIdx = startIdx * 3;
          outIn[startIdx] = slaveIdx;
          outInVector[writeVecIdx] += finalFaceVector[0];
          outInVector[writeVecIdx+1] += finalFaceVector[1];
          outInVector[writeVecIdx+2] += finalFaceVector[2];
        }
      }
      if (useSlave) {
        UInt startIdx = outInIndex[slaveIdx];
        UInt writeVecIdx = startIdx * 3;
        outIn[startIdx] = slaveIdx;
        outInVector[writeVecIdx] -= finalFaceVector[0];
        outInVector[writeVecIdx+1] -= finalFaceVector[1];
        outInVector[writeVecIdx+2] -= finalFaceVector[2];
        if (useMaster) {
          startIdx++;
          while (outIn[startIdx] != dummyIdx) {
            startIdx++;
          }
          writeVecIdx = startIdx * 3;
          outIn[startIdx] = masterIdx;
          outInVector[writeVecIdx] -= finalFaceVector[0];
          outInVector[writeVecIdx+1] -= finalFaceVector[1];
          outInVector[writeVecIdx+2] -= finalFaceVector[2];
        }
      }
    }
  }
  
  // Copy coefficients from master to slave elements
  for (UInt iElem = 0; iElem < usedElems; iElem++) {
    UInt masterElem = masterElemIdx[iElem];
    if (masterElem != iElem) {
      UInt masterIdx = outInIndex[masterElem];
      UInt iIdx = outInIndex[iElem];
      if ((outInIndex[masterElem + 1] - masterIdx) != (outInIndex[iElem + 1] - iIdx)) {
        EXCEPTION("What the fuck");
      }
      UInt vectorCount = outInIndex[masterElem + 1] - masterIdx;
      UInt* outInMaster = &outIn[masterIdx];
      UInt* outInI = &outIn[iIdx];
      for (UInt i = 0; i < vectorCount; i++) {
        outInI[i] = outInMaster[i];
      }
      vectorCount *= 3;
      Double* outInVecMaster = &outInVector[masterIdx * 3];
      Double* outInVecI = &outInVector[iIdx * 3];
      for (UInt i = 0; i < vectorCount; i++) {
        outInVecI[i] = outInVecMaster[i];
//        outInVecI[i] = 0.0;
      }
    }
  }
}

template<UInt dim>
void GaussDerivative::VolumeDivide(Double* out) {
  std::vector<bool>& hasVolume = useFvVolume_ ? usedDerivative_->hasFvVolume_ : usedDerivative_->hasFeVolume_;
  StdVector<Double>& volume = useFvVolume_ ? usedDerivative_->fvVolume_ : usedDerivative_->feVolume_;
  const UInt size = volume.GetSize();
  if (dim == 1) {
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < size; i++) {
      if (hasVolume[i]) {
        out[i] /= volume[i];
      }
    }
  } else {
    #pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < size; i++) {
      if (hasVolume[i]) {
        for (UInt d = 0; d < dim; d++) {
          out[i * dim + d] /= volume[i];
        }
      }
    }
  }
}


  
bool GaussDerivative::UpdateResults(std::set<uuids::uuid>& upResults) {
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId);
  Integer stepIndex = resultManager_->GetStepIndex(outId);
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(resId, stepIndex);
  Double* in = &inVec[0];
  Double* out = &returnVec[0];
  
  if (derivType_ == CURL) {
    usedDerivative_->interiorMatrix_.VectorCurl(out, in, false);
    if (!zeroBC_) {
      usedDerivative_->boundaryMatrix_.VectorCurl(out, in, true);
    }
    VolumeDivide<3>(out);
  } else if (derivType_ == DIVERGENCE) {
    if (inputType_ == ResultInfo::TENSOR) {
      usedDerivative_->interiorMatrix_.TensorDivergence(out, in, false);
      if (!zeroBC_) {
        usedDerivative_->boundaryMatrix_.TensorDivergence(out, in, true);
      }
      VolumeDivide<3>(out);
    } else if (inputType_ == ResultInfo::VECTOR) {
      usedDerivative_->interiorMatrix_.VectorDivergence(out, in, false);
      if (!zeroBC_) {
        usedDerivative_->boundaryMatrix_.VectorDivergence(out, in, true);
      }
      VolumeDivide<1>(out);
    }
  } else if (derivType_ == GRADIENT) {
    if (inputType_ == ResultInfo::VECTOR) {
      usedDerivative_->interiorMatrix_.VectorGradient(out, in, false);
      if (!zeroBC_) {
        usedDerivative_->boundaryMatrix_.VectorGradient(out, in, true);
      }
      VolumeDivide<9>(out);
    } else if (inputType_ == ResultInfo::SCALAR) {
      usedDerivative_->interiorMatrix_.ScalarGradient(out, in, false);
      if (!zeroBC_) {
        usedDerivative_->boundaryMatrix_.ScalarGradient(out, in, true);
      }
      VolumeDivide<3>(out);
    }
  }
  //VerboseSum(outId);
  return true;
}


ResultIdList GaussDerivative::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  resId = upResNameIds[*upResNames.begin()];
  outId = filterResIds[0];
  return generated;
}

void GaussDerivative::AdaptFilterResults() {
  //we would need to copy related data to our filter results
  // getting some basic information
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(resId);
  if (inInfo->definedOn != ExtendedResultInfo::ELEMENT) {
    EXCEPTION("Gauss derivative only works with data provided on elements");
  }
  grid_ = inInfo->ptGrid;
  if (grid_->GetDim() != 3) {
    EXCEPTION("Gauss derivative only works with three dimensional meshes");
  }
  regNames_ = inInfo->regNames;
  
  // evaluating derivative type
  inputType_ = inInfo->entryType;
  derivType_ = DIVERGENCE; // will be overwritten or exception occurs
  ResultInfo::EntryType resultEntryType = ResultInfo::VECTOR;
  if (params_->Has("type")) {
    std::string typeString = params_->Get("type")->As<std::string>();
    if (typeString == "Curl") {
      if (inputType_ != ResultInfo::VECTOR) {
        EXCEPTION("Curl only works for vector data");
      }
      derivType_ = CURL;
      resultEntryType = ResultInfo::VECTOR;
    } else if (typeString == "Divergence") {
      if (inputType_ == ResultInfo::SCALAR) {
        EXCEPTION("Divergence does not work for scalars");
      } else if (inputType_ == ResultInfo::VECTOR) {
        resultEntryType = ResultInfo::SCALAR;
      } else {
        resultEntryType = ResultInfo::VECTOR;
      }
      derivType_ = DIVERGENCE;
    } else if (typeString == "Gradient") {
      if (inputType_ == ResultInfo::TENSOR || inputType_ == ResultInfo::TENSOR) {
        EXCEPTION("Gradient is not implemented for tensors");
      } else if (inputType_ == ResultInfo::VECTOR) {
        resultEntryType = ResultInfo::TENSOR;
      } else {
        resultEntryType = ResultInfo::VECTOR;
      }
      derivType_ = GRADIENT;
    } else {
      EXCEPTION("Unknown type given for Gauss derivative, please use either Curl, Divergence or Gradient");
    }
  } else {
    if (inputType_ == ResultInfo::SCALAR) {
      derivType_ = GRADIENT;
      resultEntryType = ResultInfo::VECTOR;
    } else {
      derivType_ = DIVERGENCE;
      if (inputType_ == ResultInfo::VECTOR) {
        resultEntryType = ResultInfo::SCALAR;
      } else {
        resultEntryType = ResultInfo::VECTOR;
      }
    }
  }
  zeroBC_ = params_->Get("zeroBC")->As<bool>();
  useFvVolume_ = params_->Get("useFvVolume")->As<bool>();
  
  resultManager_->CopyResultData(resId,outId);
  resultManager_->SetEntryType(outId,resultEntryType);
  StdVector<std::string> newDofNames;
  if (resultEntryType == ResultInfo::TENSOR) {
    newDofNames.Push_back("e11");
    newDofNames.Push_back("e12");
    newDofNames.Push_back("e13");
    newDofNames.Push_back("e21");
    newDofNames.Push_back("e22");
    newDofNames.Push_back("e23");
    newDofNames.Push_back("e31");
    newDofNames.Push_back("e32");
    newDofNames.Push_back("e33");
  } else if (resultEntryType == ResultInfo::VECTOR) {
    newDofNames.Push_back("x");
    newDofNames.Push_back("y");
    newDofNames.Push_back("z");
  } else if (resultEntryType == ResultInfo::SCALAR) {
    newDofNames.Push_back("-");
  }
  resultManager_->SetDofNames(outId, newDofNames);
  resultManager_->SetValid(outId);
}

}

