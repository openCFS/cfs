// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CentroidInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Dec 1, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "CentroidInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>

namespace CFSDat{

CentroidInterpolator::CentroidInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;

  globalFactor_ = 1.0; // if scaling required, it must be done by dedicated binary operation - filter

  checkSum_ = false;
  if(config->Has("sourceSum")){
    checkSum_ = config->Get("sourceSum")->As<bool>();
  }


}

CentroidInterpolator::~CentroidInterpolator(){

}

bool CentroidInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  if(resultManager_->GetExtInfo(filterResIds[0])->dType == ExtendedResultInfo::COMPLEX){
    Vector<Complex>& returnVec = GetOwnResultVector<Complex>(filterResIds[0]);

    // vector, containing the source data values
    Vector<Complex>& inVec = GetUpstreamResultVector<Complex>(upResIds[0], stepIndex);

    //perform interpolation
    returnVec.Init(0.0);
    this->InterpolationMatrix->MultAdd_type(inVec,returnVec);
    returnVec.ScalarMult(globalFactor_);

    } else {
    Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);

    // vector, containing the source data values
    Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

    //perform interpolation
    returnVec.Init(0.0);
    this->InterpolationMatrix->MultAdd(inVec,returnVec);
    returnVec.ScalarMult(globalFactor_);


    // Check filter mesh and output values
    if(checkSum_ == 1){
      Double intSource = returnVec.Sum();
      std::cout<<"Sum over all sources (integrated) = "<<intSource<<std::endl;
    }
  }

  return true;
}

void CentroidInterpolator::PrepareCalculation(){

  std::cout << "\t ---> CentroidInterpolator preparing for interpolation" << std::endl;
  std::cout << "\t\t 1/3 Loading source and target elements " << std::endl;

  Grid* inGrid   = resultManager_->GetExtInfo(upResIds[0])->ptGrid;

  //lets declare source mesh variables and estimate the memory
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

  //lets declare target mesh variables
  StdVector<shared_ptr<EntityList> > lists;
  std::vector<UInt> allTrgElems;

  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType aReg = trgGrid_->GetRegion().Parse(*destRegIt);
    shared_ptr<ElemList> newList(new ElemList(trgGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);

    StdVector<UInt> curElems;
    trgGrid_->GetElemNumsByName(curElems,*destRegIt);
    allTrgElems.insert(allTrgElems.end(),curElems.Begin(),curElems.End());
  }

  std::cout << "\t\t\t Interpolator is dealing with " << allSrcElems.size() <<
               " source element centroids" << std::endl;

  std::cout << "\t\t 2/3 Creating interpolation matrix " << std::endl;
//  CreateCRS(inGrid, trgGrid_);

  std::cout << "\t\t\t Obtaining source element centroids " << std::endl;
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

  std::cout << "\t\t\t Searching for containing target elements (can take a while)..." << std::endl;
  trgGrid_->GetElemsAtGlobalCoords(elemCentroids,locPoints,trgElements,
                                   lists,1e-6, 1e-3);

  //Clean up lists
  lists.Clear(false);

  std::cout << "\t\t\t CRS Matrix setup..." << std::endl;
  //START CRS setup
  //get the equation mapping from the in out results
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(filterResIds[0]);
  str1::shared_ptr<EqnMapSimple> upMap   =   resultManager_->GetEqnMap(upResIds[0]);

  //TODO this is an open Question
  UInt numRows = downMap->GetNumEquations();
  UInt numCols = upMap->GetNumEquations();

  StdVector< std::set<UInt> > connecting(numRows);
  StdVector<UInt> tElemConnect;
  StdVector<UInt> tNodeEq;
  StdVector<UInt> sElemEq;
  UInt foundCounter = 0;

  //Get the coverage of target elements that have a scr point
  Vector<unsigned int> coverage(trgElements.GetSize());

  //make connectivity unique //TODO for all CFD values ...
  for(UInt aInfo=0;aInfo<trgElements.GetSize();++aInfo){
    if(trgElements[aInfo]!= NULL){
      const UInt& tElem = trgElements[aInfo]->elemNum;
      const UInt& sElem = inGrid->GetElem(allSrcElems[aInfo])->elemNum;
      trgGrid_->GetElemNodes(tElemConnect,tElem);
      upMap->GetEquation(sElemEq,sElem,ResultInfo::ELEMENT);
      for(UInt aNode = 0; aNode < tElemConnect.GetSize(); aNode++){
        downMap->GetEquation(tNodeEq,tElemConnect[aNode],ResultInfo::NODE);
        for(UInt aDOF = 0; aDOF < tNodeEq.GetSize(); ++aDOF){
          connecting[tNodeEq[aDOF]].insert(sElemEq[aDOF]);
        }
      }
      coverage[aInfo] = 1;
      ++foundCounter;
    }
  }

  //count NNZ
  UInt nnz = 0;
  std::cout << "\t\t\t\t Count NNZ..." << std::endl;
  for(UInt aC=0;aC<connecting.GetSize();++aC){
    nnz += connecting[aC].size();
  }

  //fill container and convert to CRS
  std::cout << "\t\t\t\t Fill container and convert to CRS..." << std::endl;
  CoordFormat<Double>* myContainer = new CoordFormat<Double>(numRows,numCols,nnz,false);
  for(UInt aC=0;aC<connecting.GetSize();++aC){
    std::set<UInt>::iterator aIter = connecting[aC].begin();
    for(;aIter!= connecting[aC].end();++aIter){
      myContainer->AddEntry(aC,*aIter,0.0);
    }
  }
  std::cout << "\t\t\t\t Finalize Assembly..." << std::endl;
  myContainer->FinaliseAssembly();

  //create CRS
  InterpolationMatrix = new CF::CRS_Matrix<Double>(*myContainer,false);

  //delete container
  connecting.Clear();
  delete myContainer;
  //END CRS setup

  std::cout << "\t\t\t Number of interpolation pairs computed: " << foundCounter << std::endl;
  std::cout << "\t\t\t Coverage of CAA cells by CFD cell: " << coverage.Sum() << " of " << coverage.GetSize() <<std::endl;


  std::cout << "\t\t 3/3 Prepare interpolation ..." << std::endl;
  InterpolationMatrix->Init();
  UInt negativeCounter = 0;
  UInt nanInfCounter = 0;
#pragma omp parallel reduction(+ : negativeCounter , nanInfCounter) num_threads(CFS_NUM_THREADS)
{
  StdVector<UInt> sElemEq;
  StdVector<UInt> tNodeEq;
  CF::shared_ptr<ElemShapeMap> eShape;
  Vector<Double> localPoint(elemCentroids[0].GetSize());
  CF::Vector<Double> shFnc;

#pragma omp for
  for(UInt aInfo=0;aInfo<trgElements.GetSize();++aInfo){
    if(trgElements[aInfo]!= NULL){
      const UInt& tElem = trgElements[aInfo]->elemNum;
      const UInt& sElem = inGrid->GetElem(allSrcElems[aInfo])->elemNum;
      upMap->GetEquation(sElemEq,sElem,ResultInfo::ELEMENT);

      shared_ptr<ElemShapeMap> eShape = inGrid->GetElemShapeMap(inGrid->GetElem(allSrcElems[aInfo]),true);
      Double vol = eShape->CalcVolume();

      //compute shape function
      const Elem* curTE = trgGrid_->GetElem(tElem);
      eShape = trgGrid_->GetElemShapeMap(curTE,true);
      const CF::StdVector<UInt>& tElemConnect = curTE->connect;

      localPoint.Init();
      eShape->Global2Local(localPoint,elemCentroids[aInfo]);

      FeH1 * myElem = dynamic_cast<FeH1*>(eShape->GetBaseFE());
      shFnc.Resize(tElemConnect.GetSize());
      shFnc.Init();
      myElem->GetShFnc(shFnc,localPoint,curTE);

      for(UInt aNode = 0; aNode < tElemConnect.GetSize(); aNode++){
        downMap->GetEquation(tNodeEq,tElemConnect[aNode],ResultInfo::NODE);
        Double curval  = shFnc[aNode] * vol;
        if(shFnc[aNode] < 0){
          negativeCounter++;
          curval = 0.0;
        }
        if(std::isnan(shFnc[aNode]) || std::isinf(shFnc[aNode])){
          nanInfCounter++;
          curval = 0.0;
        }
        for(UInt aDOF = 0; aDOF < tNodeEq.GetSize(); ++aDOF){
          InterpolationMatrix->AddToMatrixEntry(tNodeEq[aDOF],sElemEq[aDOF],curval);
        }
      }
    }
}
    if(negativeCounter > 0){
      WARN("Detected " << negativeCounter << " negative weights. This could indicate errors. Check your results!" << std::endl);
    }
    if(nanInfCounter > 0){
      WARN("Detected " << nanInfCounter << " nan/inf weights. This indicate errors. Setting those contributions to Zero!" << std::endl);
    }
  }

//Clean up
trgElements.Clear(false);
elemCentroids.Clear(false);
locPoints.Clear(false);
allSrcElems.clear();
allTrgElems.clear();

  std::cout << "\t\t Interpolation prepared!" << std::endl;
}




ResultIdList CentroidInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void CentroidInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("Centroid interpolation required input to be defined on mesh");
  }
  //require defined on elems
  if(inInfo->definedOn != ExtendedResultInfo::ELEMENT){
    EXCEPTION("Centroid interpolation can only handle element results");
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
