// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GridIntersectionFilter.cc
 *       \brief    <Description>
 *
 *       \date     Jan 6, 2016
 *       \author   ahueppe
 */
//================================================================================================

#include "GridIntersectionFilter.hh"
#include "Domain/Mesh/MeshUtils/Intersection/VolumeGridIntersection.hh"
#include "Domain/Mesh/MeshUtils/Intersection/IntersectAlgos/TriaIntersect.hh"
#include "Domain/Mesh/MeshUtils/Intersection/IntersectAlgos/ElemIntersect2D.hh"
#include "MatVec/CoordFormat.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "MatVec/CRS_Matrix.hh"
#include <cmath>



namespace CFSDat{

GridIntersectionFilter::GridIntersectionFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;

  globalFactor_ = 1.0; // if scaling required, it must be done by dedicated binary operation - filter

}

GridIntersectionFilter::~GridIntersectionFilter(){

}

bool GridIntersectionFilter::UpdateResults(std::set<uuids::uuid>& upResults) {

  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  if(resultManager_->GetExtInfo(filterResIds[0])->dType == ExtendedResultInfo::COMPLEX){
    /// this is the vector, which will be filled with the result
    Vector<Complex>& returnVec = GetOwnResultVector<Complex>(filterResIds[0]);

    // vector, containing the source data values
    Vector<Complex>& inVec = GetUpstreamResultVector<Complex>(upResIds[0], stepIndex);

    returnVec.Init(0.0);

    this->InterpolationMatrix->MultAdd_type(inVec,returnVec);

    returnVec.ScalarMult(globalFactor_);
  } else {
    /// this is the vector, which will be filled with the result
    Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);

    // vector, containing the source data values
    Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

    returnVec.Init(0.0);

    this->InterpolationMatrix->MultAdd(inVec,returnVec);

    returnVec.ScalarMult(globalFactor_);
  }

  return true;
}

void GridIntersectionFilter::PrepareCalculation(){
  std::cout << "\t ---> GridIntersection preparing for interpolation" << std::endl;


  StdVector<ElemIntersect::VolCenterInfo> intersectInfo;

#ifdef USE_CGAL
  //in this filter we only have one upstream result
  Grid* inGrid   = resultManager_->GetExtInfo(upResIds[0])->ptGrid;

  //first we create the cell intersector
  //YET MISSING make sanity check if all region have same dimension
  UInt intersetionDim = 1;
  std::set<std::string>::iterator reg = trgRegions_.begin();
  intersetionDim = this->trgGrid_->GetEntityDim(*reg);
  std::cout << "\t\t 1. Intersecting Grids" << std::endl;
  if(intersetionDim==2){
    VolumeGridIntersection<ElemIntersect2D> intersect(this->trgGrid_,inGrid,trgRegions_,srcRegions_);
    //TODO: make optional, read from file...
    boost::uuids::uuid tId = this->StartTime();
    intersectInfo = intersect.GetVolCenterInfo();
    std::cout << "\t\t\t Took: " << this->StopTime(tId) << std::endl;
  }else if(intersetionDim==3){
    VolumeGridIntersection<TetraIntersect> intersect(this->trgGrid_,inGrid,trgRegions_,srcRegions_);
    boost::uuids::uuid tId = this->StartTime();
    intersectInfo = intersect.GetVolCenterInfo();
    std::cout << "\t\t\t Took: " << this->StopTime(tId) << std::endl;
  }else{
    EXCEPTION("Unsupported dimension of intersection region")
  }
#else
    EXCEPTION("Intersection needs to be compiled with USE_CGAL=ON!");
#endif



  std::cout << "\t\t 2. Creating interpolation matrix" << std::endl;
  CreateCRS(intersectInfo);
  std::cout << "\t\t\t Interpolation matrix is " << InterpolationMatrix->GetNumRows() << "x" << InterpolationMatrix->GetNumCols()
            << " NNZ: " <<  InterpolationMatrix->GetNnz() << std::endl;

  std::cout << "\t\t 3. Computing Weights" << std::endl;
  FillInterpolationMatrix(intersectInfo);

  std::cout << "\t\t 4. Cleaning Up" << std::endl;
  intersectInfo.Clear();
}

void GridIntersectionFilter::FillInterpolationMatrix(const StdVector<ElemIntersect::VolCenterInfo> & infos){
  //get the equation mapping from the in out results
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(filterResIds[0]);
  str1::shared_ptr<EqnMapSimple> upMap   =   resultManager_->GetEqnMap(upResIds[0]);
  
  InterpolationMatrix->Init();
  UInt negativeCounter = 0;
  UInt nanInfCounter = 0;
#pragma omp parallel reduction(+ : negativeCounter , nanInfCounter) num_threads(CFS_NUM_THREADS)
{
  StdVector<UInt> sElemEq;
  StdVector<UInt> tNodeEq;
  CF::shared_ptr<ElemShapeMap> eShape;
  Vector<Double> localPoint(infos[0].center.GetSize());
  CF::Vector<Double> shFnc;

#pragma omp for
  for(UInt aInfo=0;aInfo<infos.GetSize();++aInfo){
    const UInt& tElem = infos[aInfo].targetElemNum;
    const UInt& sElem = infos[aInfo].sourceElemNum;
    upMap->GetEquation(sElemEq,sElem,ResultInfo::ELEMENT);

    //compute shape function
    const Elem* curTE = trgGrid_->GetElem(tElem);
    // eShape = trgGrid_->GetElemShapeMap(curTE,true);
    // eShape.reset(curTE->ptrShapeMap);
    // trgGrid_->UpdateIndividualElemShapeMap(curTE, true); // const
    // eShape = curTE->ptrShapeMap;
    eShape = ((curTE))->GetElemShapeMap(trgGrid_, true);
    const CF::StdVector<UInt>& tElemConnect = curTE->connect;

    localPoint.Init();
    eShape->Global2Local(localPoint,infos[aInfo].center);

    FeH1 * myElem = dynamic_cast<FeH1*>(eShape->GetBaseFE());
    shFnc.Resize(tElemConnect.GetSize());
    shFnc.Init();
    myElem->GetShFnc(shFnc,localPoint,curTE);

    for(UInt aNode = 0; aNode < tElemConnect.GetSize(); aNode++){
      downMap->GetEquation(tNodeEq,tElemConnect[aNode],ResultInfo::NODE);
      Double curval  = shFnc[aNode] * infos[aInfo].volume;
      if(shFnc[aNode] < 0){
        negativeCounter++;
      }
      if(std::isnan(shFnc[aNode]) || std::isinf(shFnc[aNode])){
        nanInfCounter++;
        shFnc[aNode] = 0.0;
      }
      for(UInt aDOF = 0; aDOF < tNodeEq.GetSize(); ++aDOF){
        InterpolationMatrix->AddToMatrixEntry(tNodeEq[aDOF],sElemEq[aDOF],curval);
      }
    }
}
    if(negativeCounter > 0){
      std::cerr << "Detected " << negativeCounter << " negative weights. This could indicate errors. Check your results!" << std::endl;
    }
    if(nanInfCounter > 0){
      std::cerr << "Detected " << nanInfCounter << " nan/inf weights. This indicate errors. Setting those contributions to Zero!" << std::endl;
    }
  }
}

void GridIntersectionFilter::CreateCRS(const StdVector<ElemIntersect::VolCenterInfo> & infos){

  //get the equation mapping from the in out results
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(filterResIds[0]);
  str1::shared_ptr<EqnMapSimple> upMap   =   resultManager_->GetEqnMap(upResIds[0]);
//  uuids::uuid upRes = upResIds[0];
//  Grid* sGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;

  //TODO this is an open Question
  UInt numRows = downMap->GetNumEquations();
  UInt numCols = upMap->GetNumEquations();

  StdVector< std::set<UInt> > connecting(numRows);
  StdVector<UInt> tElemConnect;
  StdVector<UInt> tNodeEq;
  StdVector<UInt> sElemEq;

  //make connectivity unique
  for(UInt aInfo=0;aInfo<infos.GetSize();++aInfo){
    const UInt& tElem = infos[aInfo].targetElemNum;
    const UInt& sElem = infos[aInfo].sourceElemNum;
    trgGrid_->GetElemNodes(tElemConnect,tElem);
    upMap->GetEquation(sElemEq,sElem,ResultInfo::ELEMENT);
    for(UInt aNode = 0; aNode < tElemConnect.GetSize(); aNode++){
      downMap->GetEquation(tNodeEq,tElemConnect[aNode],ResultInfo::NODE);
      for(UInt aDOF = 0; aDOF < tNodeEq.GetSize(); ++aDOF){
        connecting[tNodeEq[aDOF]].insert(sElemEq[aDOF]);
      }
    }
  }

  //count NNZ
  UInt nnz = 0;
  for(UInt aC=0;aC<connecting.GetSize();++aC){
    nnz += connecting[aC].size();
  }

  //fill container and convert to CRS
  CoordFormat<Double>* myContainer = new CoordFormat<Double>(numRows,numCols,nnz,false);
  for(UInt aC=0;aC<connecting.GetSize();++aC){
    std::set<UInt>::iterator aIter = connecting[aC].begin();
    for(;aIter!= connecting[aC].end();++aIter){
      myContainer->AddEntry(aC,*aIter,0.0);
    }
  }

  myContainer->FinaliseAssembly();

  //create CRS
  InterpolationMatrix = new CF::CRS_Matrix<Double>(*myContainer,false);

  //delete container
  connecting.Clear();
  delete myContainer;
}


ResultIdList GridIntersectionFilter::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void GridIntersectionFilter::AdaptFilterResults(){
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
