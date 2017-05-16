// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CurlDifferentiator.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include <Filters/Derivatives/CurlDifferentiator.hh>
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>

namespace CFSDat{

CurlDifferentiator::CurlDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;

  epsScal_ = params_->Get("RBF_Settings")->Get("epsilonScaling")->As<Double>();

}

CurlDifferentiator::~CurlDifferentiator(){

}

bool CurlDifferentiator::Run(){
  // we deactivate every result, except for our own
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();

  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1){
      WARN(" There are still active results when reaching the derivative filter. This indicates an unexpected use of the pipeline.")
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

  /// this is the vector, which will be filled with the derivative result
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init();

  // vector, containing the source data values
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);


  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;
  StdVector< StdVector<CF::UInt> >& sourceM = matrix.targetSourceIndex;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;

  UInt gridDim = inGrid_->GetDim();
  CalcCurl(returnVec, inVec, numEquPerEnt_, sourceM, targetSourceFactor, maxNumTrgEntities, gridDim);


  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}

CF::StdVector<CurlDifferentiator*> CurlDifferentiator::differentiators_;
CF::StdVector<CurlDifferentiator::Matrix> CurlDifferentiator::matrices_;


void CurlDifferentiator::PrepareCalculation(){
  std::cout << "\t ---> CurlDifferentiator preparing for interpolation" << std::endl;


  std::cout << "\t\t 1/3 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  scrMap_ = resultManager_->GetResultAdapter(upRes)->mapping;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  if(numEquPerEnt_ == 1){
    EXCEPTION("Curl of a scalar is not defined in this context!");
  }

  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;


  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  // the following map is needed because the inVec in Run() doesn't know
  // which nodeNumber belongs to which entry...we also can't hardcode it
  // because we have dynamical storage of the neighbours, means we can not
  // predict the size
  boost::unordered_map<UInt, UInt> sEnt;
  for(UInt i = 0; i < globSrcEntity.GetSize(); ++i){
    sEnt[globSrcEntity[i]] = i + 1;
  }

  trgMap_ = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  differentiators_.Push_back(this);
  matrixIndex_ = matrices_.GetSize();
  matrices_.Resize(matrixIndex_ + 1);
  Matrix& matrix = matrices_[matrixIndex_];

  std::cout << "\t\t 2/3 Obtaining target entities " << std::endl;

  const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);

  std::cout << "\t\t\t Differentiator is dealing with " << numSrcEntities <<
               " source " << (inElems ? "elements" : "nodes") << " and "<< numTrgEntities << " target nodes" << std::endl;


  std::cout << "\t\t 3/3 Creating interpolation matrix ... this can take quite a while ... " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  StdVector< StdVector<CF::UInt> >& sourceM = matrix.targetSourceIndex;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;

  targetSourceFactor.Resize(numTrgEntities);
  sourceM.Resize(numTrgEntities);


  StdVector<RegionIdType> rId;
  rId.Init(0);
  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType r = inGrid_->GetRegion().Parse(*destRegIt);
    rId.Push_back(r);
  }

//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    CF::UInt globEntityNumber;
        globEntityNumber = globTrgEntity[trgEnt];
        if (globEntityNumber != UnusedEntityNumber) {
          CF::Vector<Double> trgCoord;
          StdVector<UInt> nodeList;
          StdVector<CF::Elem*> elemList;
          StdVector<UInt> nList;
//#pragma omp critical
//          {
          inGrid_->GetElemCentroid(trgCoord, globEntityNumber,true);
//          }
          inGrid_->GetElemNodes(nodeList, globEntityNumber);
          if(rId.GetSize() == 0) EXCEPTION("REGION - OpenMP - Grid Problem")
          //inGrid_->GetElemsNextToNodes(elemList, nodeList, rId);
          //inGrid_->GetNodesOfElemList(nList, elemList, true);
          //for(UInt i = 0; i < nList.GetSize(); ++i){
          //  nodeList.Push_back(nList[i]);
          //}
          StdVector< CF::Vector<CF::Double> > neighbourCoords;
          StdVector<CF::Double> srcDist;
          CF::Vector<CF::Double> tmpCoords;
          StdVector<CF::UInt> sM;
          Double maxd = 0.0;
          for(UInt i = 0; i < nodeList.GetSize(); ++i){
            if(!sM.Contains(sEnt[nodeList[i]])){
              sM.Push_back(sEnt[nodeList[i]]);
              inGrid_->GetNodeCoordinate(tmpCoords, nodeList[i], false);
              neighbourCoords.Push_back(tmpCoords);
              if(tmpCoords.GetSize() == 2) tmpCoords.Push_back(0.0);
              CF::Double d = trgCoord.NormL2(tmpCoords);
              srcDist.Push_back(d);
              if(maxd < d) maxd = d;
            }
          }

          UInt numSrcPoints = srcDist.GetSize();
          CF::Matrix<CF::Double> tsF;
          while( !CalcLocCurl(tsF, trgCoord, maxd, srcDist, neighbourCoords, numSrcPoints,
                               numEquPerEnt_, inGrid_, epsScal_)){
            // find furthest point
            Double d = 0.0;
            UInt maxId = 0;
            for(UInt i = 0; i < srcDist.GetSize(); ++i){
              if(d < srcDist[i]){
                maxId = i;
              }
            }
            sM.Erase(maxId);
            //rowP[trgEnt + 1] -= 1;
            srcDist.Erase(maxId);
            numSrcPoints -= 1;
            neighbourCoords.Erase(maxId);
            if( sM.GetSize() < 3){
              std::cout<<"targetEntity:"<<globEntityNumber<<std::endl;
              std::cout<<"targetCoord: \n"<<trgCoord<<std::endl;
              std::cout<<"neighbourCoords: \n"<<neighbourCoords<<std::endl;
              std::cout<<"distances: \n"<<srcDist<<std::endl;
              std::cout<<"max distances:"<<maxd<<std::endl;
              EXCEPTION("Patch-Problem, modify epsilon!");
            }
          }// while local deriv is false

          targetSourceFactor[trgEnt] = tsF;
          sourceM[trgEnt] = sM;
        }
  }

  std::cout << "\t\t Differentiation prepared!" << std::endl;
}

ResultIdList CurlDifferentiator::SetUpstreamResults(){
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

void CurlDifferentiator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("CurlDifferentiator requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Differentiator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  //input must be node-based
  if (resultManager_->GetDefOn(upResIds[0]) == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    std::cout<<("============================================================ \n"
                "Input of CurlDifferentiator-filter is defined on elements!! \n"
                "This works but it's not very accurate, especially at the boundary \n"
                "You better interpolate the element-values to nodes (e.g. Cell2Node) \n"
                "and differentiate afterwards\n"
                "============================================================")<<std::endl;
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);

  // if we have a 2D-mesh, the output will be a scalar value (per definition, then all
  // differentiation- and aeroacoustic-filters are consistent)
  CF::StdVector<std::string> dofnames;
  if(trgGrid_->GetDim() == 2){
    dofnames.Push_back("x");
    dofnames.Push_back("y");
  }else{
    dofnames.Push_back("x");
    dofnames.Push_back("y");
    dofnames.Push_back("z");
  }
  resultManager_->SetDofNames(filterResIds[0],dofnames);

  //after this filter we have element values on different regions
  //on a different grid
  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}

