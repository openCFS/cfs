// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GradientDifferentiator.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "GradientDifferentiator.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>
#include <unordered_map>


namespace CFSDat{

GradientDifferentiator::GradientDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;

  epsScal_ = params_->Get("RBF_Settings")->Get("epsilonScaling")->As<Double>();
  if( params_->Get("RBF_Settings")->Has("betaScaling") ){
	  betaScal_ = params_->Get("RBF_Settings")->Get("betaScaling")->As<Double>();
  }else{
	  betaScal_ = 0.0;
  }

  if( params_->Get("RBF_Settings")->Has("kScaling") ){
	  kScal_ = params_->Get("RBF_Settings")->Get("kScaling")->As<Double>();
  }else{
	  kScal_ = 0.0;
  }

  logEps_ = params_->Get("RBF_Settings")->Get("logEps")->As<bool>();

}

GradientDifferentiator::~GradientDifferentiator(){

}

bool GradientDifferentiator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);

  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;
  StdVector< StdVector<CF::UInt> >& sourceM = matrix.targetSourceIndex;

  UInt tDim = trgMap_->GetNumEqnPerEnt();

  CalcGradient(returnVec, inVec, numEquPerEnt_, sourceM, targetSourceFactor, maxNumTrgEntities, tDim);

  return true;
}


CF::StdVector<GradientDifferentiator*> GradientDifferentiator::differentiators_;
CF::StdVector<GradientDifferentiator::Matrix> GradientDifferentiator::matrices_;


void GradientDifferentiator::PrepareCalculation(){
  std::cout << "\t ---> GradientDifferentiator preparing for interpolation" << std::endl;


  std::cout << "\t\t 1/3 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  scrMap_ = resultManager_->GetEqnMap(upRes);
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  if(numEquPerEnt_ > 1){
    EXCEPTION("Gradient of a vector is not implemented yet!");
  }

  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  // the following map is needed because the inVec in Run() doesn't know
  // which nodeNumber belongs to which entry...we also can't hardcode it
  // because we have dynamical storage of the neighbours, means we can not
  // predict the size
  std::unordered_map<UInt, UInt> sEnt;
  for(UInt i = 0; i < globSrcEntity.GetSize(); ++i){
    sEnt[globSrcEntity[i]] = i + 1;
  }

  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);


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
               " source " << (inElems ? "elements" : "nodes") << " and "<< numTrgEntities << " target elements" << std::endl;


  std::cout << "\t\t 3/3 Creating interpolation matrix ... this can take quite a while ... " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  StdVector< StdVector<CF::UInt> >& sourceM = matrix.targetSourceIndex;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;

  targetSourceFactor.Resize(numTrgEntities);
  sourceM.Resize(numTrgEntities);

  StdVector<RegionIdType> rId;
  rId.Init(0);
  StdVector<UInt> numNodesTrgRegs;
  numNodesTrgRegs.Init(0);
  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType r = trgGrid_->GetRegion().Parse(*destRegIt);
    rId.Push_back(r);
    UInt cache = trgGrid_->GetNumNodes(r);
    numNodesTrgRegs.Push_back(cache);
  }

  StdVector<RegionIdType> src_rId;
  src_rId.Init(0);
  StdVector<UInt> numNodesSrcRegs;
  numNodesSrcRegs.Init(0);
  UInt noRegions = 0;
  std::set<std::string>::const_iterator srcRegIt = this->srcRegions_.begin();
  for(; srcRegIt != this->srcRegions_.end(); ++srcRegIt ) {
    RegionIdType r = inGrid_->GetRegion().Parse(*srcRegIt);
    src_rId.Push_back(r);
    UInt cache = inGrid_->GetNumNodes(r);
    numNodesSrcRegs.Push_back(cache);
    noRegions += 1;
  }

  for(UInt i=0; i < noRegions; i++){
	    if(numNodesSrcRegs[i]!=numNodesTrgRegs[i]){
	      std::cout << "\t\t\t Differentiator is dealing with " << numNodesSrcRegs[i] <<
	                   " source " << (inElems ? "elements" : "nodes") << " and "<< numNodesTrgRegs[i] << " target nodes" << std::endl;
	      EXCEPTION("Source and target mesh of involved regions must be consistent. Use a subsequent interpolation filter, if required otherwise.");
	    }
  }

  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    CF::UInt globEntityNumber;
        globEntityNumber = globTrgEntity[trgEnt];
        if (globEntityNumber != UnusedEntityNumber) {
          CF::Vector<Double> trgCoord;
          StdVector<UInt> nodeList;
          StdVector<CF::Elem*> elemList;
          StdVector<UInt> nList;

          inGrid_->GetElemCentroid(trgCoord, globEntityNumber,true);
          inGrid_->GetElemNodes(nodeList, globEntityNumber);
          if(rId.GetSize() == 0) EXCEPTION("REGION - OpenMP - Grid Problem")
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
          while( !CalcLocGradient(tsF, trgCoord, maxd, srcDist, neighbourCoords, numSrcPoints,
                               numEquPerEnt_, inGrid_, epsScal_, betaScal_, kScal_, logEps_)){
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

ResultIdList GradientDifferentiator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void GradientDifferentiator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("GradientDifferentiator requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Differentiator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    std::cout<<("============================================================ \n"
                "Input of GradientDifferentiator-filter is defined on elements!! \n"
                "This works but it's not very accurate, especially at the boundary \n"
                "You better interpolate the element-values to nodes (e.g. Cell2Node) \n"
                "and differentiate afterwards\n"
                "============================================================")<<std::endl;
    EXCEPTION("GradientDifferentiator requires input to be defined on nodes");
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);

  CF::StdVector<std::string> dofnames;
  dofnames.Push_back("x");
  dofnames.Push_back("y");
  if (trgGrid_->GetDim() == 3){
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

