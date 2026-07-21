// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalInterp.cc
 *       \brief    Implementation File
 *
 *       \date     Jan. 20, 2013
 *       \author   Andreas Hueppe
 */
//================================================================================================



#include <fstream>
#include <algorithm>

#include "CoefFunctionGridNodalInterp.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/Results/ResultInfo.hh"
#include "MatVec/CoordFormat.hh"
#include "MatVec/CRS_Matrix.hh"
#include "MatVec/BaseMatrix.hh"
#include "FeBasis/FeSpace.hh"
#include "FeBasis/FeFunctions.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "General/Enum.hh"
#include "Driver/BaseDriver.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Timer.hh"

namespace CoupledField{

// declare class specific logging stream
DEFINE_LOG(coeffunctiongridnodalinterp, "coeffunctiongridnodalinterp")

template<typename DATA_TYPE>
CoefFunctionGridNodalInterp<DATA_TYPE>::
CoefFunctionGridNodalInterp(Domain* ptDomain,
                            PtrParamNode configNode, PtrParamNode curInfo, shared_ptr<RegionList> regions)
                            :CoefFunctionGridNodal<DATA_TYPE>(ptDomain, configNode, regions){
  //right now hardcoded for conservative
  this->stdInterpReady_ = false;
  this->consInterpReady_ = false;
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;
  this->extDataInfo_ = curInfo->Get("externalGrid",ParamNode::APPEND);
  this->destGrid_ = this->domain_->GetGrid();
  ReadXMLNode(configNode);

  // For not interpolated nodes
  if (configNode->Has("notInterpolatedNodesFile")) {
    notInterpolatedNodesFileName_ = configNode->Get("notInterpolatedNodesFile")->As<std::string>();
    std::cout << "Going to write not interpolated node numbers to file " << notInterpolatedNodesFileName_ << std::endl;
  }
  
  // Determine which steps are available
  this->domain_->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);

  //====================================================
  // Create CoefFunction info
  //====================================================
  if(this->resultInfo_->entryType == ResultInfo::SCALAR){
    this->dimDof_ = 1;
    this->dimType_ = CoefFunction::SCALAR;
  }else if(this->resultInfo_->entryType == ResultInfo::VECTOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::VECTOR;
  }else if(this->resultInfo_->entryType == ResultInfo::TENSOR){
    this->dimDof_ = this->resultInfo_->dofNames.GetSize();
    this->dimType_ = CoefFunction::TENSOR;
  }
  UInt gDim = this->destGrid_->GetDim();
  this->CreateOperator(gDim,this->dimDof_);
  this->SetDestRegions(regions);
  this->InitSolVec();
  this->WriteGlobalFactorsToXML(configNode);
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::PrepareInterpolation(){

  if(!this->stdInterpReady_){
    std::cout << "++ Preparing for interpolation of external data...";
    std::cout.flush();
    //====================================================
    // Create Data structures for easy solution access
    //====================================================
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    this->PrepareForStdInterp(this->myConfigNode_);
    std::cout << "Done" << std::endl;
    std::cout.flush();
  }else{
    EXCEPTION("CoefFunctionGridNodalInterp<DATA_TYPE>::PrepareInterpolation() This should not happen!!");
  }

}




template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                                               const LocPointMapped& lpm ){
  if(this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE){
    EXCEPTION("GetVector invalid for conservative interpolation");
  }else if (this->curInterpType_ == CoefFunctionGrid::NO_INTERPOLATION){
    WARN("Interpolation type was neither set by user nor internally. Setting to standard...");

    std::string interpTyStr = CoefFunctionGrid::InterpType_.ToString(this->curInterpType_);
    this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpTyStr);
    this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("CFS++");
  }
  EXCEPTION("Get Tensor is not implemented");
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                                               const LocPointMapped& lpm ){
  if(this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE){
    EXCEPTION("GetVector invalid for conservative interpolation");
  }else if (this->curInterpType_ == CoefFunctionGrid::NO_INTERPOLATION){
    WARN("Interpolation type was neither set by user nor internally. Setting to default...");
    std::string interpTyStr = CoefFunctionGrid::InterpType_.ToString(this->curInterpType_);
    this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpTyStr);
    this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("CFS++");
  }

  assert(this->dimType_ != CoefFunction::TENSOR);

  //two possibilities 1. we have already the intepolation function, we use it.
  //2. no interpolation function there, we do the old stuff
  if(this->dimType_ == CoefFunction::SCALAR)
    CoefMat.Resize(1);
  else if (this->dimType_ == CoefFunction::VECTOR)
    CoefMat.Resize(this->dimDof_);

  //if this is the first time we call the procedure, std interpolation is not ready
  //so we prepare it
  if(!this->stdInterpReady_){
    EXCEPTION("CoefFunctionGridNodalInterp<DATA_TYPE>::GetVector You are trying to read in the external grid again!!");
  }else{
    if(this->dependType_ != CoefFunction::SPACE){
      if(this->UpdateSolution())
        this->interpolFunction_->ApplyExternalData();
    }
  }

  
  //there is a special case when dealing with surface elements
  //we need to obtain a valid volume from them
  const Elem * curE = NULL;
  LocPointMapped curPt;
  if(lpm.isSurface){
    //we want to operate only on volume elements!
    curE = lpm.lpmVol->ptEl;
    curPt = *(lpm.lpmVol);
  }else{
    curE = lpm.ptEl;
    curPt = lpm;
  }
  Vector<DATA_TYPE> elemSol;
  this->interpolFunction_->GetElemSolution(elemSol,curE);
  BaseFE * ptFe = this->interpolFunction_->GetFeSpace()->GetFe(curE->elemNum);
  this->myOperator_->ApplyOp(CoefMat,curPt,ptFe,elemSol);
  
  //This was the obsolete approach of searching for each point in external grid
  //the comment will be removed soon, when the new concept has been verified
  
  //CoefMat.Init();
  ////const Elem* targetElem = lpm.ptEl;
  //Vector<Double> globCoord;
  //Vector<Double> localCoord;
  //Vector<DATA_TYPE> elemSol;
  //lpm.shapeMap->Local2Global(globCoord,lpm.lp);
  //LocPoint lp;
  ////build up set of source regions
  //std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  //std::set<RegionIdType> scrRegIds;
  //for( ; regIter != this->srcRegions_.end(); ++regIter) {
  //  RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
  //  scrRegIds.insert(curId);
  //}
  //const Elem* sourceElem = this->srcGrid_->GetElemAtGlobalCoord(globCoord,lp,scrRegIds);
  //
  //if(!sourceElem){
  //  WARN("Could not find element for coord " << globCoord);
  //  return;
  //}
  //
  //shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  //LocPointMapped lpmS;
  //lpmS.Set(lp,esm,1.0);
  //
  //this->GetElemSolution( elemSol, sourceElem->elemNum);
  //BaseFE * ptFe = esm->GetBaseFE();
  //this->myOperator_->ApplyOp(CoefMat,lpmS,ptFe,elemSol);
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                                              const LocPointMapped& lpm ){
  if(this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE){
    EXCEPTION("GetVector invalid for conservative interpolation");
  }else if (this->curInterpType_ == CoefFunctionGrid::NO_INTERPOLATION){
    WARN("Interpolation type was neither set by user nor internally. Setting to standard...");
    std::string interpTyStr = CoefFunctionGrid::InterpType_.ToString(this->curInterpType_);
    this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpTyStr);
    this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("CFS++");
  }

  assert(this->dimType_ == CoefFunction::SCALAR );

  Vector<DATA_TYPE> tmpVec;
  this->GetVector(tmpVec,lpm);
  CoefMat = tmpVec[0];
}



template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::SetDestRegions(shared_ptr<RegionList> regions){
  
  StdVector<RegionIdType> regIDs = regions->GetRegionIds();
  for (UInt i = 0; i < regIDs.GetSize(); i++) {
    this->destRegions_.insert(regIDs[i]);
  }
  
  std::stringstream ss;
  for (UInt i = 0; i < this->entities_.GetSize(); i++) {
    std::string name = this->entities_[i]->GetName();
    ss << name;
    if (i < this->entities_.GetSize() - 1) {
      ss << ",";
    }
    this->extDataInfo_->Get("DestinationRegionList")->Get("Region",ParamNode::APPEND)->Get("name")->SetValue(name);
  }
  
  this->destRegionNames_ = ss.str();
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::MapElemNodesConservative(){
  //first obtain the node coordinates of the source regions
  StdVector< Vector<Double> > nodeGlobCoords;
  StdVector<LocPoint> localCoords;
  StdVector<const Elem*> foundElements;

  StdVector<UInt> nodeNums;
  std::set<std::string>::iterator regIt = this->srcRegions_.begin(),
                                  endIt = this->srcRegions_.end();
  if (this->srcRegions_.size() == 1) {
    RegionIdType aReg = this->srcGrid_->GetRegion().Parse(*regIt);
    this->srcGrid_->GetNodesByRegion(nodeNums, aReg);
  } else {
    const UInt maxNumNode = this->srcGrid_->GetNumNodes();
    std::vector<bool> usedNum(maxNumNode,false);
    while (regIt != endIt){
      RegionIdType aReg = this->srcGrid_->GetRegion().Parse(*regIt);
      nodeNums.Clear();
      this->srcGrid_->GetNodesByRegion(nodeNums, aReg);
      const UInt regionNodeNum = nodeNums.GetSize();
      for (UInt i = 0; i < regionNodeNum; i++) {
        usedNum[nodeNums[i]] = true;
      }
      ++regIt;
    }
    nodeNums.Clear();
    for (UInt i = 0; i <= maxNumNode; i++) {
      if (usedNum[i]) {
        nodeNums.Push_back(i);
      }
    }
  }

  UInt numNodes = nodeNums.GetSize();
  UInt srcDim = this->srcGrid_->GetDim();
  UInt destDim = destGrid_->GetDim();
  nodeGlobCoords.Reserve(numNodes);
  Vector<Double> aCoord;
  for (UInt i=0; i<numNodes; ++i) {
    this->srcGrid_->GetNodeCoordinate(aCoord, nodeNums[i]);
    if (srcDim > destDim) {
      if (abs(aCoord[2]-xyPlaneAtZ_) > zTol_) {
        continue;
      }
      aCoord.Resize(destDim);
    }
    nodeGlobCoords.Push_back(aCoord);
  }
  if (numNodes == 0) {
    EXCEPTION("There are no nodes for interpolation at z = " << xyPlaneAtZ_
              << " m +/- " << zTol_ << " m.");
  }
  
  StdVector<shared_ptr<EntityList> > lists;
  std::set<RegionIdType>::const_iterator destRegIt = this->destRegions_.begin();
  for(; destRegIt != this->destRegions_.end(); ++destRegIt ) {
    shared_ptr<ElemList> newList(new ElemList(this->destGrid_));
    newList->SetRegion(*destRegIt);
    lists.Push_back(newList);
  }

  destGrid_->GetElemsAtGlobalCoords( nodeGlobCoords,
                                    localCoords,
                                    foundElements,
                                    lists,
                                    globalTol_,
                                    localTol_,
                                    false);

  //we now create some debugging information
  UInt elemCounter = 0;
  for (UInt i=0; i<numNodes; ++i) {
    if (!foundElements[i]) {
      ++elemCounter;
    }
    else {
      LOG_DBG3(coeffunctiongridnodalinterp) << "Node #" << nodeNums[i] << " get associated to element #" << foundElements[i]->elemNum << std::endl;
      LOG_DBG3(coeffunctiongridnodalinterp)<< "Local Coordinate is: " << localCoords[i] << std::endl << std::endl;
    }
  }
  if (elemCounter > 0) {
    WARN("There were " << elemCounter << " unmapped nodes from source region(s) \'" << this->allSrcRegionNames_ << "\' which are not mapped to region \'" << this->destRegionNames_ << "\'. Perhaps you should increase the tolerances!");
    if(this->verbose_ == true){
      PrintNodesToCSV(foundElements,nodeGlobCoords);
    }
  }

  if (notInterpolatedNodesFileName_ != "") {
    std::ofstream outFile;
    outFile.open(notInterpolatedNodesFileName_.c_str(), std::ofstream::out );
    if (outFile.fail()) {
      std::cout << "Writing not interpolated node numbers failed" << std::endl;
    } else {
      std::cout << "Writing not interpolated node numbers to file " << notInterpolatedNodesFileName_ << std::endl;
      for (UInt i=0; i<numNodes; ++i) {
        if (!foundElements[i]) {
          UInt notInt = i + 1;
          outFile << notInt << std::endl;
        }
      }
      outFile.close();
      std::cout << "Closed file " << notInterpolatedNodesFileName_ << std::endl;
    }
  }
  
  this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("numUnmappedNodes")->SetValue(elemCounter);
  this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("globalTol")->SetValue(globalTol_);
  this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("localTol")->SetValue(localTol_);

  //if the user wants to, we save the node->element association here
  //now we store the information in our map
  for (UInt i=0; i<numNodes; ++i) {
    //check if we have the information
    if (foundElements[i]) {
      std::vector<UInt> & curVec = this->elemNodeAssoc_[foundElements[i]->elemNum];
      curVec.push_back(nodeNums[i]);
      this->localCoordsNodeAssoc_[nodeNums[i]] = localCoords[i];
    }
  }
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::PrintNodesToCSV(const StdVector<const Elem*>& foundElements,
                                                                     const StdVector< Vector<Double> >& nodeGlobCoords){
  //RegionIdType regID = this->destGrid_->GetRegion().Parse(this->destRegionName_);
  std::string s = this->allSrcRegionNames_;
  std::replace( s.begin(), s.end(), ',', '_');
  std::string filename = "unmappedNodes_" + s + ".csv";
  std::cerr << "Printing unmapped node coordinates to file: " << filename << std::endl;
  std::ofstream outFile(filename.c_str(), std::ios::out | std::ios::trunc);
  outFile << "x coord, y coord, z coord, scalar\n";
  for (UInt i=0; i<foundElements.GetSize(); ++i) {
    if (!foundElements[i]) {
       const Vector<Double>& coord = nodeGlobCoords[i];
       outFile << coord[0] << ", " << coord[1];
       if(coord.GetSize()==3)
         outFile << ", " << coord[2];
       //outFile << ", " << regID << "\n";
    }
  }
  outFile.close();
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::ReadXMLNode(PtrParamNode configNode){
  //extract input and grid ids
  this->inputId_ = configNode->Get("inputId")->As<std::string>();
  this->gridId_ = configNode->Get("gridId")->As<std::string>();

  //this was already determined in constructor of CoefFunctionGrid
  //this->aSeqStep_ = this->domain_->GetDriver()->GetActSequenceStep();

  this->verbose_ = configNode->Get("verbose")->As<bool>();

  this->DetermineResult(this->inputId_,this->aSeqStep_);

  this->extDataInfo_->Get("interpolation")->Get("sourceGridID")->SetValue(this->gridId_);
  this->extDataInfo_->Get("interpolation")->Get("sourceInputID")->SetValue(this->inputId_);

  std::string interpStr = configNode->Get("interpolation")->As<std::string>();
  this->curInterpType_ = CoefFunctionGrid::InterpType_.Parse(interpStr);
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpStr);
  this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("XML");

  globalTol_ = 0.0;
  configNode->GetValue("globalTolerance", globalTol_, ParamNode::PASS);
  localTol_ = 1e-3;
  configNode->GetValue("localTolerance", localTol_, ParamNode::PASS);

  if ( this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE ) {
    xyPlaneAtZ_ = 0.0;
    zTol_ = 1e-12;
    if (configNode->Has("xyPlane")) {
      if (destGrid_->GetDim() == 3) {
        WARN("xyPlane directive is ignored, because target grid of interpolation is 3D.");
      }
      else {
        configNode->Get("xyPlane")->GetValue("z", xyPlaneAtZ_, ParamNode::EX);
        configNode->Get("xyPlane")->GetValue("tolerance", zTol_, ParamNode::PASS);
      }
    }
  }
  
  //obtain grid pointer
  this->srcGrid_ = this->domain_->GetGrid(this->gridId_);

  //determine source region list
  StdVector<PtrParamNode> regs = configNode->Get("regionList")->GetList("region");
  std::stringstream ss;
  for(UInt i=0;i<regs.GetSize();++i){
    std::string srcRegName = regs[i]->Get("srcRegionName")->As<std::string>();
    this->srcRegions_.insert(srcRegName);
    this->extDataInfo_->Get("SourceRegionList")->Get("Region",ParamNode::APPEND)->Get("name")->SetValue(srcRegName);
    ss << srcRegName;
    if (i < regs.GetSize() - 1) {
      ss << ",";
    }
  }
  this->allSrcRegionNames_ = ss.str();
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::MapConservative( shared_ptr<FeSpace> targetSpace,
                                                                      Vector<DATA_TYPE>& feFncVec){
  this->UpdateSolution();
  if(!this->consInterpReady_){
    std::cout << "++ Preparing for conservative interpolation of external data... ";
    std::cout.flush();
    shared_ptr<Timer> t(new Timer);
    t->Start();

    //ok so we need to gather the formation
    this->MapElemNodesConservative();
    CoordFormat<DATA_TYPE>* myContainer;

    //determine the number of rows in the CRS matrix
    StdVector<Integer> elemEqns;
    UInt nnz=0;
    UInt numRows= targetSpace->GetNumFreeEquations();
    UInt numCols = this->solVec_.GetSize();
    std::map<UInt,std::vector<UInt> >::iterator eIt = this->elemNodeAssoc_.begin();
    // determine number of non zeros (nnz)
    while(eIt != this->elemNodeAssoc_.end()){
      const Elem* curE = this->destGrid_->GetElem(eIt->first);
      targetSpace->GetElemEqns(elemEqns,curE);
      for(UInt i=0;i<elemEqns.GetSize();++i){
        if(elemEqns[i]>0){
          nnz += eIt->second.size();
        }
      }
      //nnz += fe->GetNumFncs()*eIt->second.size()*this->dimDof_;
      ++eIt;
    }
    
    // caring for verbose of interpolated sources
    if (this->verboseSum_) {
      UInt size = feFncVec.GetSize();
      this->countSolvecIndex_.resize(size,false);
      if (this->dimDof_ > 1) {
        this->solVecIndexDim_.Resize(size,0);
      }
    }

    //ok lets create the container
    myContainer = new CoordFormat<DATA_TYPE>(numRows,numCols,nnz,false);

    shared_ptr<ElemShapeMap> esm;
    Matrix<Double> opMat;
    eIt = this->elemNodeAssoc_.begin();
    LocPointMapped lpm;
    while(eIt != this->elemNodeAssoc_.end()){
      const Elem* curE = NULL;
      curE = this->destGrid_->GetElem(eIt->first);
      targetSpace->GetElemEqns(elemEqns,curE);
      esm = this->destGrid_->GetElemShapeMap( curE, true );
      //now we add for each source node the corresponding entries
      for(UInt aNode=0; aNode < eIt->second.size(); aNode++){
        UInt curNodeNum = eIt->second[aNode];
        LocPoint lp = this->localCoordsNodeAssoc_[curNodeNum];
        try{
          lpm.Set(lp,esm,1.0);
        }catch(...){
          WARN("Found negative Jacobian for current local point: " + lp.coord.ToString() + ". Setting to element mid-point.")
          lp.coord = Elem::shapes[curE->type].midPointCoord;
          lpm.Set(lp,esm,1.0);
        }
        BaseFE * fe = targetSpace->GetFe(curE->elemNum);
        this->myOperator_->CalcOpMatTransposed(opMat,lpm,fe);
        for(UInt j=0;j<fe->GetNumFncs();j++){
          for(UInt d =0;d<this->dimDof_;++d){
            if(elemEqns[j*this->dimDof_+d]>0){
              UInt curSrcEq =  curNodeNum*this->dimDof_+d;
              myContainer->AddEntry(elemEqns[j*this->dimDof_+d]-1,curSrcEq,opMat[j*this->dimDof_+d][d]);
              if (this->verboseSum_) {
                this->countSolvecIndex_[elemEqns[j*this->dimDof_+d]-1] = true;
                if (this->dimDof_ > 1) {
                  this->solVecIndexDim_[elemEqns[j*this->dimDof_+d]-1] = d;
                }
              }
            }
          }
        }//add to coordMatrix
      }//for each source node
      ++eIt;
    }//for each target element
    myContainer->FinaliseAssembly();
    //TODO: delete unnecessary data structures


    //now we create a CRS_Matrix from it
    this->consInterpMat_.reset(new CRS_Matrix<DATA_TYPE>(*myContainer,false));

    //delete coordinate matrix
    delete myContainer;
    //ready to go for conservative interpolation
    this->consInterpReady_ = true;
    this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("matrix")->Get("Created")->SetValue("yes");
    this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("matrix")->Get("nnz")->SetValue(nnz);
    this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("matrix")->Get("numRow")->SetValue(numRows);
    this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("matrix")->Get("numCol")->SetValue(numCols);
    t->Stop();
    std::stringstream elapsed;
    const int walltime((int) t->GetWallTime());
    if(walltime > 120) {
      const int wallmin((int) (walltime / 60.0));
      if(wallmin > 60){
        elapsed << wallmin/60 << "h";
      }else{
        elapsed << wallmin << "min";
      }
    }else{
      elapsed << walltime << "s";
    }

    this->extDataInfo_->Get("interpolation")->Get("conservative")->
        Get("creationTime")->SetValue(elapsed.str());

    std::cout << " done." << std::endl;
    std::cout.flush();
  }

  //here it gets simple we just take the external solution vector and multiply with matrix
  //std::cout << this->consInterpMat_->GetNumRows() << " " << this->consInterpMat_->GetNumCols() << std::endl;
  //std::cout << this->solVec_.GetSize() << " " << feFncVec.GetSize() << std::endl;
  //this->consInterpMat_->MultAdd(this->solVec_,feFncVec);
  this->consInterpMat_->Mult(this->solVec_,feFncVec);
  
  if (this->verboseSum_) {
    UInt size = feFncVec.GetSize();
    if (this->dimDof_ > 1) {
      StdVector<DATA_TYPE> sum(this->dimDof_);
      sum.Init(0.0);
#pragma omp parallel
      {
        StdVector<DATA_TYPE> tSum(this->dimDof_);
        tSum.Init(0.0);
#pragma omp for
        for (Integer i = 0; i < (Integer) size; ++i) {
          if (this->countSolvecIndex_[i]) {
            tSum[this->solVecIndexDim_[i]] += feFncVec[i];
          }
        }
#pragma omp critical
        {
          for (UInt d = 0; d < this->dimDof_; ++d) {
            sum[d] += tSum[d];
          }
        }
      }
      std::cout << "Sum of " << this->solName_ << " interpolated from " << this->allSrcRegionNames_ << " to " 
        << this->destRegionNames_ << ": ";
      for (UInt d = 0; d < this->dimDof_; ++d) {
        std::cout << sum[d];
        if (d < this->dimDof_ - 1) {
          std::cout << ", ";
        }
      }
      std::cout << std::endl;
    } else {
      DATA_TYPE sum = 0.0;
#pragma omp parallel
      {
        DATA_TYPE tSum = 0.0;
#pragma omp for
        for (Integer i = 0; i < (Integer) size; ++i) {
          if (this->countSolvecIndex_[i]) {
            tSum += feFncVec[i];
          }
        }
#pragma omp critical
        {
          sum += tSum;
        }
      }
      std::cout << "Sum of " << this->solName_ << " interpolated from " << this->allSrcRegionNames_ << " to " 
        << this->destRegionNames_ << ": " << sum << std::endl;
    }
  }
}

//Prepare for standard interpolation
template<class DATA_TYPE>
  void CoefFunctionGridNodalInterp<DATA_TYPE>::PrepareForStdInterp(PtrParamNode motherNode){
  /* We create two paramnodes for configuration of the interpolation space:
   FOR POLYNOMIALS
   <Lagrange id="<random>_polyId" useGridOrder="true" spectral="false">
       <isoOrder>1</isoOrder>
   </Lagrange>

   For INTEGRATION
   <scheme id=<random>_integId">
     <method>GAUSS</method>
     <order>1</order>
     <mode>absolute</mode>
   </scheme>
  */
  
  //now we generate a fake integID random string such that nobody will ever use this...
  std::string polyId;
  std::string integId;
  
  std::string charset = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
  
  std::stringstream result;
  
  //boost::random::random_device rng;
  boost::random::mt19937 gen;
  boost::random::uniform_int_distribution<unsigned int> index_dist(0, charset.size() - 1);
  for(int i = 0; i < 12; ++i) {
    result  << charset[index_dist(gen)];
  }
  
  polyId = result.str() + "_poly";
  integId = result.str() + "_integ";
  
  //now we obtain the paramnode pointer for the xml and add our new nodes
  //in future versions we will read those parameters from xml node
  PtrParamNode polyMotherNode = this->myConfigNode_->GetRoot()->
      Get("fePolynomialList", ParamNode::INSERT );
  polyNode_ = polyMotherNode->GetByVal("Lagrange","id",polyId, ParamNode::INSERT);
  
  polyNode_->Get("spectral", ParamNode::INSERT)->SetValue("false");
  polyNode_->Get("gridOrder",ParamNode::INSERT);
  //polyNode_->Get("isoOrder", ParamNode::INSERT)->SetValue("2");
  //comment: actually this would only make sense for a field interpolation,
  // for a conservative interpolation this does not result in a better quality
  
  PtrParamNode integMotherNode = this->myConfigNode_->GetRoot()->
      Get("integrationSchemeList", ParamNode::INSERT );
  integNode_ = integMotherNode->GetByVal("scheme","id",integId,ParamNode::INSERT);
  integNode_->Get("method", ParamNode::INSERT)->SetValue("Lobatto");
  integNode_->Get("order", ParamNode::INSERT)->SetValue("1");
  integNode_->Get("mode", ParamNode::INSERT)->SetValue("absolute");
  
  StdVector<PtrParamNode> regs = motherNode->Get("regionList")->GetList("region");
  for(UInt i=0;i<regs.GetSize();++i){
    regs[i]->Get("polyId",ParamNode::INSERT)->SetValue(polyId);
    regs[i]->Get("integId",ParamNode::INSERT)->SetValue(integId);
  }
  
  //with this information we can now create a nice FeSpace
  PtrParamNode interpolSpaceNode = this->extDataInfo_->Get("interpolation")->Get("standard")->Get("interpolationSpace");
  shared_ptr<FeSpace> interpolSpace;
  interpolSpace = FeSpace::CreateInstance(motherNode,interpolSpaceNode,FeSpace::H1,this->destGrid_);
  
  shared_ptr<SolStrategy> solStrat(new SolStrategyStd(motherNode));
  
  interpolSpace->Init(solStrat);
  this->interpolFunction_.reset(new FeFunction<DATA_TYPE>(this->domain_->GetMathParser()));
  this->interpolFunction_->SetFeSpace(interpolSpace);
  this->interpolFunction_->SetGrid(this->destGrid_);
  interpolSpace->AddFeFunction(interpolFunction_);
  this->interpolFunction_->SetResultInfo(this->resultInfo_);
  
  //loop over the regions of the coeffunction
  //now we do a little trick, we pass ourselves to the FeFunction and let the function fill itselve
  if(this->entities_.GetSize() > 1)
    EXCEPTION("Right now, we only support a single entity list here");
  
  for(UInt aReg=0;aReg<this->entities_.GetSize();aReg++){
    this->interpolFunction_->AddEntityList(this->entities_[aReg]);
    interpolSpace->SetRegionApproximation(this->entities_[aReg]->GetRegion(),polyId,integId);
    this->interpolFunction_->AddExternalDataSource(this->shared_from_this(), this->entities_);
    //now check if this is a surface region
    //if so, we add the associated volume region
    if(this->entities_[aReg]->GetType() == EntityList::SURF_ELEM_LIST ||
        this->entities_[aReg]->GetType() == EntityList::NC_ELEM_LIST ){
      EntityIterator it = this->entities_[aReg]->GetIterator();
      RegionIdType volRegionId = it.GetSurfElem()->ptVolElems[0]->regionId;
      std::string volReg = this->destGrid_->GetRegion().ToString(volRegionId);
      shared_ptr<EntityList> VolList = this->destGrid_->GetEntityList( EntityList::ELEM_LIST, volReg);
      this->interpolFunction_->AddEntityList(VolList);
      interpolSpace->SetRegionApproximation(volRegionId,polyId,integId);
    }
  }
  this->interpolFunction_->SetFctId(PSEUDO_FCT_ID);
  
  interpolSpace->Finalize();
  this->interpolFunction_->Finalize();
  
  //for time varing grid Data, we would need to do this in each timestep...
  this->interpolFunction_->ApplyExternalData();
  
  this->stdInterpReady_ = true;
}



template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                                                      StdVector< Vector<DATA_TYPE> >& values,
                                                                      Grid* ptGrid,
                                                                      const StdVector<shared_ptr<EntityList> >& srcEntities,
                                                                      bool updatedGeo )
{
  if(!this->stdInterpReady_){
    std::cout << "Preparing for interpolation of external data...";
    std::cout.flush();
    //====================================================
    // Create Data structures for easy solution access
    //====================================================
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    //    this->PrepareForStdInterp(this->myConfigNode_);
    std::cout << "Done" << std::endl;
    std::cout.flush();
  }
  
  if(localCoords_.GetSize() == 0){
    this->srcGrid_->GetElemsAtGlobalCoords( globCoord,
                                             localCoords_,
                                             foundElements_,
                                             StdVector<shared_ptr<EntityList> >(),this->globalTol_,this->localTol_);
  }

  Vector<DATA_TYPE> eSol;
  Matrix<DATA_TYPE> opMat;
  LocPointMapped lpm;
  shared_ptr<ElemShapeMap> esm;
  values.Resize(globCoord.GetSize(), Vector<DATA_TYPE>(this->dimDof_));
  for(UInt i=0;i<foundElements_.GetSize();++i){
    const Elem* curE = foundElements_[i];
    if(!curE){
      continue;
    }
    this->GetElemSolution(eSol,curE->elemNum);
    esm = this->srcGrid_->GetElemShapeMap( curE, true );
    LocPoint lp = localCoords_[i];
    lpm.Set(lp,esm,1.0);
    BaseFE * fe = esm->GetBaseFE();
    this->myOperator_->CalcOpMat(opMat,lpm,fe);
    values[i] = opMat * eSol;
  }
  //release memory in case of constant data
  if(this->dependType_ == CoefFunction::SPACE){
    localCoords_.Clear();
    foundElements_.Clear();
  }
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                                                      StdVector< DATA_TYPE >& values,
                                                                      Grid* ptGrid,
                                                                      const StdVector<shared_ptr<EntityList> >& srcEntities )
{
  StdVector< Vector<DATA_TYPE> > vecValues;
  this->GetVectorValuesAtCoords(globCoord,vecValues, ptGrid, srcEntities);
  values.Resize(globCoord.GetSize(),0.0);
  for(UInt i=0;i<vecValues.GetSize();++i){
    values[i] = vecValues[i][0];
  }
}

template class CoefFunctionGridNodalInterp<Double>;
template class CoefFunctionGridNodalInterp<Complex>;

} // namespace CoupledField
