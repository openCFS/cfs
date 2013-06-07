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

#include <def_expl_templ_inst.hh>

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
DECLARE_LOG(coeffunctiongridnodalinterp)
DEFINE_LOG(coeffunctiongridnodalinterp, "coeffunctiongridnodalinterp")

template<typename DATA_TYPE>
CoefFunctionGridNodalInterp<DATA_TYPE>::
CoefFunctionGridNodalInterp(Domain* ptDomain,
                            PtrParamNode configNode, PtrParamNode curInfo)
                            :CoefFunctionGridNodal<DATA_TYPE>(ptDomain, configNode){
  //right now hardcoded for conservative

  this->stdInterpReady_ = false;
  this->consInterpReady_ = false;
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;
  this->extDataInfo_ = curInfo->Get("externalGrid");
  ReadXMLNode(configNode);

  //obtain grid pointer and store its dimension
  this->srcGrid_ = this->domain_->GetGrid(this->gridId_);

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
  UInt gDim = this->domain_->GetGrid()->GetDim();
  UInt dDim = this->resultInfo_->dofNames.GetSize();
  this->CreateOperator(gDim,dDim);
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
    std::cout << "Preparing for interpolation of external data...";
    std::cout.flush();
    //====================================================
    // Create Data structures for easy solution access
    //====================================================
    //in this special class we start with the equation mapping right
    //after the first call to get entities as this method should not be called twice!
    this->MapEqns();
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    this->PrepareForStdInterp(this->myConfigNode_);
    std::cout << "Done" << std::endl;
    std::cout.flush();
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
void CoefFunctionGridNodalInterp<DATA_TYPE>::AddEntityList(shared_ptr<EntityList> ents){
  if(this->destRegions_.size() > 1){
    EXCEPTION("CoefFunctionGridNodalInterp::AddEntityList called twice. Currently this is forbidden!");
  }

  if(ents->GetDefineType() == EntityList::REGION){
    this->destRegions_.insert( ents->GetRegion() );
    this->entities_.Push_back(ents);
  }else{
    EXCEPTION("CoefFunction Grid can only be called with a region.")
  }

  this->extDataInfo_->Get("DestinationRegionList")->Get("Region")->Get("name")->SetValue(ents->GetName());


}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::MapElemNodesConservative(){

  //first obtain the node coordinates of the source regions


  StdVector<Vector<Double> > nodeGlobCoords;
  StdVector<LocPoint> localCoords;
  StdVector< const Elem* > foundElements;

  Grid* destGrid = this->domain_->GetGrid();


  StdVector<UInt> allNodes;
  std::set<std::string>::iterator regIt = this->srcRegions_.begin();
  bool multipleRegions = (this->srcRegions_.size()>1);
  if(!multipleRegions){
     RegionIdType aReg = this->srcGrid_->GetRegion().Parse(*regIt);
     this->srcGrid_->GetNodesByRegion(allNodes,aReg);
  }else{
     StdVector<UInt> curNodes;
     while(regIt != this->srcRegions_.end()){
       RegionIdType aReg = this->srcGrid_->GetRegion().Parse(*regIt);
       this->srcGrid_->GetNodesByRegion(curNodes,aReg);
       allNodes.Reserve(allNodes.GetSize()+curNodes.GetSize());
       for(UInt i=0;i<curNodes.GetSize();i++){
         if(allNodes.Find(curNodes[i])==-1)
           allNodes.Push_back(curNodes[i]);
       }
       regIt++;
     }
  }
  UInt dim = this->srcGrid_->GetDim();
  nodeGlobCoords.Resize(allNodes.GetSize(),Vector<Double>(dim));
  Vector<Double> aCoord;
  for(UInt i=0;i<allNodes.GetSize();++i){
    this->srcGrid_->GetNodeCoordinate(aCoord,allNodes[i]);
    nodeGlobCoords[i] = aCoord;
  }

  destGrid->GetElemsAtGlobalCoords( nodeGlobCoords,
                                    localCoords,
                                    foundElements,
                                    this->destRegions_,
                                    false);

  //we now create some debugging information
  UInt elemCounter = 0;
  for(UInt i=0;i<foundElements.GetSize();++i){
    if(!foundElements[i]){
      elemCounter++;
    }else{
      LOG_DBG3(coeffunctiongridnodalinterp) << "Node #" << allNodes[i] << " get associated to element #" << foundElements[i]->elemNum << std::endl;
      LOG_DBG3(coeffunctiongridnodalinterp)<< "Local Coordinate is: " << localCoords[i] << std::endl << std::endl;
    }
  }
  if(elemCounter>0)
    WARN("There were " << elemCounter << " unmapped nodes. Perhaps you should think about tolerances!");

  this->extDataInfo_->Get("interpolation")->Get("conservative")->Get("numUnmappedNodes")->SetValue(elemCounter);

  //if the user wants to, we save the node->element association here
  //now we store the information in our map
  for(UInt i=0;i<foundElements.GetSize();++i){
    //check if we have the information
    if(foundElements[i]){
      std::vector<UInt> & curVec = this->elemNodeAssoc_[foundElements[i]->elemNum];
      curVec.push_back(allNodes[i]);
      this->localCoordsNodeAssoc_[allNodes[i]] = localCoords[i];
    }
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::ReadXMLNode(PtrParamNode configNode){
  //extract input and grid ids

  this->inputId_ = configNode->Get("inputId")->As<std::string>();
  this->gridId_ = configNode->Get("gridId")->As<std::string>();
  this->aSeqStep_ = this->domain_->GetDriver()->GetActSequenceStep();

  this->DetermineResult(this->inputId_,this->aSeqStep_);

  this->extDataInfo_->Get("interpolation")->Get("sourceGridID")->SetValue(this->gridId_);
  this->extDataInfo_->Get("interpolation")->Get("sourceInputID")->SetValue(this->inputId_);

  std::string interpStr = configNode->Get("interpolation")->As<std::string>();
  this->curInterpType_ = CoefFunctionGrid::InterpType_.Parse(interpStr);
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue(interpStr);
  this->extDataInfo_->Get("interpolation")->Get("setBy")->SetValue("XML");

  //obtain grid pointer
  this->srcGrid_ = this->domain_->GetGrid(this->gridId_);

  //check what kind of dependency is there...
  //hardcode to general even we have static data. perhaps we cover this case...
  this->dependType_ = CoefFunction::GENERAL;

  //TODO: Cover the case of static external data. i.e. do not bind to time/freq parameter

  //determine source region list
  StdVector<PtrParamNode> regs = configNode->Get("regionList")->GetList("region");
  for(UInt i=0;i<regs.GetSize();++i){
    std::string srcRegName = regs[i]->Get("srcRegionName")->As<std::string>();
    this->srcRegions_.insert(srcRegName);
    this->extDataInfo_->Get("SourceRegionList")->Get("Region",ParamNode::INSERT)->Get("name")->SetValue(srcRegName);
  }

  this->curStep_ = 1;
  this->curTStep_ = 0;
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::MapConservative( shared_ptr<FeSpace> targetSpace,
                                                                      Vector<DATA_TYPE>& feFncVec){

  if(!this->consInterpReady_){
    std::cout << "Preparing for conservative interpolation of external data...";
    std::cout.flush();
    boost::shared_ptr<Timer> t(new Timer);
    t->Start();
    this->MapEqns();
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);

    //ok so we need to gather the formation
    this->MapElemNodesConservative();
    CoordFormat<DATA_TYPE>* myContainer;

    //now we know which node is in which element, lets loop over the entity list and
    //create the matrix
    UInt dDim = this->resultInfo_->dofNames.GetSize();
    //determine the number of rows in the CRS matrix
    StdVector<Integer> elemEqns;
    UInt nnz=0;
    UInt numRows= targetSpace->GetNumFreeEquations();
    UInt numCols = this->solVec_.GetSize();
    std::map<UInt,std::vector<UInt> >::iterator eIt = this->elemNodeAssoc_.begin();
    while(eIt != this->elemNodeAssoc_.end()){
      const Elem* curE = this->domain_->GetGrid()->GetElem(eIt->first);
      targetSpace->GetElemEqns(elemEqns,curE);
      for(UInt i=0;i<elemEqns.GetSize();++i){
        if(elemEqns[i]>0){
          nnz += eIt->second.size();
        }
      }
      //nnz += fe->GetNumFncs()*eIt->second.size()*dDim;
      eIt++;
    }

    //ok lets create the container
    myContainer = new CoordFormat<DATA_TYPE>(numRows,numCols,nnz,false);

    shared_ptr<ElemShapeMap> esm;
    Matrix<Double> opMat;
    eIt = this->elemNodeAssoc_.begin();
    LocPointMapped lpm;
    while(eIt != this->elemNodeAssoc_.end()){
      const Elem* curE = NULL;
      curE = this->domain_->GetGrid()->GetElem(eIt->first);
      targetSpace->GetElemEqns(elemEqns,curE);
      esm = this->domain_->GetGrid()->GetElemShapeMap( curE, true );
      //now we add for each source node the corresponding entries
      for(UInt aNode=0; aNode < eIt->second.size(); aNode++){
        UInt curNodeNum = eIt->second[aNode];
        LocPoint lp = this->localCoordsNodeAssoc_[curNodeNum];
        lpm.Set(lp,esm,1.0);
        BaseFE * fe = targetSpace->GetFe(curE->elemNum);
        this->myOperator_->CalcOpMatTransposed(opMat,lpm,fe);
        for(UInt j=0;j<fe->GetNumFncs();j++){
          UInt idx = this->nodeIdxMap_[curNodeNum];
          for(UInt d =0;d<dDim;++d){
            if(elemEqns[j*dDim+d]>0){
              UInt curSrcEq =  this->eqnNumbers_[idx][d];
              myContainer->AddEntry(elemEqns[j*dDim+d]-1,curSrcEq,opMat[j*dDim+d][d]);
            }
          }
        }//add to coordMatrix
      }//for each source node
      eIt++;
    }//for each target element
    myContainer->FinaliseAssembly();
    //TODO: delete unnecessary data structures


    //now we create a CRS_Matrix from it
    this->consInterpMat_.reset(new CRS_Matrix<DATA_TYPE>(*myContainer));

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

    std::cout << "Done" << std::endl;
    std::cout.flush();
  }
  //perhaps we need to reread the solution vector from file
  if(this->GetDependency() != CoefFunction::CONSTANT){
    this->UpdateSolution();
  }

  //here it gets simple we just take the external solution vector and multiply with matrix
  this->consInterpMat_->MultAdd(this->solVec_,feFncVec);
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
   polyNode_->Get("useGridOrder", ParamNode::INSERT)->SetValue("true");
   polyNode_->Get("spectral", ParamNode::INSERT)->SetValue("true");
   polyNode_->Get("isoOrder", ParamNode::INSERT)->SetValue("1");

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
   interpolSpace = FeSpace::CreateInstance(motherNode,interpolSpaceNode,FeSpace::H1,this->domain_->GetGrid());


   shared_ptr<SolStrategy> solStrat(new SolStrategyStd(motherNode));

   interpolSpace->Init(solStrat);
   this->interpolFunction_.reset(new FeFunction<DATA_TYPE>(this->domain_->GetMathParser()));
   this->interpolFunction_->SetFeSpace(interpolSpace);
   this->interpolFunction_->SetGrid(this->domain_->GetGrid());
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
                                                                      Grid* ptGrid ){

  StdVector<LocPoint> localCoords;
  StdVector< const Elem* > foundElements;
  //build up set of source regions
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  std::set<RegionIdType> scrRegIds;
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    scrRegIds.insert(curId);
  }

  this->srcGrid_->GetElemsAtGlobalCoords( globCoord,
                                          localCoords,
                                          foundElements,
                                          scrRegIds);

  Vector<DATA_TYPE> eSol;
  Matrix<DATA_TYPE> opMat;
  LocPointMapped lpm;
  shared_ptr<ElemShapeMap> esm;
  values.Resize(globCoord.GetSize(), Vector<DATA_TYPE>(this->dimDof_));
  for(UInt i=0;i<foundElements.GetSize();++i){
    const Elem* curE = foundElements[i];
    if(!curE){
      continue;
    }
    this->GetElemSolution(eSol,curE->elemNum);
    esm = this->srcGrid_->GetElemShapeMap( curE, true );
    LocPoint lp = localCoords[i];
    lpm.Set(lp,esm,1.0);
    BaseFE * fe = esm->GetBaseFE();
    this->myOperator_->CalcOpMat(opMat,lpm,fe);
    values[i] = opMat * eSol;
  }
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalInterp<DATA_TYPE>::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                                                      StdVector< DATA_TYPE >& values,
                                                                      Grid* ptGrid ){
  StdVector< Vector<DATA_TYPE> > vecValues;
  this->GetVectorValuesAtCoords(globCoord,vecValues, ptGrid);
  values.Resize(globCoord.GetSize(),0.0);
  for(UInt i=0;i<vecValues.GetSize();++i){
    values[i] = vecValues[i][0];
  }
}

}
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionGridNodalInterp<Double>;
  template class CoefFunctionGridNodalInterp<Complex>;
#endif
