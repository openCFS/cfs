// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MeshBasedDerivative.cc
 *       \brief    <Description>
 *
 *       \date     Oct 4, 2016
 *       \author   kroppert
 */
//================================================================================================

#include "MeshBasedDerivative.hh"
#include "DataInOut/DefineInOutFiles.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

namespace CFSDat{

MeshBasedDerivative::MeshBasedDerivative(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                      :BaseDerivativeFilter(numWorkers,config,resMan){

  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  std::string outRes = params_->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outRes);
  upResNames.insert(inRes);


  ParamNodeList srcList =  params_->Get("regions")->Get("sourceRegions")->GetList("region");
  for(UInt iP = 0; iP < srcList.GetSize(); ++iP){
    srcRegions_.insert(srcList[iP]->Get("name")->As<std::string>());
  }

  ParamNodeList trgList =  params_->Get("regions")->Get("targetRegions")->GetList("region");
  for(UInt iP = 0; iP < trgList.GetSize(); ++iP){
    trgRegions_.insert(trgList[iP]->Get("name")->As<std::string>());
  }

  //Now we read the input grid. Another possibility would be to give the user the opportunity
  //to perform interpolation onto a grid which already has results
  //a grid only input filter is not directly possible due to the concept of a result driven pipeline.

  //create grid
  CreateDummyCfsParamNode();
  PtrParamNode infoNode;
  std::string filename = params_->Get("targetMesh")->GetChild()->Get("fileName")->As<std::string>();
  trgMeshInp_ = CoupledField::DefineInOutFiles::CreateSingleInputFileObject(filename,filterId_,params_->Get("targetMesh")->GetChild(),infoNode);
  trgMeshInp_->InitModule();
  UInt maxDim = trgMeshInp_->GetDim();

  trgGrid_ = new CF::GridCFS(maxDim,dummyXMLNode_,infoNode,filterId_);


  trgMeshInp_->ReadMesh(trgGrid_);
  //it would be nice not to finish the grid here
  //in order to let other filters add some entities
  //unfortunately this is not possible as we can not access anything
  //without it another question, how can two inputs share a common grid?
  trgGrid_->FinishInit();

}

void MeshBasedDerivative::FinishInit(){
  //first go up
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->FinishInit();
  }
  PrepareDifferentiation();
}

//TODO for now copy paste code from inputfilter... not very nice
void MeshBasedDerivative::CreateDummyCfsParamNode(){

  PtrParamNode meshInputNode = params_->Get("targetMesh")->GetChild();
  dummyXMLNode_.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));
  CoupledField::PtrParamNode iNode = dummyXMLNode_->Get("fileFormats",ParamNode::APPEND)->Get("input",ParamNode::INSERT);
  iNode->AddChildNode(meshInputNode);

  //create domain node
  //UInt dim = inFile_->GetDim();
  CoupledField::PtrParamNode dNode = dummyXMLNode_->Get("domain",ParamNode::APPEND);
  //dNode->Get("geometryType",ParamNode::APPEND)->SetValue(itoa(dim));

}

}


