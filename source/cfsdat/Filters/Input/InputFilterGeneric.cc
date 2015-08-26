// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     GenericInputFilter.cc
 *       \brief    <Description>
 *
 *       \date     Aug 13, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "cfsdat/Filters/Input/InputFilterGeneric.hh"
#include "DataInOut/DefineInOutFiles.hh"

namespace CFSDat{

InputFilterGeneric::InputFilterGeneric(){

  //Utilize the CFS way to create input files.
  CoupledField::DefineInOutFiles * createInput = new CoupledField::DefineInOutFiles();

  //Starting from here, we would (of course) depend on user input
  PtrParamNode configNode;
  PtrParamNode infoNode;
  configNode.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));
  infoNode.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));

  PtrParamNode inputNode = configNode->Get("fileFormats",ParamNode::APPEND);
  PtrParamNode hdfNode = inputNode->Get("input",ParamNode::INSERT)->Get("hdf5",ParamNode::INSERT);
  hdfNode->Get("fileName",ParamNode::INSERT)->SetValue("test.h5");
  hdfNode->Get("id",ParamNode::INSERT)->SetValue("default");
  hdfNode->Get("gridId",ParamNode::INSERT)->SetValue("default");

  //create input files (dynamic IDs or, as usual, from XML)
  inFile_ = CoupledField::DefineInOutFiles::CreateSingleInputFileObject("test.h5",hdfNode,infoNode);

  //delete the input again
  delete createInput;

}

InputFilterGeneric::~InputFilterGeneric(){

}

void InputFilterGeneric::Init(){
  inFile_->InitModule();
}

void InputFilterGeneric::Run(){

}

}
