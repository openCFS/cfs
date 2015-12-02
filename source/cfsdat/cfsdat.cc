// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     cfsdat.cc
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================

//System includes
#include <iostream>
#include <stdio.h>


//CFS includes
#include <def_cfs_stats.hh>
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "PDE/BasePDE.hh"

//CFSDatIncludes
#include "cfsdat/Filters/BaseFilter.hh"
#include "cfsdat/Utils/CFSDatProgramOptions.hh"
#include "cfsdat/Utils/ResultManager.hh"
#include "cfsdat/Utils/DataStructs.hh"
#include "cfsdat/Utils/Defines.hh"



namespace CFSDat {

  namespace po = boost::program_options;

  void GenerateParamNode(int argc, char** argv, PtrParamNode param){

  }

}



int main(int argc, const char** argv)
{



  CFSDat::CFSDatProgramOptions * options = new CFSDat::CFSDatProgramOptions(argc,argv);
  // Write information to command line
  std::cout << "--- Reading parameter file " << std::endl;

  options->ParseData();

  SetEnvironmentEnums();
  BasePDE::SetEnums();
  EntityList::SetEnums();
  ElemShape::Initialize();


  // this is the new param stuff which replaces the old params - delete this comment finally
  std::string schema = options->GetSchemaPathStr();
  schema += "/CFS-Dat/CFS_Dat.xsd";

  // Initialize our xerces dom parser to handle the cfs xml file
  Xerces xerces(schema);

  xerces.SetFile(options->GetParamFileStr());

  CoupledField::PtrParamNode configNode = xerces.CreateParamNodeInstance();

  CFSDat::PtrResultManager resMan(new CFSDat::ResultManager());

  PtrParamNode pipelineNode = configNode->Get("pipeline");

  //loop over all elements, create filters
  ParamNodeList filters = pipelineNode->GetChildren();

  //create filters
  std::cout << "--- Creating Filters" << std::endl;
  std::map<std::string,CFSDat::FilterPtr> allFilters;
  std::set<CFSDat::FilterPtr> outputs;
  for(UInt i=0;i<filters.GetSize();++i){
    CFSDat::FilterPtr newFilt = CFSDat::BaseFilter::Generate(filters[i],resMan);
    if(newFilt){
      std::cout << "\t---> Adding Filter type \"" << filters[i]->GetName() << "\" with ID \"" << newFilt->GetId() << "\"" << std::endl;
      allFilters[newFilt->GetId()] = newFilt;
      if(newFilt->IsOutput()){
        outputs.insert(newFilt);
      }
    }
  }

  std::cout << "--- Creating filter connections" << std::endl;
  //create connections
  std::map<std::string,CFSDat::FilterPtr>::iterator fIter = allFilters.begin();
  for(;fIter != allFilters.end();++fIter){
    std::set<std::string>::iterator inIter = fIter->second->GetInputIds().begin();
    for(;inIter != fIter->second->GetInputIds().end();++inIter){
      if(*inIter != "none"){
        std::cout << "\t ---> Adding filter \"" << allFilters[*inIter]->GetId() << "\" as source for filter \"" <<  fIter->second->GetId() << "\"" << std::endl;
        CFSDat::BaseFilter::Connect(allFilters[*inIter],fIter->second);
      }
    }
  }

  std::cout << "--- Initializing filters" << std::endl;
  //for now we just loop over output filters
  std::set<CFSDat::FilterPtr>::iterator outIter =  outputs.begin();
  for(;outIter != outputs.end();++outIter){
    (*outIter)->InitResults();
  }

  std::cout << "--- Finalizing filters" << std::endl;
  resMan->Finalize();

  outIter =  outputs.begin();
  for(;outIter != outputs.end();++outIter){
    (*outIter)->FinishInit();
  }

  //after that no invalid result definition should be remaining
  std::set<boost::uuids::uuid> list = resMan->GetActiveResults();
  if(list.size() != 0){
    CoupledField::Exception("there are still Active results! This may not happen!");
  }


  UInt counter = 1;
  bool allFinished = true;
  std::cout << "--- Initiating filter traversal" << std::endl;
  do{
    std::cout << "\t---> Traversing outputs: step #" << counter++ << std::endl;
    allFinished = true;
    outIter =  outputs.begin();
    for(;outIter != outputs.end();++outIter){
      allFinished &= (*outIter)->Run();
    }
  } while(!allFinished);

  std::cout << "--- Computation done." << std::endl;
  delete options;

  return 0;
}
