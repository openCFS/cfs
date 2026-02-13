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

//boost parameter handling and command line
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <filesystem>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#define  BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>

//CFS includes
#include <def_cfs_stats.hh>
#include <def_use_openmp.hh>
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "PDE/BasePDE.hh"
#include "Utils/Timer.hh"
#include <sstream>

#include <boost/version.hpp>
#include <boost/asio/ip/host_name.hpp>



//CFSDatIncludes
#include "cfsdat/Filters/BaseFilter.hh"
#include "cfsdat/DatUtils/CFSDatProgramOptions.hh"
#include "cfsdat/DatUtils/ResultManager.hh"
#include "cfsdat/DatUtils/DataStructs.hh"
#include "cfsdat/DatUtils/Defines.hh"

#ifdef USE_OPENMP
#include <omp.h>
#endif


namespace CFSDat {

  namespace po = boost::program_options;

  void GenerateParamNode(int argc, char** argv, PtrParamNode param){

  }

}



int main(int argc, const char** argv)
{

  CFSDat::CFSDatProgramOptions * options = new CFSDat::CFSDatProgramOptions(argc,argv);

  options->ParseData();

  SetEnvironmentEnums();


#ifdef USE_OPENMP
  //now set the number of threads from the commandline
  SetNumberOfThreads(options->GetNumThreads());
  ////SetNumberOfThreads(omp_get_max_threads());
#else
  SetNumberOfThreads(1);
#endif
  BasePDE::SetEnums();
  EntityList::SetEnums();
  ElemShape::Initialize();

  // Log program startup
  options->GetHeaderString( std::cout );

  // this is our hostname. Empty if it cannot be determined */
  std::string hostname = boost::asio::ip::host_name();

  // This is a string for output with the start time */
  std::string start_time;

  // Print information about program start time and host
  std::cout << "Simulation run started at " << Timer::TimeStamp() << std::endl;
  if(!hostname.empty()) std::cout<< "on " << hostname << std::endl;





  // this is the new param stuff which replaces the old params - delete this comment finally
  std::string schema = options->GetSchemaPathStr();
  schema += "/CFS-Dat/CFS_Dat.xsd";

  //start the timer
  CoupledField::shared_ptr<CoupledField::Timer> datTimer(new CoupledField::Timer);
  datTimer->Start();

  // Write information to command line
  std::cout << "--- Reading parameter file " << std::endl;
  PtrParamNode configNode = XmlReader::ParseFile(options->GetParamFileStr(), schema,
                                                 "http://www.cfs++.org/simulation");

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
        std::cout << "\t---> Adding filter \"" << allFilters[*inIter]->GetId() << "\" as source for filter \"" <<  fIter->second->GetId() << "\"" << std::endl;
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
    std::cout << "\t---> Processing output filters for step #" << counter++ << std::endl;
    allFinished = true;
    outIter =  outputs.begin();
    for(;outIter != outputs.end();++outIter){
      std::cout << "\t     Filter ID: " << (*outIter)->GetId() << std::endl;
      allFinished &= (*outIter)->Run();
    }
  } while(!allFinished);

  datTimer->Stop();
  std::stringstream elapsed;
  const int walltime((int) datTimer->GetWallTime());


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


  std::cout << std::endl << "---> COMPUTATION DONE. Time elapsed: " << elapsed.str() << std::endl << std::endl;

  delete options;

  return 0;
}
