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

//boost parameter handling and command line
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>

//CFS includes
#include <def_cfs_stats.hh>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"


#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "General/defs.hh"
#include "General/Environment.hh"

using namespace CoupledField;
namespace CoupledField
{
  class SimInput;
  class SimOutput;
}
//! Global parameter node and info node instance
PtrParamNode param;
PtrParamNode info;

namespace CFSDat {

  namespace po = boost::program_options;

  void GenerateParamNode(int argc, char** argv, PtrParamNode param){
    //lets test with some options

    po::options_description desc("Supported Modes");
    desc.add_options()
            ("help,h", "Produce help message")
            ("version,v", "Print version information");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    po::notify(vm);

    //1. check for config files
    //   - try to validate. But how?
    //2.
  }

}

int main(int argc, char** argv)
{
  std::string infoFileName;
  SetEnvironmentEnums();
  BasePDE::SetEnums();
  EntityList::SetEnums();
  ElemShape::Initialize();

  domain = NULL;

  try
   {
     param.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));
     CFSDat::GenerateParamNode(argc, argv, param);

   } catch(std::exception& ex) {
     std::cerr << "The following error occured during cfsdat execution:\n\n" << ex.what();
     if (info != NULL)
     {
       info->Get(ParamNode::PN_ERROR)->SetValue(ex.what());
       info->ToFile(infoFileName);
     }
     return -1;
   }
   //info->ToFile(infoFileName);

}
