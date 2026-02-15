// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     ProgramOptions.cc
 *       \brief    <Description>
 *
 *       \date     Nov 19, 2015
 *       \author   ahueppe
 */
//================================================================================================


#include "CFSDatProgramOptions.hh"
#include "DataInOut/ColoredConsole.hh"
#include "DatUtils/Defines.hh"
#include <def_cfs_stats.hh>
#include <exception>

#include <boost/bind/bind.hpp>

namespace CFSDat{

namespace po  = boost::program_options;
namespace fs  = std::filesystem;
using namespace boost::placeholders;

CFSDatProgramOptions::CFSDatProgramOptions( CF::Integer argc, const char **argv )
                     : CF::ProgramOptions(argc,argv){

}

CFSDatProgramOptions::~CFSDatProgramOptions(){

}

void CFSDatProgramOptions::ParseData(){
    // define header std::string
    std::stringstream os;
    os << CF::fg_blue << "\n Usage: " << CF::fg_reset
       << "cfsdat [param] <name of processing run>\n\n"
       << " where [param] is one or more of the following parameters";


    // 1) Define "visible" commandLine options
    // -----------------------------------------
    po::options_description cmdVisible( os.str() );

    cmdVisible.add_options()
      ( "help,h",
        "display this usage information" )

      ( "version,v",
        "information about the CFS-Dat executable" )

      ( "paramFile,p", po::value<std::string>(),
        "name of XML parameter file for the simulation" )

      ( "numThreads,t", po::value<UInt>()->default_value(1),
        "number of threads used in CFSDat run. Default 1." )

      ( "schemaRoot,s", po::value<std::string>(),
        "path to XML schema definitions (env CFS_SCHEMA_ROOT)")

      ( "forceSegFault,f",
        "force a segmentation fault at exceptions")

      ( "noColor",
        "turn off colored output")

      ( "quiet",
        "shall reduce output, currently only for cfs command line compatibility")

      ;


    // Store help description in variable
    std::stringstream helpStream;
    helpStream << cmdVisible << std::endl;
    helpMsg_ = helpStream.str();


    // 2) Define "invisible" commandLine options
    // -----------------------------------------
    po::options_description cmdInvisible( os.str() );

    cmdInvisible.add_options()
      ("simName", po::value<std::string>(),
       "simulation name");

    // make the input file the only allowed positional parameter
    po::positional_options_description p;
    p.add( "simName", 1 );


    // 4) Combine visible and invisble options to commandLine options
    // -----------------------------------------
    po::options_description cmdLineOptions;
    cmdLineOptions.add(cmdVisible).add(cmdInvisible);

    try {

      // 5) Parse command line
      // -----------------------------------------
      po::store( po::command_line_parser( args_ ).
                 options( cmdLineOptions ).positional( p ).run(),
                 varMap_ );

      // 6) Parse environment
      // -----------------------------------------

      // first, define function pointer, which maps environment variables
      // (e.g. CFS_SCHEMA_ROOT) to local string paramter representation
      // (e.g. schemaRoot)
      boost::function1<std::string, std::string> name_mapper;

      name_mapper = boost::bind(&CFSDatProgramOptions::EnvironmentNameMapper, *this, _1);

      po::store( po::parse_environment( cmdLineOptions, name_mapper ),
                 varMap_ );

    } catch (std::exception& e) {
      EXCEPTION( "Wrong command line arguments: " << e.what() )
    }

    po::notify( varMap_ );

    // Check if colored output should be switched off
    if( varMap_.count("noColor") != 0) {
      CF::ColoredConsole::suppressed = true;
    }

    // Check for version
    if( varMap_.count("version") != 0  )
    {
      GetHeaderString(std::cout);
      PrintVersion( std::cout, true );
      exit( EXIT_SUCCESS );
    }

    // Check for help
    if( varMap_.count("help") != 0 || varMap_.count("simName") == 0)
    {
      GetHeaderString(std::cout);
      std::cout << helpMsg_;
      exit(EXIT_SUCCESS);
    }

    // obtain schema path once
    schemaPath_ = FindSchemaPath();
}

CoupledField::UInt CFSDatProgramOptions::GetNumThreads() const
{
  if( varMap_.count( "numThreads") != 0 ) {
    return varMap_["numThreads"].as<UInt>();
  }else{
    return 1;
  }
}


void CFSDatProgramOptions::GetHeaderString(std::ostream & out)
{
  if(IsQuiet())
  {
    out << ">> CFDat '" << CFS_VERSION << " " << CFS_NAME << "'"
        << " Compiled: '" << __DATE__ << "'"
        << " Build: '" << CMAKE_BUILD_TYPE << "'" << std::endl;
  }
  else
  {

    // CFS_VERSION and CFS_NAME are to be set in source/CMakeLists.txt
    out << std::endl
        << "============================================================"
        << "===========" << std::endl;
    out << " openCFS Data Processing Tool" << std::endl << std::endl
        << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
        << " (rev " << CFS_GIT_COMMIT << ")" << std::endl
        << " compiled " << __DATE__
        << " as " << CMAKE_BUILD_TYPE << std::endl
        << " CFSDat routines use " << CoupledField::CFS_NUM_THREADS << " threads for this run" << std::endl;
    out << "============================================================"
        << "==========="
        << std::endl << std::endl;
  }
}

std::string CFSDatProgramOptions::EnvironmentNameMapper(const std::string& var )  {
  std::string ret;

  if(var == "CFS_SCHEMA_ROOT") ret = "schemaRoot";
  if(var == "CFS_NO_COLOR")    ret = "noColor";
  if(var == "CFS_QUIET")       ret = "quiet";

  return ret;
}




}

