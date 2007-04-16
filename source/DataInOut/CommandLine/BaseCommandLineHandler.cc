// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_cfs_stats.hh>
#include <def_use_blas.hh>
#include <def_use_hdf5.hh>
#include <def_use_mesh.hh>
#include <def_use_ansysrst.hh>
#include <def_use_pardiso.hh>
#include <def_use_unv.hh>
#include <def_cfs_stats.hh>
#include <def_use_gidpost.hh>
#include <def_use_ilupack.hh>
#include <def_use_metis.hh>
#include <def_use_python.hh>
#include <def_use_xerces.hh>
#include <def_use_arpack.hh>
#include <def_use_gmv.hh>
#include <def_use_lapack.hh>
#include <def_use_mpcci.hh>
#include <def_use_tcl.hh>
#include <def_xmlschema.hh>

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

// Be a little bit "colorful"
#ifndef SUPPRESS_COLORED_OUTPUT
// #define COLOR_INIT "\033[1m"
#define COLOR_INIT "\033[34m"
#define COLOR_STOP "\033[0m"
#else
#define COLOR_INIT ""
#define COLOR_STOP ""
#endif

namespace CoupledField {

  // ------------------------------------------------
  //  Definition of static data members of the class
  // ------------------------------------------------

  // Help strings
  const std::string BaseCommandLineHandler::helpSimName_       =
  "name of simulation run";

  const std::string BaseCommandLineHandler::helpParamFile_     =
  "name of XML parameter file for the simulation";

  const std::string BaseCommandLineHandler::helpSchemaPath_    =
  "path to directory containing the XML schema file";

  const std::string BaseCommandLineHandler::helpMeshFile_      =
  "name of mesh file for the simulation";  
#ifdef USE_SCRIPTING
  const std::string BaseCommandLineHandler::helpScriptFileName_ =
  "name of script file to be evaluated";  
#endif

  const std::string BaseCommandLineHandler::helpTraceDepth_    =
#ifdef TRACE
  "depth of function tracing";
#else
  "depth of function tracing \033[1m(de-activated in this version)\033[0m";
#endif

  const std::string BaseCommandLineHandler::helpWriteSkeleton_ =
  "write skeleton of XML file for subsequent simulation";

  const std::string BaseCommandLineHandler::helpPrintGrid_     =
  "only read grid from mesh-file and write it to output file";

  const std::string BaseCommandLineHandler::helpShowEqnMap_    =
  "print node to equation number mapping to OLAS report file";

  const std::string BaseCommandLineHandler::helpDoProfile_     =
  "turns on generation of profiling information";

  const std::string BaseCommandLineHandler::helpRestart_     =
  "read restart file";

  const std::string BaseCommandLineHandler::helpHelp_          =
  "print this usage information";

  const std::string BaseCommandLineHandler::helpDumpStats_          =
  "dump information about the CFS++ executable";

#ifdef DEBUG
  const std::string BaseCommandLineHandler::helpForceSegFault_      =
    "force a segmentation fault instead of throwing an exception";
#endif

  // Short forms of markers
  const std::string BaseCommandLineHandler::markerParamFile_       = "-p";
  const std::string BaseCommandLineHandler::markerMeshFile_        = "-m";
#ifdef USE_SCRIPTING
  const std::string BaseCommandLineHandler::markerScriptFileName_  = "-e";
#endif
  const std::string BaseCommandLineHandler::markerTraceDepth_      = "-t";
  const std::string BaseCommandLineHandler::markerWriteSkeleton_   = "-w";
  const std::string BaseCommandLineHandler::markerPrintGrid_       = "-g";
  const std::string BaseCommandLineHandler::markerHelp_            = "-h";
  const std::string BaseCommandLineHandler::markerSchemaPath_      = "-s";
  const std::string BaseCommandLineHandler::markerDoProfile_       = "-d";
  const std::string BaseCommandLineHandler::markerRestart_         = "-r";
  const std::string BaseCommandLineHandler::markerDumpStats_       = "-u";
#ifdef DEBUG
  const std::string BaseCommandLineHandler::markerForceSegFault_   = "-f";
#endif

  // Long forms of markers
  const std::string BaseCommandLineHandler::markerLongParamFile_     =
  "--paramFile";
  const std::string BaseCommandLineHandler::markerLongMeshFile_      =
  "--meshFile";
#ifdef USE_SCRIPTING
  const std::string BaseCommandLineHandler::markerLongScriptFileName_ =
  "--scriptFile";
#endif
  const std::string BaseCommandLineHandler::markerLongTraceDepth_    =
  "--traceDepth";
  const std::string BaseCommandLineHandler::markerLongWriteSkeleton_ =
  "--writeSkeleton";
  const std::string BaseCommandLineHandler::markerLongPrintGrid_     =
  "--printGrid";
  const std::string BaseCommandLineHandler::markerLongHelp_          =
  "--help";
  const std::string BaseCommandLineHandler::markerLongSchemaPath_    =
  "--schemaPath";
  const std::string BaseCommandLineHandler::markerLongShowEqnMap_    =
  "--showEqnMap";
  const std::string BaseCommandLineHandler::markerLongDoProfile_     =
  "--doProfile";
  const std::string BaseCommandLineHandler::markerLongRestart_     =
  "--restart";
  const std::string BaseCommandLineHandler::markerLongDumpStats_     =
  "--dumpStats";
#ifdef DEBUG
  const std::string BaseCommandLineHandler::markerLongForceSegFault_ =
  "--forceSegFault";
#endif


  // ---------------------------------------------------------------
  //  Implementation of some class methods too large for the header
  // ---------------------------------------------------------------


  // **************
  //   PrintUsage
  // **************
  void BaseCommandLineHandler::PrintUsage() {

    ENTER_FCN( "BaseCommandLineHandler::PrintUsage" );

    std::ostream &os = std::cout;

    os << COLOR_INIT << "\n Usage: " << COLOR_STOP

       << "./cfs [param] <name of simulation run>\n\n"
       << " where [param] is one or more of the following parameters\n\n"

      // --paramFile
       << COLOR_INIT
       << " " << markerParamFile_ << ", " << markerLongParamFile_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpParamFile_ << "\n\n"

      // --schemaPath
       << COLOR_INIT
       << " " << markerSchemaPath_ << ", " << markerLongSchemaPath_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpSchemaPath_ << "\n\n"

      // --meshFile
       << COLOR_INIT
       << " " << markerMeshFile_ << ", " << markerLongMeshFile_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpMeshFile_ << "\n\n"

#ifdef USE_SCRIPTING
      // --scriptFile
       << COLOR_INIT
       << " " << markerScriptFileName_ << ", " << markerLongScriptFileName_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpScriptFileName_ << "\n\n"
#endif

      // --traceDepth
       << " " << COLOR_INIT
       << markerTraceDepth_ << ", " << markerLongTraceDepth_
       << COLOR_STOP
       << " = <non-negative integer>\n"
       << " " << helpTraceDepth_ << "\n\n"

      // --printGrid
       << " " << COLOR_INIT
       << markerPrintGrid_ << ", " << markerLongPrintGrid_
       << COLOR_STOP
       << '\n'
       << " " << helpPrintGrid_ << "\n\n"

      // --WriteSkeleton
       << " " << COLOR_INIT
       << markerWriteSkeleton_ << ", " << markerLongWriteSkeleton_
       << COLOR_STOP
       << '\n'
       << " " << helpWriteSkeleton_ << "\n\n"

      // --showEqnMap
       << " " << COLOR_INIT
       << " "  << ", " << markerLongShowEqnMap_
       << COLOR_STOP
       << '\n'
       << " " << helpShowEqnMap_ << "\n\n"

      // --doProfile
       << " " << COLOR_INIT
       << markerDoProfile_ << ", " << markerLongDoProfile_
       << COLOR_STOP
       << '\n'
       << " " << helpDoProfile_ << "\n\n"

      // --Restart
       << " " << COLOR_INIT
       << markerRestart_ << ", " << markerLongRestart_
       << COLOR_STOP
       << '\n'
       << " " << helpRestart_ << "\n\n"

      // --help
       << " " << COLOR_INIT
       << markerHelp_ << ", " << markerLongHelp_
       << COLOR_STOP
       << '\n'
       << " " << helpHelp_ << "\n\n"

#ifdef DEBUG

      // --forceSegFault
       << " " << COLOR_INIT
       << markerForceSegFault_ << ", " << markerLongForceSegFault_
       << COLOR_STOP
       << '\n'
       << " " << helpForceSegFault_ << "\n\n"
#endif
      // --dumpStats
       << " " << COLOR_INIT
       << markerDumpStats_ << ", " << markerLongDumpStats_
       << COLOR_STOP
       << '\n'
       << " " << helpDumpStats_ << "\n\n";


  }


  // ***************
  //   PrintParams
  // ***************
  void BaseCommandLineHandler::PrintParams( std::ostream &out,
                                            bool colorise ) {

    ENTER_FCN( "BaseCommandLineHandler::PrintParams" );

    std::string colorInit = "";
    std::string colorStop = "";

    if ( colorise == true ) {
      colorInit = COLOR_INIT;
      colorStop = COLOR_STOP;
    }

    out << colorInit
        << "\n\nValues of command line parameters (including defaults):\n\n"
        << colorStop

        << " name of simulation run = "
        << colorInit
        << GetSimName()
        << colorStop << '\n'

        << ' ' << markerLongParamFile_ << " = "
        << colorInit
        << GetParamFile()
        << colorStop << '\n'

        << ' ' << markerLongSchemaPath_ << " = "
        << colorInit
        << GetSchemaPath()
        << colorStop << '\n'

        << ' ' << markerLongMeshFile_ << " = "
        << colorInit
        << GetMeshFile()
        << colorStop << '\n'

#ifdef USE_SCRIPTING
        << ' ' << markerLongScriptFileName_ << " = "
        << colorInit
        << GetScriptFileName()
        << colorStop << '\n'
#endif

        << ' ' << markerLongTraceDepth_ << " = "
        << colorInit
        << GetTraceDepth()
        << colorStop << '\n'

        << ' ' << markerLongPrintGrid_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetPrintGrid() == true )
        << colorStop << '\n'

        << ' ' << markerLongShowEqnMap_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetShowEqnMap() == true )
        << colorStop << '\n'

        << ' ' << markerLongDoProfile_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetDoProfile() == true )
        << colorStop << '\n'

        << ' ' << markerLongRestart_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetRestart() == true )
        << colorStop << '\n'

        << ' ' << markerLongWriteSkeleton_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetWriteSkeleton() == true )
        << colorStop << "\n"

#ifdef DEBUG
        << ' ' << markerLongForceSegFault_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetForceSegFault() == true )
        << colorStop << "\n"
#endif

        << ' ' << markerLongDumpStats_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetDumpStats() == true )
        << colorStop << "\n\n";




  }

  void BaseCommandLineHandler::GetDumpString(std::string& dumpStr)
  {
    std::stringstream outstr;
    
    std::string build_type = CMAKE_BUILD_TYPE;

    outstr << "CFS_VERSION:           " << CFS_VERSION << std::endl
           << "CFS_BUILD_HOST:        " << CFS_BUILD_HOST << std::endl
           << "CFS_CONF_DATE:         " << CFS_CONF_DATE << std::endl << std::endl
           << "CFS_SUBVERSION_REV:    " << CFS_SUBVERSION_REV << std::endl
           << "CFS_SUBVERSION_REPOS:  " << CFS_SUBVERSION_REPOS << std::endl << std::endl
           << "CFS_CXX_COMPILER_NAME: " << CFS_CXX_COMPILER_NAME << std::endl
           << "CFS_CXX_COMPILER_VER:  " << CFS_CXX_COMPILER_VER << std::endl << std::endl
           << "CFS_FORTRAN_COMPILER_NAME: " << CFS_FORTRAN_COMPILER_NAME << std::endl
           << "CFS_FORTRAN_COMPILER_VER:  " << CFS_FORTRAN_COMPILER_VER << std::endl << std::endl
           << "CFS_DISTRO:            " << CFS_DISTRO << std::endl
           << "CFS_DISTRO_VER:        " << CFS_DISTRO_VER << std::endl
           << "CFS_ARCH:              " << CFS_ARCH << std::endl << std::endl
           << "CMAKE_BUILD_TYPE:      " << build_type << std::endl;
    if(build_type == "DEBUG") 
    {
      outstr << "COMPILE_FLAGS:       " << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_DEBUG << std::endl;
      outstr << "LINK_FLAGS:          " << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_DEBUG
             << std::endl << std::endl;
    }
    else
    {
      outstr << "COMPILE_FLAGS:       " << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_RELEASE << std::endl;
      outstr << "LINK_FLAGS:          " << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_RELEASE
             << std::endl << std::endl;
    }

 #ifdef USE_ARPACK
    outstr << "USE_ARPACK:            YES" << std::endl;
 #else
    outstr << "USE_ARPACK:            NO" << std::endl;
 #endif

 #ifdef USE_BLAS
    outstr << "USE_BLAS:              YES" << std::endl;
 #else
    outstr << "USE_BLAS:              NO" << std::endl;
 #endif

 #ifdef USE_LAPACK
    outstr << "USE_LAPACK:            YES" << std::endl;
 #else
    outstr << "USE_LAPACK:            NO" << std::endl;
 #endif

 #ifdef USE_ILUPACK
    outstr << "USE_ILUPACK:           YES" << std::endl;
 #else
    outstr << "USE_ILUPACK:           NO" << std::endl;
 #endif

 #ifdef USE_PARDISO
    outstr << "USE_PARDISO:           YES" << std::endl;
 #else
    outstr << "USE_PARDISO:           NO" << std::endl;
 #endif

 #ifdef USE_METIS
    outstr << "USE_METIS:             YES" << std::endl;
 #else
    outstr << "USE_METIS:             NO" << std::endl;
 #endif

 #ifdef MpCCI
    outstr << "MpCCI:                 YES" << std::endl;
    outstr << "MpCCI_RELEASE:         " << MpCCI_RELEASE
           << std::endl;
 #else
    outstr << "MpCCI:                 NO" << std::endl;
 #endif

    outstr << std::endl;

 #ifdef USE_SCRIPTING_TCL
    outstr << "USE_SCRIPTING_TCL:     YES" << std::endl;
 #else
    outstr << "USE_SCRIPTING_TCL:     NO" << std::endl;
 #endif

 #ifdef USE_SCRIPTING_PYTHON
    outstr << "USE_SCRIPTING_PYTHON:  YES" << std::endl;
 #else
    outstr << "USE_SCRIPTING_PYTHON:  NO" << std::endl;
 #endif

    outstr << std::endl;

 #ifdef USE_GIDPOST
    outstr << "USE_GIDPOST:           YES" << std::endl;
 #else
    outstr << "USE_GIDPOST:           NO" << std::endl;
 #endif

 #ifdef USE_GMV_INPUT
    outstr << "USE_GMV_INPUT:         YES" << std::endl;
 #else
    outstr << "USE_GMV_INPUT:         NO" << std::endl;
 #endif

 #ifdef USE_GMV_OUTPUT
    outstr << "USE_GMV_OUTPUT:        YES" << std::endl;
 #else
    outstr << "USE_GMV_OUTPUT:        NO" << std::endl;
 #endif

 #ifdef USE_HDF5
    outstr << "USE_HDF5:              YES" << std::endl;
 #else
    outstr << "USE_HDF5:              NO" << std::endl;
 #endif

 #ifdef USE_MESH
    outstr << "USE_MESH:              YES" << std::endl;
 #else
    outstr << "USE_MESH:              NO" << std::endl;
 #endif

 #ifdef USE_UNV
    outstr << "USE_UNV:               YES" << std::endl;
 #else
    outstr << "USE_UNV:               NO" << std::endl;
 #endif

 #ifdef USE_ANSYSRST
    outstr << "USE_ANSYSRST:          YES" << std::endl;
    outstr << "ANSYS_VERSION:         " << ANSYS_VERSION << std::endl;
 #else
    outstr << "USE_ANSYSRST:          NO" << std::endl;
 #endif

    outstr << std::endl;

 #ifdef USE_XERCES
    outstr << "USE_XERCES:            YES" << std::endl;
    outstr << "XERCES_VERSION:        " << XERCES_VERSION << std::endl;
    outstr << "XMLSCHEMA:             " << XMLSCHEMA << std::endl;
 #else
    outstr << "USE_XERCES:            NO" << std::endl;
 #endif
    
    dumpStr = outstr.str();
  }
  

  void BaseCommandLineHandler::DumpStats() {
    std::string dumpStr;
    
    GetDumpString(dumpStr);
    
    std::cout << dumpStr;
  }
    
}
