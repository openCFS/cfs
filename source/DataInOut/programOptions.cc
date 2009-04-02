// include general defines
#include <def_cfs_stats.hh>
#include <def_use_blas.hh>
#include <def_use_hdf5.hh>
#include <def_use_ansysrst.hh>
#include <def_use_mesh.hh>
#include <def_use_pardiso.hh>
#include <def_use_unv.hh>
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
#include <def_use_ansysrst.hh>

#include "programOptions.hh"
#include "coloredConsole.hh"

#include <sstream>

#include <boost/bind.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>


#ifdef USE_MKL
#include <mkl_service.h>
#endif

#ifdef USE_ACML
#include <acml.h>
#endif

#include "DataInOut/ParamHandling/InfoNode.hh"
#include "General/environment.hh"
#include "General/exception.hh"
#include "Utils/tools.hh"

// Lapack version function interface
extern "C" void ilaver_(int*, int*, int*);

namespace CoupledField {

  // extern definition
  ProgramOptions * progOpts = NULL;

  ProgramOptions::ProgramOptions( Integer argc,
                                  const char **argv )
  {

    // copy command line into vector
    args_.resize( argc-1 );
    for( Integer i = 1; i < argc; i++ ) {
      args_[i-1] = std::string( argv[i] );
    }

  }

  ProgramOptions::~ProgramOptions()
  {
  }

  void ProgramOptions::ParseData()
  {
    // define header string
    std::stringstream os;
    os << fg_blue << "\n Usage: " << fg_reset
       << "cfs [param] <name of simulation run>\n\n"
       << " where [param] is one or more of the following parameters";


    // 1) Define "visible" commandLine options
    // -----------------------------------------
    po::options_description cmdVisible( os.str() );

    cmdVisible.add_options()
      ( "help,h",
        "display this usage information" )

      ( "version,v",
        "information about the CFS++ executable" )

      ( "history,H",
        "history of revisions" )

      ( "meshFile,m", po::value<std::string>(),
        "name of mesh file for the simulation" )

      ( "paramFile,p", po::value<std::string>(),
        "name of XML parameter file for the simulation" )

      ( "schemaRoot,s", po::value<std::string>(),
        "path to XML schema definitions (env CFS_SCHEMA_ROOT)")

      ( "restart,r",
        "read restart file of previous simulation run" )

      ( "forceSegFault,f",
        "force a segmentation fault at exceptions")

      ( "printGrid,g",
        "read grid from input and write it to output file" )

      ( "writeSkeleton,w",
        "write skeleton of XML file for subsequent simulation" )

#ifdef USE_SCRIPTING
      ( "scriptFile,e", po::value<std::string>(),
        "name of script file to be evaluated" )
#endif
      ( "doProfile,d",
        "turns on generation of profiling information" )

      ( "listMapping,l",
        "add equation and local/global mapping to info.xml")

      ( "quiet,q",
        "more compressed console output (env CFS_QUIET)")

      ( "noColor",
        "turn off colored output")
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
       "simulation name")
      ;

    // make the input file the only allowed positional parameter
    po::positional_options_description p;
    p.add( "simName", 1 );


    // Explicit definition of environment options is not necessary
    // as we are using a mapping function anyway
//    // 3) Define environment options
//    // -----------------------------------------
//    po::options_description envOptions("" );
//    envOptions.add_options()
//      ("schemaRoot", po::value<std::string>(),
//       "Path to schema definitions of CFS++ installation")
//      ;

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

      name_mapper = boost::bind(&ProgramOptions::EnvironmentNameMapper, *this, _1);

      po::store( po::parse_environment( cmdLineOptions, name_mapper ),
                 varMap_ );

    } catch (std::exception& e) {
      EXCEPTION( "Wrong command line arguments: " << e.what() );
    }

    po::notify( varMap_ );

    // Check if colored output should be switched off
    if( varMap_.count("noColor") != 0) {
      ColoredConsole::suppressed = true;
    }

    // Check for version
    if( varMap_.count("version") != 0  ) 
    {
      GetHeaderString(std::cout);      
      GetVersionString( std::cout, true );
      exit( EXIT_SUCCESS );
    }

    // Check for history
    if( varMap_.count("history") != 0  ) 
    {
      GetHeaderString(std::cout);      
      GetHistoryString(std::cout);
      exit( EXIT_SUCCESS );
    }

    // Check for help
    if( varMap_.count("help") != 0) 
    {
      GetHeaderString(std::cout);      
      std::cout << helpMsg_;
      exit(EXIT_SUCCESS);
    }

    // If no argument was given, print additional information
    if( varMap_.count("simName") == 0 )
    {
      GetHeaderString(std::cout);
      std::cout << "cfs: no input files. Pleas run with --help for help\n";
      exit(EXIT_SUCCESS);
    }
  }

  std::string ProgramOptions::EnvironmentNameMapper(const std::string& var )  {
    std::string ret;

    if(var == "CFS_SCHEMA_ROOT") ret = "schemaRoot";
    if(var == "CFS_QUIET")       ret = "quiet";
    
    return ret;
  }

  std::string ProgramOptions::GetSimName() const
  {
    if( varMap_.count( "simName") != 0 ) {

      // get complete path
      fs::path simPath = fs::path( varMap_["simName"].as<std::string>(),
                                  fs::native );

      // return only file name without path information
      return simPath.leaf();
    } else {
      return "";
    }
  }

  fs::path ProgramOptions::GetSimPath() const
   {
     if( varMap_.count( "simName") != 0 ) {

       // get complete path
       fs::path simPath ( varMap_["simName"].as<std::string>(),
                          fs::native );

       fs::complete( simPath.branch_path()).native_directory_string();

       // return path to simulation
       return fs::complete( simPath.branch_path());
     } else {
       return fs::initial_path().native_file_string();
     }
   }

  std::string ProgramOptions::GetSimPathStr() const
   {
     // get complete path
     fs::path simPath = GetSimPath();

     return simPath.native_directory_string();
   }

  fs::path ProgramOptions::GetParamFile() const
  {
    if( varMap_.count( "paramFile" ) == 0 ) {
      return GetSimPath() / fs::path(GetSimName()+".xml" );
    } else {
      fs::path paramPath( varMap_["paramFile"].as<std::string>(),
                          fs::native);
      return fs::system_complete( paramPath );
    }
  }

  std::string ProgramOptions::GetParamFileStr() const
  {
    fs::path paramPath = GetParamFile();

    return paramPath.native_file_string();
  }

#ifdef USE_SCRIPTING
  fs::path ProgramOptions::GetScriptFile() const
  {
    if( varMap_.count("scriptFile") > 0 ) {
      fs::path scriptPath( varMap_["scriptFile"].as<std::string>() );
      return fs::system_complete( scriptPath);
    } else {
      return fs::path();
    }
  }

  std::string ProgramOptions::GetScriptFileStr() const
  {
    fs::path scriptPath = GetScriptFile();

    return scriptPath.native_file_string();
  }
#endif

  fs::path ProgramOptions::GetSchemaPath() const
  {
    std::string schema;
    fs::path schemaPath;

    // If the user specified a path on the command line use it instead.
    if( varMap_.count( "schemaRoot" ) ) {
      schema = varMap_[ "schemaRoot" ].as<std::string>();
    } else {
      schema = XMLSCHEMA;
    }

    // If none of the above paths exists, raise exception
    if(!fs::exists(schema)) {
      EXCEPTION( "Schema path '" << schema << "' does not exist!" );
    }

    try
    {
      schemaPath = fs::system_complete(schema);
      schemaPath.normalize();
    } catch (fs::filesystem_error& ex)
    {
      EXCEPTION(ex.what());
    }

    return schemaPath;
  }

  std::string ProgramOptions::GetSchemaPathStr() const
  {
    fs::path schemaPath = GetSchemaPath();

    return schemaPath.native_directory_string();
  }

  fs::path ProgramOptions::GetMeshFile() const
  {

    if( varMap_.count( "meshFile") != 0 ) {
      fs::path meshPath( varMap_["meshFile"].as<std::string>(),
                         fs::native );
      return fs::system_complete( meshPath );
    } else {
      return fs::path();
    }
  }

  std::string ProgramOptions::GetMeshFileStr() const
  {
    fs::path meshFile = GetMeshFile();

    return meshFile.native_file_string();
  }

  bool ProgramOptions::GetPrintGrid() const
  {
    return (varMap_.count( "printGrid") > 0);
  }

  bool ProgramOptions::GetDoProfile() const
  {
    return (varMap_.count( "doProfile") > 0);
  }


  bool ProgramOptions::GetRestart() const
  {
    return (varMap_.count( "restart") > 0);
  }


  bool ProgramOptions::GetWriteSkeleton() const
  {
    return (varMap_.count( "writeSkeleton") > 0);
  }

  bool ProgramOptions::GetForceSegFault() const
  {
    return (varMap_.count( "forceSegFault") > 0);
  }

  bool ProgramOptions::DoListMapping() const
  {
    return varMap_.count("listMapping") > 0;
  }

  bool ProgramOptions::IsQuiet() const
  {
    return varMap_.count("quiet") > 0;
  }
  
  void ProgramOptions::PrintHelp( std::ostream& out )
  {
    out << helpMsg_;
  }

  void ProgramOptions::ToInfo(InfoNode* in) const
  {
    in->SetComment("values of command line parameters (including defaults)");
    in->Get("problem", "name of simulation run")->SetValue(GetSimName());
    in->Get("parameterFile")->SetValue(GetParamFileStr());
    in->Get("schemaPath")->SetValue(GetSchemaPathStr());
    in->Get("meshFile")->SetValue(GetMeshFileStr());
#ifdef USE_SCRIPTING
    in->Get("scriptFile")->SetValue(GetScriptFileStr());
#endif
    in->Get("printGrid")->SetValue(GetPrintGrid());
    in->Get("doProfile")->SetValue(GetDoProfile());
    in->Get("restart")->SetValue(GetRestart());
    in->Get("writeSkeleton")->SetValue(GetWriteSkeleton());
    in->Get("listMapping")->SetValue(DoListMapping());
    in->Get("forceSegFault")->SetValue(GetForceSegFault());
  }

  
  void ProgramOptions::GetVersionString( std::ostream & outstr,
                                         bool colorise)
  {
    bool colTmp = ColoredConsole::colorise;
    ColoredConsole::colorise = colorise;

    std::string build_type = CMAKE_BUILD_TYPE;

    outstr << "CFS_VERSION:           "
           << fg_blue << CFS_VERSION << fg_reset << std::endl

           << "CFS_NAME:              "
           << fg_blue << CFS_NAME << fg_reset << std::endl

           << "CFS_BUILD_HOST:        "
           << fg_blue << CFS_BUILD_HOST << fg_reset << std::endl

           << "CFS_BUILD_USER:        "
           << fg_blue << CFS_BUILD_USER << fg_reset << std::endl

           << "CFS_SUBVERSION_REV:    "
           << fg_blue << CFS_SUBVERSION_REV << fg_reset << std::endl

           << "CFS_SUBVERSION_REPOS:  "
           << fg_blue << CFS_SUBVERSION_REPOS
           << fg_reset << std::endl << std::endl

           << "CFS_CXX_COMPILER_NAME: "
           << fg_blue << CFS_CXX_COMPILER_NAME << fg_reset << std::endl

           << "CFS_CXX_COMPILER_VER:  "
           << fg_blue << CFS_CXX_COMPILER_VER << fg_reset
           << std::endl << std::endl

           << "CFS_FORTRAN_COMPILER_NAME: "
           << fg_blue << CFS_FORTRAN_COMPILER_NAME << fg_reset<< std::endl

           << "CFS_FORTRAN_COMPILER_VER:  "
           << fg_blue << CFS_FORTRAN_COMPILER_VER
           << fg_reset<< std::endl << std::endl

           << "CFS_DISTRO:            "
           << fg_blue << CFS_DISTRO << fg_reset << std::endl

           << "CFS_DISTRO_VER:        "
           << fg_blue << CFS_DISTRO_VER << fg_reset << std::endl

           << "CFS_ARCH:              "
           << fg_blue << CFS_ARCH
           << fg_reset << std::endl << std::endl

           << "CMAKE_BUILD_TYPE:      "
           << fg_blue  << build_type << fg_reset << std::endl;
#ifndef NDEBUG
      outstr << "COMPILE_FLAGS:       "
             << fg_blue  << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_DEBUG << fg_reset << std::endl;

      outstr << "LINK_FLAGS:          "
             << fg_blue  << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_DEBUG
             << fg_reset
             << std::endl << std::endl;
#else
      outstr << "COMPILE_FLAGS:       "
             << fg_blue  << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_RELEASE << fg_reset << std::endl;

      outstr << "LINK_FLAGS:          "
             << fg_blue << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_RELEASE
             << fg_reset << std::endl << std::endl;
#endif
 #ifdef USE_ARPACK
    outstr << "USE_ARPACK:            "
           << fg_blue  << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_ARPACK:            "
           << fg_blue  << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_BLAS
    outstr << "USE_BLAS:              "
           << fg_blue  << "YES" << fg_reset << std::endl;
    outstr << "BLAS_IMPLEMENTATION:   "
           << fg_blue  << CFS_BLAS_LAPACK << fg_reset << std::endl;
 #ifdef USE_MKL
    MKLVersion ver;

    MKLGetVersion(&ver);

    outstr << "MKL_VERSION:           " << fg_blue
           << ver.MajorVersion << "."
           << ver.MinorVersion << "."
           << ver.BuildNumber << fg_reset
           << std::endl;
    outstr << "MKL_PRODSTAT:          " << fg_blue
           << ver.ProductStatus << fg_reset
           << std::endl;
    outstr << "MKL_BUILD:             " << fg_blue
           << ver.Build << fg_reset
           << std::endl;

    MKL_FreeBuffers();
 #endif
 #ifdef USE_ACML
    Integer acml_major, acml_minor, acml_patch;

    acmlversion(&acml_major, &acml_minor, &acml_patch);

    outstr << "ACML_VERSION:          " << fg_blue
           << acml_major << "."
           << acml_minor << "."
           << acml_patch << fg_reset
           << std::endl;
    acmlinfo();
 #endif
#else
    outstr << "USE_BLAS:              "
           << fg_blue  << "NO" << fg_reset << std::endl;
#endif
    outstr << std::endl;

 #ifdef USE_LAPACK
    outstr << "USE_LAPACK:            "
           << fg_blue << "YES" << fg_reset << std::endl;
    Integer major, minor, rev;
    ilaver_(&major, &minor, &rev);
    outstr << "LAPACK_VERSION:        "
           << fg_blue << major << "." << minor << "." << rev
           << fg_reset << std::endl;
#else
    outstr << "USE_LAPACK:            "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_ILUPACK
    outstr << std::endl;
    outstr << "USE_ILUPACK:           "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_ILUPACK:           "
           << fg_blue  << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_PARDISO
    outstr << std::endl;
    outstr << "USE_PARDISO:           "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "PARDISO_IMPL:          "
           << fg_blue  << CFS_PARDISO << fg_reset << std::endl;
 #else
    outstr << "USE_PARDISO:           "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_METIS
    outstr << std::endl
           << "USE_METIS:             "
           << fg_blue << "YES" << fg_reset << std::endl;
#else
    outstr << "USE_METIS:             "
           << fg_blue << "NO" << fg_reset<< std::endl;
#endif

#ifdef MpCCI
    outstr << "MpCCI:                 "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "MpCCI_RELEASE:         "
           << fg_blue << MpCCI_RELEASE
           << fg_reset << std::endl;
#else
    outstr << "MpCCI:                 "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

    outstr << std::endl;

 #ifdef USE_SCRIPTING_TCL
    outstr << "USE_SCRIPTING_TCL:     "
           << fg_blue  << "YES" << fg_reset << std::endl;
    outstr << "CFS_TCL_VERSION:       "
           << fg_blue << CFS_TCL_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_SCRIPTING_TCL:     "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_SCRIPTING_PYTHON
    outstr << "USE_SCRIPTING_PYTHON:  "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_PYTHON_VERSION:    "
           << fg_blue << CFS_PYTHON_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_SCRIPTING_PYTHON:  "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

    outstr << std::endl;

 #ifdef USE_GIDPOST
    outstr << "USE_GIDPOST:           "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_GIDPOST_VERSION:   "
           << fg_blue << CFS_GIDPOST_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_GIDPOST:           "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_GMV_INPUT
    outstr << "USE_GMV_INPUT:         "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_GMV_INPUT:         "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_GMV_OUTPUT
    outstr << "USE_GMV_OUTPUT:        "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_GMV_OUTPUT:        "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_HDF5
    outstr << "USE_HDF5:              "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_HDF5_VERSION:      "
           << fg_blue << CFS_HDF5_VERSION << fg_reset << std::endl;
    outstr << "CFS_ZLIB_VERSION:      "
           << fg_blue << CFS_ZLIB_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_HDF5:              "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_MESH
    outstr << "USE_MESH:              "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_MESH:              "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_UNV
    outstr << "USE_UNV:               "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_UNV:               "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

    outstr << std::endl;
    outstr << "CFS_BOOST_VERSION:     "
           << fg_blue << CFS_BOOST_VERSION << fg_reset << std::endl;
    outstr << std::endl;
    outstr << "CFS_METIS_VERSION:     "
           << fg_blue << CFS_METIS_VERSION << fg_reset << std::endl;
    outstr << std::endl;


#ifdef USE_XERCES
    outstr << "USE_XERCES:            "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_XERCES_VERSION:    "
           << fg_blue << CFS_XERCES_VERSION << fg_reset << std::endl;
    outstr << "XMLSCHEMA:             "
           << fg_blue << XMLSCHEMA << fg_reset << std::endl;
    outstr << "USE_XERCES:            "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

    outstr << "USE_ANSYSRST:          "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_ANSYS_VERSION:    "
           << fg_blue << CFS_ANSYS_VERSION << fg_reset << std::endl;

  ColoredConsole::colorise = colTmp;
  }

  void ProgramOptions::GetHistoryString( std::ostream & out)
  {
    out << "This is a incomplete revision history. It starts with the new versioning schema:" << std::endl
        << " * CFS-Version: <yy>.<mm> with two digits for year and month" << std::endl
        << " * CFS-Name:    Changing name on every major release" << std::endl
        << std::endl
        << "08.12, (pre) Kerniger Kaerntner" << std::endl
        << "  Sneak preview of the 'Kerniger Kaerntner' release which will bring the" << std::endl
        << "  brand new OLAS from the MACHETE branch. But here only first steps of" << std::endl
        << "  the info.xml file." << std::endl
        << std::endl
        << "09.03, Maehende Machete" << std::endl
        << "  This finally is the famous MACHETE branch with a 0-based OLAS, unified" << std::endl
        << "  matrix and vector classes and the brand new CFSDEPS building matching libs." << std::endl
        << std::endl;
  }

  void ProgramOptions::GetHeaderString(std::ostream & out)
  {  
    if(IsQuiet())
    {
      out << ">> CFS++ '" << CFS_VERSION << " " << CFS_NAME << "'"
          << " Compiled: '" << __DATE__ << "'"
          << " Build: '" << CMAKE_BUILD_TYPE << "'" << std::endl;
    }
    else
    {
      // CFS_VERSION and CFS_NAME are to be set in source/CMakeLists.txt
      out << std::endl
          << "============================================================"
          << "===========" << std::endl;
      out << " CFS++ - Coupled Field Simulation" << std::endl << std::endl
          << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
          << " (rev " << CFS_SUBVERSION_REV << ")" << std::endl
          << " compiled " << __DATE__
          << " as " << CMAKE_BUILD_TYPE << std::endl;
      out << "============================================================"
          << "==========="
          << std::endl << std::endl;
    }
  }
}
