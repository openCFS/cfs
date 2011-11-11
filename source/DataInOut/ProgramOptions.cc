// include general defines
#include <def_cfs_stats.hh>
#include <def_cplreader.hh>
#include <def_use_blas.hh>
#include <def_use_hdf5.hh>
#include <def_use_ansysrst.hh>
#include <def_use_mesh.hh>
#include <def_use_pardiso.hh>
#include <def_use_unv.hh>
#include <def_use_gidpost.hh>
#include <def_use_ilupack.hh>
#include <def_use_cholmod.hh>
#include <def_use_metis.hh>
#include <def_use_python.hh>
#include <def_use_xerces.hh>
#include <def_use_arpack.hh>
#include <def_use_gmv.hh>
#include <def_use_gmsh.hh>
#include <def_use_lapack.hh>
#include <def_use_mpcci.hh>
#include <def_use_interpolation.hh>
#include <def_use_tcl.hh>
#include <def_xmlschema.hh>

#include "ProgramOptions.hh"
#include "ColoredConsole.hh"

#include <sstream>
#include <cstdlib>

#include <boost/bind.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>

#ifdef USE_MKL
#include <mkl_service.h>
#endif

#ifdef USE_ACML
#include <acml.h>
#endif

#include <bzlib.h>
#include <zlib.h>

#ifdef USE_ILUPACK
#include <amd.h>
#endif

#ifdef USE_CHOLMOD
#include <cholmod.h> 
#endif

#ifdef USE_SCRIPTING_PYTHON
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#include <Python.h>
#endif

#ifdef USE_SCRIPTING_TCL
#include <tcl.h>
#endif

#ifdef USE_ARPACK
#include <arpack_version.h>
#endif

#ifdef USE_INTERPOLATION
#include <CGAL/version.h>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Utils/tools.hh"

#define _QUOTEME(x) #x
#define QUOTEME(x) _QUOTEME(x)


using std::string;
using std::endl;
using std::cout;



// Lapack version function interface
extern "C" void ilaver_(int*, int*, int*);

// CFX IO library version interface
extern "C" const char* io_get_version();


namespace CoupledField {

  // extern definition
  ProgramOptions * progOpts = NULL;

  ProgramOptions::ProgramOptions( Integer argc,
                                  const char **argv )
  {

    // copy command line into vector
    args_.resize( argc-1 );
    for( Integer i = 1; i < argc; i++ ) {
      args_[i-1] = string( argv[i] );
    }

    exe_ = argv[0];
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

      ( "meshFile,m", po::value<string>(),
        "name of mesh file for the simulation" )

      ( "paramFile,p", po::value<string>(),
        "name of XML parameter file for the simulation" )

      ( "schemaRoot,s", po::value<string>(),
        "path to XML schema definitions (env CFS_SCHEMA_ROOT)")

      ( "ersatz,x", po::value<string>(),
        "name of ersatz material density file")

      ( "restart,r",
        "read restart file of previous simulation run" )

      ( "forceSegFault,f",
        "force a segmentation fault at exceptions")

      ( "printGrid,g",
        "read grid from input and write it to output file" )

      ( "exportGrid,G",
        "export the grid to the info.xml file, best with -g")

      ( "writeSkeleton,w",
        "write skeleton of XML file for subsequent simulation" )

#ifdef USE_SCRIPTING
      ( "scriptFile,e", po::value<string>(),
        "name of script file to be evaluated" )
#endif
      ( "mapping,l",
        "equation and constraints to info.xml")

      ( "detailed,d",
        "detailed output to info.xml")

      ( "quiet,q",
        "more compressed console output (env CFS_QUIET)")

      ( "noColor",
        "turn off colored output")
      ;


    // Store help description in variable
    std::stringstream helpStream;
    helpStream << cmdVisible << endl;
    helpMsg_ = helpStream.str();


    // 2) Define "invisible" commandLine options
    // -----------------------------------------
    po::options_description cmdInvisible( os.str() );

    cmdInvisible.add_options()
      ("simName", po::value<string>(),
       "simulation name");

    // make the input file the only allowed positional parameter
    po::positional_options_description p;
    p.add( "simName", 1 );


    // Explicit definition of environment options is not necessary
    // as we are using a mapping function anyway
//    // 3) Define environment options
//    // -----------------------------------------
//    po::options_description envOptions("" );
//    envOptions.add_options()
//      ("schemaRoot", po::value<string>(),
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
      boost::function1<string, string> name_mapper;

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
      GetHeaderString(cout);      
      GetVersionString( cout, true );
      exit( EXIT_SUCCESS );
    }

    // Check for history
    if( varMap_.count("history") != 0  ) 
    {
      GetHeaderString(cout);      
      GetHistoryString(cout);
      exit( EXIT_SUCCESS );
    }

    // Check for help
    if( varMap_.count("help") != 0 || varMap_.count("simName") == 0) 
    {
      GetHeaderString(cout);      
      cout << helpMsg_;
      exit(EXIT_SUCCESS);
    }

    /*
    // If no argument was given, print additional information
    if( varMap_.count("simName") == 0 )
    {
      GetHeaderString(cout);
      cout << "cfs: no input files. Please run with --help for help\n";
      exit(EXIT_SUCCESS);
    }
    */
  }

  string ProgramOptions::EnvironmentNameMapper(const string& var )  {
    string ret;

    if(var == "CFS_SCHEMA_ROOT") ret = "schemaRoot";
    if(var == "CFS_NO_COLOR")    ret = "noColor";
    if(var == "CFS_QUIET")       ret = "quiet";
    
    return ret;
  }

  string ProgramOptions::GetSimName() const
  {
    if( varMap_.count( "simName") != 0 ) {

      // get complete path
      fs::path simPath = fs::path( varMap_["simName"].as<string>(),
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
       fs::path simPath ( varMap_["simName"].as<string>(),
                          fs::native );

       fs::complete( simPath.branch_path()).native_directory_string();

       // return path to simulation
       return fs::complete( simPath.branch_path());
     } else {
       return fs::initial_path().native_file_string();
     }
   }

  string ProgramOptions::GetSimPathStr() const
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
      fs::path paramPath( varMap_["paramFile"].as<string>(),
                          fs::native);
      return fs::system_complete( paramPath );
    }
  }

  string ProgramOptions::GetParamFileStr() const
  {
    fs::path paramPath = GetParamFile();

    return paramPath.native_file_string();
  }

#ifdef USE_SCRIPTING
  fs::path ProgramOptions::GetScriptFile() const
  {
    if( varMap_.count("scriptFile") > 0 ) {
      fs::path scriptPath( varMap_["scriptFile"].as<string>() );
      return fs::system_complete( scriptPath);
    } else {
      return fs::path();
    }
  }

  string ProgramOptions::GetScriptFileStr() const
  {
    fs::path scriptPath = GetScriptFile();

    return scriptPath.native_file_string();
  }
#endif

  string ProgramOptions::GetErsatzMaterialStr() const
  {
    if(varMap_.count("ersatz") != 0)
      return varMap_["ersatz"].as<string>();
    else
      return "";
  }

  fs::path ProgramOptions::GetSchemaPath() const
  {
    string schema;
    fs::path schemaPath;

    // If the user specified a path on the command line use it instead.
    if( varMap_.count( "schemaRoot" ) ) {
      schema = varMap_[ "schemaRoot" ].as<string>();
    } else {
      schema = XMLSCHEMA;
    }

    //      fs::path fn = fs::system_complete(progOpts->exe_);
    //      fn.normalize();

    try
    {
      schemaPath = fs::system_complete(schema);
      schemaPath.normalize();
    } catch (fs::filesystem_error& ex)
    {
      EXCEPTION(ex.what());
    }

    // If none of the above paths exists, raise exception
    if(!fs::exists(schemaPath)) {
      EXCEPTION( "Schema path '" << schemaPath << "' does not exist!" );
    }

    return schemaPath;
  }

  string ProgramOptions::GetSchemaPathStr() const
  {
    fs::path schemaPath = GetSchemaPath();

    return schemaPath.native_directory_string();
  }

  fs::path ProgramOptions::GetMeshFile() const
  {

    if( varMap_.count( "meshFile") != 0 ) {
      fs::path meshPath( varMap_["meshFile"].as<string>(),
                         fs::native );
      return fs::system_complete( meshPath );
    } else {
      return fs::path();
    }
  }

  string ProgramOptions::GetMeshFileStr() const
  {
    fs::path meshFile = GetMeshFile();

    return meshFile.native_file_string();
  }

  bool ProgramOptions::GetPrintGrid() const
  {
    return (varMap_.count( "printGrid") > 0);
  }

  bool ProgramOptions::DoExportGrid() const
  {
    return varMap_.count("exportGrid") > 0;
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
    return varMap_.count("mapping") > 0;
  }

  bool ProgramOptions::DoDetailedInfo() const
  {
    return varMap_.count("detailed") > 0;
  }


  bool ProgramOptions::IsQuiet() const
  {
    return varMap_.count("quiet") > 0;
  }
  
  void ProgramOptions::PrintHelp( std::ostream& out )
  {
    out << helpMsg_;
  }

  void ProgramOptions::ToInfo(PtrParamNode in) const
  {
    in->SetComment("values of command line parameters (including defaults)");
    in->Get("problem")->SetValue(GetSimName());
    in->Get("parameterFile")->SetValue(GetParamFileStr());
    in->Get("schemaPath")->SetValue(GetSchemaPathStr());
    in->Get("meshFile")->SetValue(GetMeshFileStr());
#ifdef USE_SCRIPTING
    in->Get("scriptFile")->SetValue(GetScriptFileStr());
#endif
    in->Get("mapping")->SetValue(DoListMapping());
    in->Get("detailed")->SetValue(DoDetailedInfo());
    in->Get("MKL_NUM_THREADS")->SetValue(getenv("MKL_NUM_THREADS") != NULL ? getenv("MKL_NUM_THREADS") : "-");

    // cfs information
    in = in->Get("cfs");
    in->Get("version")->SetValue(CFS_VERSION);
    in->Get("name")->SetValue(CFS_NAME);
    in->Get("build")->SetValue(CMAKE_BUILD_TYPE);
    in->Get("svn_revision")->SetValue(CFS_SUBVERSION_REV);
  }

  
  void ProgramOptions::GetVersionString( std::ostream & out,
                                         bool colorise)
  {
    bool colTmp = ColoredConsole::colorise;
    ColoredConsole::colorise = colorise;
    
    string build_type = CMAKE_BUILD_TYPE;
    string cxxFlags = CMAKE_CXX_FLAGS;
    string ldFlags = CMAKE_EXE_LINKER_FLAGS;
    boost::trim(cxxFlags);
    boost::trim(ldFlags);    
    
    out << "CFS_VERSION:           "
        << fg_blue << CFS_VERSION << fg_reset << endl

        << "CFS_NAME:              "
        << fg_blue << CFS_NAME << fg_reset << endl;

    if(progOpts) 
    {      
      fs::path fn = fs::system_complete(progOpts->exe_);
      fn.normalize();
      out << "CFS_EXECUTABLE:        "
          << fg_blue << fn.native_directory_string() << fg_reset << endl;
    }
    
    out << "CFS_BUILD_HOST:        "
        << fg_blue << CFS_BUILD_HOST << fg_reset << endl

        << "CFS_BUILD_USER:        "
        << fg_blue << CFS_BUILD_USER << fg_reset << endl

        << "CFS_SUBVERSION_REV:    "
        << fg_blue << CFS_SUBVERSION_REV << fg_reset << endl

        << "CFS_SUBVERSION_REPOS:  "
        << fg_blue << CFS_SUBVERSION_REPOS
        << fg_reset << endl << endl

        << "CFS_CXX_COMPILER_NAME: "
        << fg_blue << CFS_CXX_COMPILER_NAME << fg_reset << endl

        << "CFS_CXX_COMPILER_VER:  "
        << fg_blue << CFS_CXX_COMPILER_VER << fg_reset
        << endl << endl

        << "CFS_FORTRAN_COMPILER_NAME: "
        << fg_blue << CFS_FORTRAN_COMPILER_NAME << fg_reset<< endl

        << "CFS_FORTRAN_COMPILER_VER:  "
        << fg_blue << CFS_FORTRAN_COMPILER_VER
        << fg_reset<< endl << endl

        << "CFS_DISTRO:            "
        << fg_blue << CFS_DISTRO << fg_reset << endl

        << "CFS_DISTRO_VER:        "
        << fg_blue << CFS_DISTRO_VER << fg_reset << endl

        << "CFS_ARCH:              "
        << fg_blue << CFS_ARCH
        << fg_reset << endl << endl

        << "CMAKE_BUILD_TYPE:      "
        << fg_blue  << build_type << fg_reset << endl;
#ifndef NDEBUG
    out << "COMPILE_FLAGS:         "
        << fg_blue << cxxFlags << " " << CMAKE_CXX_FLAGS_DEBUG
        << fg_reset << endl;

    out << "LINK_FLAGS:            "
        << fg_blue  << CMAKE_EXE_LINKER_FLAGS
        << " " << CMAKE_EXE_LINKER_FLAGS_DEBUG
        << fg_reset
        << endl << endl;
#else
    out << "COMPILE_FLAGS:         "
        << fg_blue << cxxFlags
        << " " << CMAKE_CXX_FLAGS_RELEASE << fg_reset << endl;
    
    out << "LINK_FLAGS:            "
        << fg_blue << ldFlags
        << " " << CMAKE_EXE_LINKER_FLAGS_RELEASE
        << fg_reset << endl << endl;
#endif
 #ifdef USE_ARPACK    
    out << "USE_ARPACK:            "
        << fg_blue  << "YES" << fg_reset << endl;
    out << "ARPACK_VERSION:        "
        << fg_blue  << ARPACK_VERSION_NUMBER << " "
        << ARPACK_VERSION_DATE << fg_reset << endl;
 #else
    out << "USE_ARPACK:            "
        << fg_blue  << "NO" << fg_reset << endl;
#endif

 #ifdef USE_BLAS
    out << "USE_BLAS:              "
        << fg_blue  << "YES" << fg_reset << endl;
    out << "BLAS_IMPLEMENTATION:   "
        << fg_blue  << CFS_BLAS_LAPACK << fg_reset << endl;
 #ifdef USE_MKL
    CFSMKLVersion ver;

    MKLGetVersion(reinterpret_cast<MKLVersion*>(&ver));

    out << "MKL_VERSION:           " << fg_blue
        << ver.MajorVersion << "."
        << ver.MinorVersion << "."
        << ver.BuildNumber 
        << fg_reset
        << endl;
    out << "MKL_PRODSTAT:          " << fg_blue
        << ver.ProductStatus << fg_reset
        << endl;
    out << "MKL_BUILD:             " << fg_blue
        << ver.Build << fg_reset
        << endl;

    MKL_FreeBuffers();

    out << "MKL_NUM_THREADS:       "
        << fg_blue << (getenv("MKL_NUM_THREADS") != NULL ? getenv("MKL_NUM_THREADS") : "-")
        << fg_reset << endl;
 #endif
 #ifdef USE_ACML
    Integer acml_major, acml_minor, acml_patch;

    acmlversion(&acml_major, &acml_minor, &acml_patch);

    out << "ACML_VERSION:          " << fg_blue
        << acml_major << "."
        << acml_minor << "."
        << acml_patch << fg_reset
        << endl;
    acmlinfo();
 #endif
#else
    out << "USE_BLAS:              "
        << fg_blue  << "NO" << fg_reset << endl;
#endif
    out << endl;

 #ifdef USE_LAPACK
    out << "USE_LAPACK:            "
        << fg_blue << "YES" << fg_reset << endl;
    Integer major, minor, rev;
    ilaver_(&major, &minor, &rev);
    out << "LAPACK_VERSION:        "
        << fg_blue << major << "." << minor << "." << rev
        << fg_reset << endl;
#else
    out << "USE_LAPACK:            "
        << fg_blue << "NO" << fg_reset << endl;
#endif

 #ifdef USE_ILUPACK
    out << endl;
    out << "USE_ILUPACK:           "
        << fg_blue << "YES" << fg_reset << endl;
    out << "AMD_VERSION:           "
        << fg_blue << AMD_MAIN_VERSION << "." << AMD_SUB_VERSION << "."
        << AMD_SUBSUB_VERSION << " " << AMD_DATE
        << fg_reset << endl;
 #else
    out << "USE_ILUPACK:           "
        << fg_blue  << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_CHOLMOD
    out << endl;
    out << "USE_CHOLMOD:           "
        << fg_blue << "YES" << fg_reset << endl;
    out << "SUITESPARSE_VERSION:   "
        << fg_blue << SUITESPARSE_MAIN_VERSION << "." << SUITESPARSE_SUB_VERSION << "."
        << SUITESPARSE_SUBSUB_VERSION << " " << SUITESPARSE_DATE
        << fg_reset << endl;
 #else
    out << "USE_CHOLMOD:           "
        << fg_blue  << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_PARDISO
    out << endl;
    out << "USE_PARDISO:           "
        << fg_blue << "YES" << fg_reset << endl;
    out << "PARDISO_IMPL:          "
        << fg_blue  << CFS_PARDISO << fg_reset << endl;
 #else
    out << "USE_PARDISO:           "
           << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_METIS
    out << endl
        << "USE_METIS:             "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_METIS_VERSION:     "
        << fg_blue << CFS_METIS_VERSION
        << fg_reset << endl;
#else
    out << "USE_METIS:             "
        << fg_blue << "NO" << fg_reset<< endl;
#endif

    out << endl;

 #ifdef USE_SCRIPTING_TCL
    out << "USE_SCRIPTING_TCL:     "
        << fg_blue  << "YES" << fg_reset << endl;
    Integer tcl_major, tcl_minor, tcl_rev, tcl_patch;
    Tcl_GetVersion(&tcl_major, &tcl_minor, &tcl_rev, &tcl_patch);
    out << "CFS_TCL_VERSION:       "
        << fg_blue << tcl_major << "." << tcl_minor << "." << tcl_rev
        << fg_reset << endl;
 #else
    out << "USE_SCRIPTING_TCL:     "
        << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_SCRIPTING_PYTHON
    out << "USE_SCRIPTING_PYTHON:  "
        << fg_blue << "YES" << fg_reset << endl;
    string pyVersion = Py_GetVersion();
    UInt i=0;
    for(i=0; i < pyVersion.length(); i++) {
      if(pyVersion[i] == ' ')
        break;
    }
    pyVersion = pyVersion.substr(0, i);
    out << "CFS_PYTHON_VERSION:    "
        << fg_blue << pyVersion << fg_reset << endl;
 #else
    out << "USE_SCRIPTING_PYTHON:  "
        << fg_blue << "NO" << fg_reset << endl;
#endif

    out << endl;

#ifdef USE_ANSYSRST
    out << "USE_ANSYSRST:          "
        << fg_blue << "YES" << fg_reset << endl;
#else
    out << "USE_ANSYSRST:          "
        << fg_blue << "NO" << fg_reset << endl;
#endif

#ifdef USE_GIDPOST
    out << "USE_GIDPOST:           "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_GIDPOST_VERSION:   "
        << fg_blue << CFS_GIDPOST_VERSION << fg_reset << endl;
 #else
    out << "USE_GIDPOST:           "
        << fg_blue << "NO" << fg_reset << endl;
 #endif

#ifdef USE_GMSH
    out << "USE_GMSH:              "
        << fg_blue << "YES" << fg_reset << endl;
#else
    out << "USE_GMSH:              "
        << fg_blue << "NO" << fg_reset << endl;
#endif

#ifdef USE_GMV_INPUT
    out << "USE_GMV_INPUT:         "
        << fg_blue << "YES" << fg_reset << endl;
 #else
    out << "USE_GMV_INPUT:         "
        << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_GMV_OUTPUT
    out << "USE_GMV_OUTPUT:        "
        << fg_blue << "YES" << fg_reset << endl;
 #else
    out << "USE_GMV_OUTPUT:        "
        << fg_blue << "NO" << fg_reset << endl;
#endif

 #ifdef USE_HDF5
    out << "USE_HDF5:              "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_HDF5_VERSION:      "
        << fg_blue << CFS_HDF5_VERSION << fg_reset << endl;
 #else
    out << "USE_HDF5:              "
        << fg_blue << "NO" << fg_reset << endl;
#endif

 #ifdef USE_MESH
    out << "USE_MESH:              "
        << fg_blue << "YES" << fg_reset << endl;
 #else
    out << "USE_MESH:              "
        << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_UNV
    out << "USE_UNV:               "
        << fg_blue << "YES" << fg_reset << endl;
 #else
    out << "USE_UNV:               "
        << fg_blue << "NO" << fg_reset << endl;
#endif

    out << endl;

#ifdef USE_XERCES
    out << "USE_XERCES:            "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_XERCES_VERSION:    "
        << fg_blue << CFS_XERCES_VERSION << fg_reset << endl;
    out << "XMLSCHEMA:             ";
    if(progOpts)
      out << fg_blue << progOpts->GetSchemaPath() << fg_reset << endl;
    else
      out << fg_blue << XMLSCHEMA << fg_reset << endl;
    
#else
    out << "USE_XERCES:            "
        << fg_blue << "NO" << fg_reset << endl;
#endif

    out << endl;

#ifdef MpCCI
    out << "MpCCI:                 "
        << fg_blue << "YES" << fg_reset << endl;
    out << "MpCCI_RELEASE:         "
        << fg_blue << MpCCI_RELEASE
        << fg_reset << endl;
#else
    out << "MpCCI:                 "
        << fg_blue << "NO" << fg_reset << endl;
#endif
#ifdef USE_INTERPOLATION
    out << "USE_INTERPOLATION:     "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_CGAL_VERSION:      "
        << fg_blue << QUOTEME(CGAL_VERSION)
        << " (" << CGAL_VERSION_NR
        << ", SVN rev. " << CGAL_SVN_REVISION << ")"
        << fg_reset << endl;
#else
    out << "USE_INTERPOLATION:     "
        << fg_blue << "NO" << fg_reset << endl;
#endif
    
    out << endl;
    out << "CFS_BOOST_VERSION:     "
        << fg_blue << CFS_BOOST_VERSION << fg_reset << endl;
    out << "CFS_ZLIB_VERSION:      "
        << fg_blue << zlibVersion() << fg_reset << endl;
    out << "CFS_BZIP2_VERSION:     "
        << fg_blue << BZ2_bzlibVersion() << fg_reset << endl;
    out << "CFS_MUPARSER_VERSION:  "
        << fg_blue << CFS_MUPARSER_VERSION << fg_reset << endl;

   
#ifdef CPLREADER_CFX
    out << "CFX_IO_VERSION:        "
        << fg_blue << io_get_version() << fg_reset << endl;
#endif

    out << endl;
    out << "sizeof(int)            "  << fg_blue << sizeof(int) << fg_reset << endl;
    out << "sizeof(Integer)        "  << fg_blue << sizeof(Integer) << fg_reset << endl;
    out << "sizeof(size_t)         "  << fg_blue << sizeof(size_t) << fg_reset << endl;
    out << "sizeof(long)           "  << fg_blue << sizeof(long) << fg_reset << endl;    
    out << "sizeof(long int)       "  << fg_blue << sizeof(long int) << fg_reset << endl;
    out << "sizeof(double)         "  << fg_blue << sizeof(double) << fg_reset << endl;    
    out << "sizeof(void*)          "  << fg_blue << sizeof(void*) << fg_reset << endl;
    out << "CLOCKS_PER_SEC         "  << fg_blue << CLOCKS_PER_SEC << fg_reset << endl;
    
  ColoredConsole::colorise = colTmp;
  }

  void ProgramOptions::GetHistoryString( std::ostream & out)
  {
    out << "This is a incomplete revision history. It starts with the new versioning schema:" << endl
        << " * CFS-Version: <yy>.<mm> with two digits for year and month" << endl
        << " * CFS-Name:    Changing name on every major release" << endl
        << endl
        << "08.12, (pre) Kerniger Kaerntner" << endl
        << "  Sneak preview of the 'Kerniger Kaerntner' release which will bring the" << endl
        << "  brand new OLAS from the MACHETE branch. But here only first steps of" << endl
        << "  the info.xml file." << endl
        << endl
        << "09.03, Maehende Machete" << endl
        << "  This finally is the famous MACHETE branch with a 0-based OLAS, unified" << endl
        << "  matrix and vector classes and the brand new CFSDEPS building libs." << endl
        << endl
        << "09.08, Rasanter Rechner" << endl
        << "  ILUPACK by M. Bollhoefer works again (iterative preconditioner/ solver package)" << endl
        << "  CHOLMOD (as in MATLAB) is a s.p.d. direct solver, >= 20% faster than PARDISO" << endl
        << endl
        << "10.07, Optimal Orgasm" << endl
        << "  This shall reflect that CFS++ has become also a strong structural optimization" << endl
        << "  tool - yet the acceptance by engineers could still be optimized (e.g. by more" << endl
        << "  realistic optimization problems)." << endl
        << endl
        << "10.10, Prickly Porcupine" << endl
        << "  No major advance beside taking care that young students are not confronted with" << endl
        << "  X-rated terms. Luckily they don't hear developers conversation with compilers!." << endl;

  }

  void ProgramOptions::GetHeaderString(std::ostream & out)
  {  
    if(IsQuiet())
    {
      out << ">> CFS++ '" << CFS_VERSION << " " << CFS_NAME << "'"
          << " Compiled: '" << __DATE__ << "'"
          << " Build: '" << CMAKE_BUILD_TYPE << "'" << endl;
    }
    else
    {
      // CFS_VERSION and CFS_NAME are to be set in source/CMakeLists.txt
      out << endl
          << "============================================================"
          << "===========" << endl;
      out << " CFS++ - Coupled Field Simulation" << endl << endl
          << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'"
          << " (rev " << CFS_SUBVERSION_REV << ")" << endl
          << " compiled " << __DATE__
          << " as " << CMAKE_BUILD_TYPE << endl;
      out << "============================================================"
          << "==========="
          << endl << endl;
    }
  }
}
