// include general defines
#include <def_cfs_stats.hh>
#include <def_use_blas.hh>
#include <def_use_hdf5.hh>
#include <def_use_ansysrst.hh>
#include <def_use_mesh.hh>
#include <def_use_pardiso.hh>
#include <def_use_unv.hh>
#include <def_use_comsol.hh>
#include <def_use_gidpost.hh>
#include <def_use_ilupack.hh>
#include <def_use_suitesparse.hh>
#include <def_use_lis.hh>
#include <def_use_metis.hh>
#include <def_use_xerces.hh>
#include <def_use_arpack.hh>
#include <def_use_gmv.hh>
#include <def_use_gmsh.hh>
#include <def_use_cgns.hh>
#include <def_use_ccmio.hh>
#include <def_use_lapack.hh>
#include <def_use_cgal.hh>
#include <def_use_libfbi.hh>
#include <def_use_flann.hh>
#include <def_xmlschema.hh>
#include <def_use_openmp.hh>
#include <def_disable_optimization.hh>


#include <def_cfs_fortran_interface.hh>

#include "ProgramOptions.hh"
#include "ColoredConsole.hh"

#include <sstream>
#include <cstdlib>

#include <boost/bind.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/trim.hpp>

#ifdef USE_MKL
#include <mkl_service.h>
#ifndef mkl_get_version
#define MKL_Get_Version MKLGetVersion
#define MKL_Free_Buffers MKL_FreeBuffers
#endif
#endif

#ifdef USE_ACML
#include <acml.h>
#endif
#include <bzlib.h>
#include <zlib.h>

#include <H5public.h>

#ifdef USE_METIS
#include <defs.h> 
#endif

#ifdef USE_SUITESPARSE
#include <cholmod.h> 
#include <umfpack.h> 
#include <amd.h> 
#endif

#ifdef USE_LIS
#include <lis.h> 
#endif

#ifdef USE_ARPACK
#include <arpack_version.h>
#endif

#ifdef USE_CGAL
#include <CGAL/version.h>
#endif

#ifdef USE_FLANN
#include <flann/flann.hpp>
#endif

#include <xercesc/util/XercesVersion.hpp>

#ifdef USE_GIDPOST
#include <gidpost.h>
#endif

#ifdef USE_CGNS
#include <cgnslib.h>
// The NO_ERROR and NO_DATA symbols are defined in some windows headers
// and conflict with ADF headers...
#define CFS_DUMMY_NO_ERROR NO_ERROR
#define CFS_DUMMY_NO_DATA NO_DATA
#undef NO_ERROR
#undef NO_DATA
#include <adf/ADF.h>
#include <adfh/ADFH.h>
#undef NO_ERROR
#undef NO_DATA
#define NO_ERROR CFS_DUMMY_NO_ERROR
#define NO_DATA CFS_DUMMY_NO_DATA
#undef CFS_DUMMY_NO_ERROR
#undef CFS_DUMMY_NO_DATA
#endif

#ifdef USE_CCMIO
#include <libccmio/ccmioversion.h>
#endif

#include <boost/version.hpp>

#include <muParserBase.h>

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
extern "C" void ilaver(int*, int*, int*);

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

      ( "numThreads,t", po::value<UInt>()->default_value(1),
        "number of threads used in CFS run. Default 1." )

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

      ( "logConfFile,l", po::value<string>(),
        "name of configuration file for logging streams" )

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
      fs::path simPath = fs::path( varMap_["simName"].as<string>());

      // return only file name without path information
      return simPath.filename().string();
    } else {
      return "";
    }
  }

  fs::path ProgramOptions::GetSimPath() const
   {
     if( varMap_.count( "simName") != 0 ) {

       // get complete path
       //std::cout << "sN='" << varMap_["simName"].as<string>() << "'\n";
       fs::path simPath ( varMap_["simName"].as<string>());
       //std::cout << "abs='" << boost::filesystem::absolute( simPath.parent_path()).string() << "'\n";

       // return path to simulation
       return fs::absolute( simPath.parent_path());
     } else {
       return fs::initial_path().string();
     }
   }

  string ProgramOptions::GetSimPathStr() const
   {
     // get complete path
     fs::path simPath = GetSimPath();

     return simPath.string();
   }

  fs::path ProgramOptions::GetParamFile() const
  {
    if( varMap_.count( "paramFile" ) == 0 ) {
      return GetSimPath() / fs::path(GetSimName()+".xml" );
    } else {
      fs::path paramPath( varMap_["paramFile"].as<string>());
      return fs::complete( paramPath ); //
    }
  }

  string ProgramOptions::GetParamFileStr() const
  {
    fs::path paramPath = GetParamFile();

    return paramPath.string();
  }

  fs::path ProgramOptions::GetLogConfFile() const
  {
    if( varMap_.count("logConfFile") > 0 ) {
      fs::path filePath( varMap_["logConfFile"].as<string>() );
      return fs::system_complete( filePath);
    } else {
      return fs::path();
    }
  }

  string ProgramOptions::GetLogConfFileStr() const
  {
    fs::path filePath = GetLogConfFile();

    return filePath.string();
  }

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
#if defined(WIN32) || defined(__MINGW32__)
      // watch out for leading and closing " in schema string
      boost::trim_if( schema, boost::is_any_of(" \t\"\'") );

      // std::cout << "schema = " << schema << std::endl;
#endif      
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

    return schemaPath.string();
  }

  fs::path ProgramOptions::GetMeshFile() const
  {

    if( varMap_.count( "meshFile") != 0 ) {
      fs::path meshPath( varMap_["meshFile"].as<string>() );
      return fs::system_complete( meshPath );
    } else {
      return fs::path();
    }
  }

  string ProgramOptions::GetMeshFileStr() const
  {
    fs::path meshFile = GetMeshFile();

    return meshFile.string();
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

  bool ProgramOptions::DoDetailedInfo() const
  {
    return varMap_.count("detailed") > 0;
  }


  bool ProgramOptions::IsQuiet() const
  {
    return varMap_.count("quiet") > 0;
  }
  
  UInt ProgramOptions::GetNumThreads() const
  {
#ifdef USE_OPENMP
    if( varMap_.count( "numThreads") != 0 ) {
      return varMap_["numThreads"].as<UInt>();
    }else{
      return 1;
    }
#else
    return 1;
#endif
  }

  void ProgramOptions::PrintHelp( std::ostream& out )
  {
    out << helpMsg_;
  }

  void ProgramOptions::ToInfo(PtrParamNode in) const
  {
    in->Get("problem")->SetValue(GetSimName());
    in->Get("parameterFile")->SetValue(GetParamFileStr());
    in->Get("schemaPath")->SetValue(GetSchemaPathStr());
    in->Get("meshFile")->SetValue(GetMeshFileStr());
    in->Get("logConfFile")->SetValue(GetLogConfFileStr());
    in->Get("detailed")->SetValue(DoDetailedInfo());
    in->Get("MKL_NUM_THREADS")->SetValue(getenv("MKL_NUM_THREADS") != NULL ? getenv("MKL_NUM_THREADS") : "-");

    // cfs information
    in = in->Get("cfs");
    in->Get("version")->SetValue(CFS_VERSION);
    in->Get("name")->SetValue(CFS_NAME);
    in->Get("build")->SetValue(CMAKE_BUILD_TYPE);
    in->Get("svn_revision")->SetValue(CFS_WC_REVISION);
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
          << fg_blue << fn.string() << fg_reset << endl;
    }
    
    out << "CFS_BUILD_HOST:        "
        << fg_blue << CFS_BUILD_HOST << fg_reset << endl

        << "CFS_BUILD_USER:        "
        << fg_blue << CFS_BUILD_USER << fg_reset << endl

        << "CFS_BUILD_DISTRO:      "
        << fg_blue << CFS_BUILD_DISTRO << fg_reset << endl

        << "CFS_WC_REVISION:       "
        << fg_blue << CFS_WC_REVISION << fg_reset << endl

        << "CFS_WC_URL:            "
        << fg_blue << CFS_WC_URL
        << fg_reset << endl;
    
    if( std::string(CFS_WC_TYPE) == "Git" ) 
    {
      out << "CFS_GIT_COMMIT:        "
          << fg_blue << CFS_GIT_COMMIT << fg_reset << endl
        
          << "CFS_GIT_BRANCH:        "
          << fg_blue << CFS_GIT_BRANCH << fg_reset << endl;
    }

    out << endl;

    out << "CFS_CXX_COMPILER_NAME: "
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

#ifndef __MINGW32__
        << "CFS_DISTRO_VER:        "
        << fg_blue << CFS_DISTRO_VER << fg_reset << endl
#endif

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

#ifdef USE_OPENMP
    out << "USE_OPENMP:            "
        << fg_blue  << "YES" << fg_reset << endl;
#else
    out << "USE_OPENMP:            "
        << fg_blue  << "NO" << fg_reset << endl;
#endif

#ifdef DISABLE_OPTIMIZATION
    out << "DISABLE_OPTIMIZATION:  "
        << fg_blue  << "YES" << fg_reset << endl;
#else
    out << "DISABLE_OPTIMIZATION:  "
        << fg_blue  << "NO" << fg_reset << endl;
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

    MKL_Get_Version(reinterpret_cast<MKLVersion*>(&ver));

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

    MKL_Free_Buffers();

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
    ilaver(&major, &minor, &rev);
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
 #else
    out << "USE_ILUPACK:           "
        << fg_blue  << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_SUITESPARSE
    out << endl;
    out << "USE_SUITESPARSE:       "
        << fg_blue << "YES" << fg_reset << endl;
    out << "SUITESPARSE_VERSION:   "
        << fg_blue << SUITESPARSE_MAIN_VERSION << "." << SUITESPARSE_SUB_VERSION << "."
        << SUITESPARSE_SUBSUB_VERSION << " (" << SUITESPARSE_DATE << ") "
        << fg_reset << endl;
    out << "AMD_VERSION:           "
        << fg_blue << AMD_MAIN_VERSION << "." << AMD_SUB_VERSION << "."
        << AMD_SUBSUB_VERSION << " (" << AMD_DATE << ") "
        << fg_reset << endl;
    out << "CHOLMOD_VERSION:       "
        << fg_blue << CHOLMOD_MAIN_VERSION << "." << CHOLMOD_SUB_VERSION << "."
        << CHOLMOD_SUBSUB_VERSION << " (" << CHOLMOD_DATE << ") "
        << fg_reset << endl;
    out << "UMFPACK_VERSION:       "
        << fg_blue << UMFPACK_MAIN_VERSION << "." << UMFPACK_SUB_VERSION << "."
        << UMFPACK_SUBSUB_VERSION << " (" << UMFPACK_DATE << ") "
        << fg_reset << endl;
 #else
    out << "USE_SUITESPARSE:       "
        << fg_blue  << "NO" << fg_reset << endl;
 #endif

    out << endl;
 #ifdef USE_PARDISO
    out << "USE_PARDISO:           "
        << fg_blue << "YES" << fg_reset << endl;
    out << "PARDISO_IMPL:          "
        << fg_blue  << CFS_PARDISO << fg_reset << endl;
 #else
    out << "USE_PARDISO:           "
           << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_LIS
    out << endl;
    out << "USE_LIS:               "
        << fg_blue << "YES" << fg_reset << endl;
    out << "LIS_VERSION:           "
        << fg_blue  << LIS_VERSION << fg_reset << endl;
 #else
    out << "USE_LIS:               "
        << fg_blue  << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_METIS
    std::string metistitle(METISTITLE);
    boost::trim(metistitle);

    out << endl
        << "USE_METIS:             "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_METIS_VERSION:     "
        << fg_blue << metistitle
        << fg_reset << endl;
#else
    out << "USE_METIS:             "
        << fg_blue << "NO" << fg_reset<< endl;
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
    out << "GIDPOST_VERSION:       "
        << fg_blue << GIDPOST_VERSION << fg_reset << endl;
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

#ifdef USE_GMV
    out << "USE_GMV:               "
        << fg_blue << "YES" << fg_reset << endl;
 #else
    out << "USE_GMV:               "
        << fg_blue << "NO" << fg_reset << endl;
 #endif

 #ifdef USE_HDF5
    out << "USE_HDF5:              "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_HDF5_VERSION:      "
        << fg_blue 
        << H5_VERS_MAJOR << "."
        << H5_VERS_MINOR << "."
        << H5_VERS_RELEASE << fg_reset << endl;
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
#ifdef USE_COMSOL
    out << "USE_COMSOL:            "
        << fg_blue << "YES" << fg_reset << endl;
    out << "MINIZIP_VERSION:       "
        << fg_blue << MINIZIP_VERSION 
        << " (for reading zipped .mph file)"
        << fg_reset << endl;
#else
    out << "USE_COMSOL:            "
        << fg_blue << "NO" << fg_reset << endl;
#endif
#ifdef USE_CGNS
    out << "USE_CGNS:              "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CGNS_VERSION:          "
        << fg_blue << (CGNS_VERSION/1000) 
        << "." << ((CGNS_VERSION%1000)/100)
        << ((CGNS_VERSION%100)/10)
        << fg_reset << endl;
#if defined(CGNS_COMPATVERSION)
    out << "CGNS_COMPATVERSION:    "
        << fg_blue << (CGNS_COMPATVERSION/1000) 
        << "." << ((CGNS_COMPATVERSION%1000)/100)
        << ((CGNS_COMPATVERSION%100)/10)
        << fg_reset << endl;
#endif
    char version[1024];
    int error_return = 0;
    ADF_Library_Version(version, &error_return ) ;
    out << "ADF_VERSION:           "
        << fg_blue << version << fg_reset << endl;
    ADFH_Library_Version(version, &error_return ) ;
    out << "ADFH_VERSION:          "
        << fg_blue << version << fg_reset << endl;
#else
    out << "USE_CGNS:              "
        << fg_blue << "NO" << fg_reset << endl;
#endif
#ifdef USE_CCMIO
    out << "USE_CCMIO:             "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CCMIO_VERSION:         "
        << fg_blue << kCCMIOVersionStr << fg_reset << endl;
#else
    out << "USE_CCMIO:             "
        << fg_blue << "NO" << fg_reset << endl;
#endif

    out << endl;

    out << "CFS_XERCES_VERSION:    "
        << fg_blue << XERCES_FULLVERSIONDOT << fg_reset << endl;
    out << "XMLSCHEMA:             ";
    if(progOpts)
      out << fg_blue << progOpts->GetSchemaPath() << fg_reset << endl;
    else
      out << fg_blue << XMLSCHEMA << fg_reset << endl;


    out << endl;

#ifdef USE_CGAL
    out << "USE_CGAL:              "
        << fg_blue << "YES" << fg_reset << endl;
    out << "CFS_CGAL_VERSION:      "
        << fg_blue << QUOTEME(CGAL_VERSION)
        << " (" << CGAL_VERSION_NR
        << ", SVN rev. " << CGAL_SVN_REVISION << ")"
        << fg_reset << endl;
#else
    out << "USE_CGAL:              "
        << fg_blue << "NO" << fg_reset << endl;
#endif
#ifdef USE_LIBFBI
    out << "USE_LIBFBI:            "
        << fg_blue << "YES" << fg_reset << endl;
#else
    out << "USE_LIBFBI:            "
        << fg_blue << "NO" << fg_reset << endl;
#endif
#ifdef USE_FLANN
    out << "USE_FLANN:             "
        << fg_blue << "YES" << fg_reset << endl;
    out << "FLANN_VERSION:         "
        << fg_blue << FLANN_VERSION_
        << fg_reset << endl;
#else
    out << "USE_FLANN:             "
        << fg_blue << "NO" << fg_reset << endl;
#endif
    
    out << endl;
    out << "CFS_BOOST_VERSION:     "
        << fg_blue << 
      (BOOST_VERSION / 100000) << "." <<
      (BOOST_VERSION / 100 % 1000) << "." <<
      (BOOST_VERSION % 100) << fg_reset << endl;
    out << "CFS_ZLIB_VERSION:      "
        << fg_blue << zlibVersion() << fg_reset << endl;
    out << "CFS_BZIP2_VERSION:     "
        << fg_blue << BZ2_bzlibVersion() << fg_reset << endl;
    out << "CFS_MUPARSER_VERSION:  "
        << fg_blue << MUP_VERSION << " " MUP_VERSION_DATE << fg_reset << endl;

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
        << "  X-rated terms. Luckily they don't hear developers conversation with compilers!." << endl
        << endl
        << "11.06, Querulous Quokka" << endl
        << endl
        << "11.12, Scintillant Scapolite" << endl
        << endl
        << "12.02, Branchy Burbank" << endl
        << "  c.f. http://web.mit.edu/invent/iow/burbank.html," << endl
        << "  http://en.wikipedia.org/wiki/Russet_Burbank_potato" << endl
        << endl
        << "13.02, Cocky Coriolis" << endl
        << endl
        << "13.06, Nonmatching Noolbenger" << endl
        << "  CFS++ now has a framework for non-conforming interfaces. Coupling can be done" << endl
        << "  using Mortar or Nitsche formulation." << endl
        << endl
        << "13.10, Rotating Rhino" << endl
        << "  Non-conforming interfaces can now cope with moving grids, especially rotation." << endl
        << "  The animal part of the name shall express that CFS++ is heavier than ever," << endl
        << "  despite lots of code optimization." << endl
        << endl
        << "14.08, Maximale Verwirrung" << endl
        << "  The FE-Space branch is the new trunk and the optimization group starts to add "
        << "  its stuff." << endl
        << endl
        << "15.07, Verlorene Soehne" << endl
        << "  The optimization group is finaly back on the FE-space trunk. Not everyone yet, " << endl
        << "  not all features yet, but the the boys are back in town! :)" << endl
        << endl
        << "15.11, Back To The Future" << endl
        << "  Precompiled CFSDEPS are back and eamc080 is a new mirror server for CFSDEPS." << endl
        << "  Tests now are able to compare info.xml files." << endl
        << endl
        << "16.1, Concurrent Monorail" << endl
        << "  Starting point of making classes thread safe in preparation to parallelize assembly loop." << endl
        << "  Introducing CFSDat program for lightweight, pipeline based data processing." << endl;
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
          << " (rev " << CFS_WC_REVISION << ")" << endl
          << " compiled " << __DATE__
          << " as " << CMAKE_BUILD_TYPE << endl
          << " CFS++ routines use " << NUM_CFS_THREADS << " threads for this run" << endl;
      out << "============================================================"
          << "==========="
          << endl << endl;
    }
  }
}
