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


#ifdef USE_MKL_LIB
#include <mkl_service.h>
#endif

#include "General/environment.hh"
#include "General/exception.hh"

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
       << "./nacs [param] <name of simulation run>\n\n"
       << " where [param] is one or more of the following parameters";
    
    
    // 1) Define "visible" commandLine options
    // -----------------------------------------
    po::options_description cmdVisible( os.str() );
    
    cmdVisible.add_options()
      ( "help,h", 
        "print this usage information" )
      
      ( "printGrid,g", 
        "only read grid from mesh-file and write it to output file" )

      ( "paramFile,p", po::value<std::string>(), 
        "name of XML parameter file for the simulation" )
      
      ( "meshFile,m", po::value<std::string>(),
        "name of mesh file for the simulation" )
      
#ifdef USE_SCRIPTING
      ( "scriptFile,e", po::value<std::string>(), 
        "name of script file to be evaluated" )
#endif

#ifdef DEBUG
      ( "traceDepth,t", po::value<UInt>(),
        "depth of function tracing" )
#endif
      
      ( "doProfile,d", 
        "turns on generation of profiling information" )
      
      ( "restart,r", 
        "read restart file of previous simulation run" )
      
      ( "dumpStats,u",
        "dump information about the CFS++ executable and license" )
      ;

    
    // Store help description in variable
    std::stringstream helpStream;
    helpStream << cmdVisible << std::endl;
    helpMsg_ = helpStream.str();


    // 2) Define "invisible" commandLine options
    // -----------------------------------------
    po::options_description cmdInvisible( os.str() );
    
    cmdInvisible.add_options()
      ("forceSegFault","force a segmentation fault instead of throwing an exception")

      ("simName", po::value<std::string>(), 
       "simulation name")

      ( "writeSkeleton,w",
        "write skeleton of XML file for subsequent simulation" )
      ;

    // make the input file the only allowed positional parameter
    po::positional_options_description p;
    p.add( "simName", 1 );


    // 3) Define environment options
    // -----------------------------------------
    po::options_description envOptions("" );
    envOptions.add_options()
      ("root",po::value<std::string>(), 
       "Path to base directory of CFS++ installation")
      
      ("schemaRoot", po::value<std::string>(),
       "Path to schema definitions of CFS++ installation")
       
      ("ompNumThreads", po::value<int>(),
       "Number of processors/threads for OpenMP")
      ;

    // 4) Combine visible and invisble options to commandLine options
    // -----------------------------------------
    po::options_description cmdLineOptions;
    cmdLineOptions.add(cmdVisible).add(cmdInvisible);

    
    // 5) Parse command line
    // -----------------------------------------
    po::store( po::command_line_parser( args_ ).
               options( cmdLineOptions ).positional( p ).run(), 
               varMap_ );
    
    // 6) Parse environment
    // -----------------------------------------

    // first, define function pointer, which maps environment variables
    // (e.g. CFS_ROOT) to local string paramter representation (e.g. root)
    boost::function1<std::string, std::string> name_mapper;
    name_mapper = boost::bind(&ProgramOptions::EnvironmentNameMapper, *this, _1);
    po::store( po::parse_environment( envOptions, name_mapper ),
               varMap_ );

    po::notify( varMap_ );

    // Check for dumpStats
    if( varMap_.count("dumpStats") != 0  ) {
      DumpStats( std::cout, true );
      exit( EXIT_SUCCESS );
    }

    // Check for help
    if( varMap_.count("help") != 0) {
      std::cout << helpMsg_;
      exit(EXIT_SUCCESS);
    }
    
    // If no argument was given, print additional information
    if( varMap_.count("simName") == 0 ) {
      std::cout << "\nPlease specify a simulation file (without extension)\n\n\n";
      exit(EXIT_SUCCESS);
    }

    // Print all found information
//     po::variables_map::iterator it = varMap_.begin();
//     for( ; it != varMap_.end(); it++ ) {
//       std::cerr << "key: " << it->first ;
//       try {
//         std::cerr << ", value: " << any_cast<std::string>((it->second).value()) << std::endl;
//       } catch( std::exception& ex ) {}
      
//     }
  }
    
  std::string ProgramOptions::EnvironmentNameMapper( std::string envVarName )  {
    std::string ret;
    
    if( envVarName == "CFS_ROOT" )  {
      ret = "root";
    }
    if( envVarName == "CFS_SCHEMA_ROOT" ) {
      ret = "schemaRoot";
    }
    
    if( envVarName == "OMP_NUM_THREADS" ) {
      ret = "ompNumThreads";
    }

    return ret;
  }
  
  
  bool ProgramOptions::GetPrintHelp() 
  { 
    return (varMap_.count("help") > 0);
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
    

    if( varMap_.count( "schemaRoot" ) ) {
      schema = varMap_[ "schemaRoot" ].as<std::string>();
    }
    else if( varMap_.count( "root" ) ) {
      std::string pathsep;
      pathsep = fs::path("/").native_directory_string();

      schema.append(varMap_[ "root" ].as<std::string>());
      schema.append(pathsep);
      schema.append("xml");
    }
    else
      schema = XMLSCHEMA;
      
    try 
    {
      fs::path schemaPath = fs::system_complete(schema);
      schemaPath.normalize();
      return schemaPath;
    } catch (fs::filesystem_error& ex)
    {
      EXCEPTION(ex.what());
    }    
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
  
   
  UInt ProgramOptions::GetTraceDepth() const 
  {

    if( varMap_.count( "traceDepth" ) > 0 ) {
      return varMap_["traceDepth"].as<UInt>();
    } else {
      return 0;
    }
  }

  
  bool ProgramOptions::GetDumpStats() const 
  {

    return (varMap_.count( "dumpStats") > 0);
  }
  

#ifdef DEBUG
  
  bool ProgramOptions::GetForceSegFault() const 
  {

    return (varMap_.count( "forceSegFault") > 0);
  }
#endif
   
  void ProgramOptions::PrintHelp( std::ostream& out ) 
  {

    out << helpMsg_;
  }
  
  void ProgramOptions::PrintParams( std::ostream &out, bool colorise) {
    bool colTmp = ColoredConsole::colorise;
    ColoredConsole::colorise = colorise;

    // Print knwon parameter to stream
    out << fg_blue
        << "\n\nValues of command line / environment "
        << "parameters (including defaults):\n\n"
        << fg_reset

        << " name of simulation run = "
        << fg_blue
        << GetSimName()
        << fg_reset << '\n'

        << " parameter file = "
        << fg_blue
        << GetParamFile().string()
        << fg_reset << '\n'

        << " path to XSD schema definition = "
        << fg_blue
        << GetSchemaPath().string()
        << fg_reset << '\n'
      
        << " mesh file = "
        << fg_blue
        << GetMeshFile().string()
        << fg_reset << '\n'

#ifdef USE_SCRIPTING
        << " scriptFileName = "
        << fg_blue
        << GetScriptFile().string()
        << fg_reset << '\n'
#endif

        << " function trace depth = "
        << fg_blue
        << GetTraceDepth()
        << fg_reset << '\n'

        << " print grid only = "
        << fg_blue
        << std::boolalpha
        << ( GetPrintGrid() == true )
        << fg_reset << '\n'

        << " perform profiling = "
        << fg_blue
        << std::boolalpha
        << ( GetDoProfile() == true )
        << fg_reset << '\n'

        << " perform restart = "
        << fg_blue
        << std::boolalpha
        << ( GetRestart() == true )
        << fg_reset << '\n'

     

#ifdef DEBUG
        << " force creation of segfault = "
        << fg_blue
        << std::boolalpha
        << ( GetForceSegFault() == true )
        << fg_reset << "\n"
#endif

        << " dump xml parameters = "
        << fg_blue
        << std::boolalpha
        << ( GetDumpStats() == true )
        << fg_reset << "\n\n";

    ColoredConsole::colorise = colTmp;
  }

  void ProgramOptions::DumpStats( std::ostream & outstr, bool colorise) {
    GetDumpString(outstr, colorise);
  }

  void ProgramOptions::GetDumpString( std::ostream & outstr,
                                      bool colorise)
  {
    bool colTmp = ColoredConsole::colorise;
    ColoredConsole::colorise = colorise;

    std::string build_type = CMAKE_BUILD_TYPE;

    outstr << "CFS_VERSION:           " 
           << fg_blue << CFS_VERSION << fg_reset << std::endl

           << "CFS_BUILD_HOST:        " 
           << fg_blue << CFS_BUILD_HOST << fg_reset << std::endl

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
    
    if(build_type == "DEBUG") 
    {
      outstr << "COMPILE_FLAGS:       " 
             << fg_blue  << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_DEBUG << fg_reset << std::endl;
      
      outstr << "LINK_FLAGS:          " 
             << fg_blue  << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_DEBUG
             << fg_reset
             << std::endl << std::endl;
    }
    else
    {
      outstr << "COMPILE_FLAGS:       " 
             << fg_blue  << CMAKE_CXX_FLAGS
             << " " << CMAKE_CXX_FLAGS_RELEASE << fg_reset << std::endl;

      outstr << "LINK_FLAGS:          " 
             << fg_blue << CMAKE_EXE_LINKER_FLAGS
             << " " << CMAKE_EXE_LINKER_FLAGS_RELEASE
             << fg_reset << std::endl << std::endl;
    }

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
#else
    outstr << "USE_BLAS:              "
           << fg_blue  << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_LAPACK
    outstr << "USE_LAPACK:            "
           << fg_blue << "YES" << fg_reset << std::endl;
#else
    outstr << "USE_LAPACK:            "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

 #ifdef USE_ILUPACK
    outstr << "USE_ILUPACK:           "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_ILUPACK:           "
           << fg_blue  << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_PARDISO
    outstr << "USE_PARDISO:           "
           << fg_blue << "YES" << fg_reset << std::endl;
 #else
    outstr << "USE_PARDISO:           "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_MKL_LIB
    outstr << std::endl
           << "USE_MKL_LIB:           "
           << fg_blue << "YES" << fg_reset << std::endl;

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
 #else
    outstr << std::endl
           << "USE_MKL_LIB:           "
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
    outstr << "CFS_TCL_VERSION:      " 
           << fg_blue << CFS_TCL_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_SCRIPTING_TCL:     "
           << fg_blue << "NO" << fg_reset << std::endl;
 #endif

 #ifdef USE_SCRIPTING_PYTHON
    outstr << "USE_SCRIPTING_PYTHON:  "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_PYTHON_VERSION:   " 
           << fg_blue << CFS_PYTHON_VERSION << fg_reset << std::endl;
 #else
    outstr << "USE_SCRIPTING_PYTHON:  "
           << fg_blue << "NO" << fg_reset << std::endl;
#endif

    outstr << std::endl;
    
 #ifdef USE_GIDPOST
    outstr << "USE_GIDPOST:           "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_GIDPOST_VERSION:  " 
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
    outstr << "CFS_HDF5_VERSION:     " 
           << fg_blue << CFS_HDF5_VERSION << fg_reset << std::endl;
    outstr << "CFS_ZLIB_VERSION:     " 
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
    outstr << "CFS_BOOST_VERSION:    " 
           << fg_blue << CFS_BOOST_VERSION << fg_reset << std::endl;
    outstr << std::endl;
    outstr << "CFS_METIS_VERSION:    " 
           << fg_blue << CFS_METIS_VERSION << fg_reset << std::endl;
    outstr << std::endl;

    
#ifdef USE_XERCES
    outstr << "USE_XERCES:            "
           << fg_blue << "YES" << fg_reset << std::endl;
    outstr << "CFS_XERCES_VERSION:   " 
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
  

}
