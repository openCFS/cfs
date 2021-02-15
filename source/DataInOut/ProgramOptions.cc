// include general defines
#include <def_cfs_stats.hh>

#include <def_cfs_fortran_interface.hh>

#include "ProgramOptions.hh"
#include "ColoredConsole.hh"

#include <sstream>
#include <cstdlib>

#include <boost/bind/bind.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Utils/tools.hh"

using std::string;
using std::endl;
using std::cout;
using namespace boost::placeholders;

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
        "number of threads used in CFS run, default = 1 (env CFS_NUM_THREADS)" )

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

      ( "detailed,d",
        "detailed output to info.xml")

      ( "printGrid,g",
        "read grid from input and write it to output file" )

      ( "exportGrid,G",
        "export the grid to the info.xml file, best with -g")

      ( "equationMap,M",
        "create a .map file with details about the equation mapping")

      ( "dependencies,c", po::value<string>(),
        "write CMake file with dependencies to given file")

      ( "logConfFile,l", po::value<string>(),
        "name of configuration file for logging streams, only works for DEBUG builds" )

      ( "forceSegFault,f",
        "force a segmentation fault at exceptions")

      ( "quiet,q",
        "more compressed console output (env CFS_QUIET)")

      ( "id", po::value<string>(),
        "set the provided value in info.xml as cfsInfo/header/@id")

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

    // obtain schema path once and prior --version prints the information
    schemaPath_ = FindSchemaPath();

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
    if(varMap_.count("history") != 0)
    {
      GetHeaderString(cout);      
      GetHistoryString(cout);
      exit( EXIT_SUCCESS );
    }

    // do we want to output dependency as cmake file? Exit afterwards
    if(varMap_.count("dependencies"))
    {
      GetHeaderString(cout);
      deps_.WriteCMakeUSE(varMap_["dependencies"].as<string>());
      exit(EXIT_SUCCESS);
    }

    // Check for help
    if( varMap_.count("help") != 0 || varMap_.count("simName") == 0) 
    {
      GetHeaderString(cout);      
      cout << helpMsg_;
      exit(EXIT_SUCCESS);
    }
  }

  string ProgramOptions::EnvironmentNameMapper(const string& var )
  {
    string ret;

    if(var == "CFS_SCHEMA_ROOT") ret = "schemaRoot";
    if(var == "CFS_NO_COLOR")    ret = "noColor";
    if(var == "CFS_QUIET")       ret = "quiet";
    if(var == "CFS_NUM_THREADS") ret = "numThreads";
    
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
#ifdef NDEBUG
      WARN("logging only works for DEBUG builds");
#endif
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

  string ProgramOptions::GetId() const
  {
    return varMap_.count("id") != 0 ? varMap_["id"].as<string>() : "";
  }

  string ProgramOptions::GetErsatzMaterialStr() const
  {
    if(varMap_.count("ersatz") != 0)
      return varMap_["ersatz"].as<string>();
    else
      return "";
  }

  fs::path ProgramOptions::FindSchemaPath() const
  {
    string schema;

    // If the user specified a path on the command line or the CFS_SCHEMA_ROOT variable, use it instead.
    if(varMap_.count("schemaRoot"))
    {
      schema = varMap_["schemaRoot"].as<string>();
#if defined(WIN32)
      // watch out for leading and closing " in schema string
      boost::trim_if(schema, boost::is_any_of(" \t\"\'"));
#endif      
    }
    else
    {
      // from CMakeLists.txt with ${CFS_SOURCE_DIR}/share/xml as default
      // is invalid if the binaries are copied to another location (e.g. binary distribution)
      schema = XMLSCHEMA;
    }

    // make a normalized schemaPath - works also in the non-existing case
    fs::path schemaPath = fs::system_complete(schema); // if it did not start with root inserts current working directory, which is clearly nonsense but does not throw an error
    schemaPath.normalize(); // resolves stuff like bla/../fasel. Is depreciated and should work without

    // try to obtain the schema root from argv[0]. Use it only if schema is invalid
    // This is the full path for Windows but on Linux and macOS this is either the full path or the pure link name (e.g. "cfs_rel")
    // For the full path we are in .../<cfs_root>/bin/cfsbin or .../<cfs_root>/bin/MACOSX_11.1_X86_64/cfsbin
    // for the link case, parent_path() is ".." and parent_path().parent_path() is "", all expection save
    fs::path exe_schema_root;
    if(fs::path(exe_).parent_path().filename() == "bin")
      exe_schema_root = fs::path(exe_).parent_path().parent_path();
    else if(fs::path(exe_).parent_path().parent_path().filename() == "bin") // remove when getting rid of CFS_ARCH_STR folders
      exe_schema_root = fs::path(exe_).parent_path().parent_path().parent_path();
    exe_schema_root.normalize(); // shall be save to remove when the depreciaton hurts
    exe_schema_root = exe_schema_root.string() + "/share/xml"; // shall work also on Windows

    if(fs::exists(schemaPath))
      return schemaPath;
    if(fs::exists(exe_schema_root)) {
      std::cout << "Warning: given xml schema path '" + schemaPath.string() + "' is invalid, use extraction from executable.\n"; // given in .info.xml
      return exe_schema_root;
    }

    EXCEPTION( "Schema path " << schemaPath << " does not exist. Use -s or environment variable CFS_SCHEMA_ROOT." );
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

  bool ProgramOptions::DoEquationMapping() const
  {
    return varMap_.count("equationMap") > 0;
  }

  bool ProgramOptions::GetRestart() const
  {
    return (varMap_.count( "restart") > 0);
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

    // cfs information
    PtrParamNode in_cfs = in->Get("cfs");
    in_cfs->Get("version")->SetValue(CFS_VERSION);
    in_cfs->Get("name")->SetValue(CFS_NAME);
    in_cfs->Get("build")->SetValue(CMAKE_BUILD_TYPE);
    in_cfs->Get("git_commit")->SetValue(CFS_GIT_COMMIT);
    in_cfs->Get("git_branch")->SetValue(CFS_GIT_BRANCH);
    in_cfs->Get("exe")->SetValue(exe_);

    // openmp information
    PtrParamNode in_omp = in->Get("openmp");
    #ifdef USE_OPENMP
      in_omp->Get("CFS_NUM_THREADS")->SetValue(CFS_NUM_THREADS);
      in_omp->Get("MKL_NUM_THREADS")->SetValue(getenv("MKL_NUM_THREADS") != NULL ? getenv("MKL_NUM_THREADS") : "-");
      in_omp->Get("OMP_NUM_THREADS")->SetValue(getenv("OMP_NUM_THREADS") != NULL ? getenv("OMP_NUM_THREADS") : "-");
    #endif
    #ifndef USE_OPENMP
      in_omp->SetValue("compiled without");
    #endif
  }

  void ProgramOptions::WriteColoredString(std::ostream& out, int trim_size, const string& head, const string& data, const string& opt1, const string& opt2, const string& opt3)
  {
    out << std::left << std::setw(trim_size) << (head + ": ");
    out << fg_blue << data;
    if(opt1 != "")
      out << ", " << opt1;
    if(opt2 != "")
      out << ", " << opt2;
    if(opt3 != "")
      out << ", " << opt3;
    out << fg_reset << endl;
  }

  void ProgramOptions::WriteColoredString(std::ostream& out, int trim_size, const string& head, int data)
  {
    WriteColoredString(out, trim_size, head, std::to_string(data));
  }
  
  void ProgramOptions::GetVersionString(std::ostream & out, bool colorise)
  {
    // local instance as GetVersionString() is static.
    // If you want to know deps.IsDistributable() outside, simply create your own object
    Dependencies deps;

    bool colTmp = ColoredConsole::colorise;
    ColoredConsole::colorise = colorise;
    
    const int trim_size = 26;
    WriteColoredString(out, trim_size, "CFS_VERSION", CFS_VERSION);
    WriteColoredString(out, trim_size, "CFS_NAME", CFS_NAME);

    if(progOpts) 
    {      
      fs::path fn = fs::system_complete(progOpts->exe_);
      fn.normalize();
      WriteColoredString(out, trim_size, "CFS_EXECUTABLE", fn.string());
      WriteColoredString(out, trim_size, "XMLSCHEMA", progOpts->GetSchemaPathStr());
    }
    else
      WriteColoredString(out, trim_size, "XMLSCHEMA", XMLSCHEMA);
    
    WriteColoredString(out, trim_size, "CFS_BUILD_HOST", CFS_BUILD_HOST);
    WriteColoredString(out, trim_size, "CFS_BUILD_USER", CFS_BUILD_USER);
    WriteColoredString(out, trim_size, "CFS_BUILD_DISTRO", CFS_BUILD_DISTRO);
    WriteColoredString(out, trim_size, "CFS_GIT_COMMIT", CFS_GIT_COMMIT);
    WriteColoredString(out, trim_size, "CFS_GIT_BRANCH", CFS_GIT_BRANCH);
    WriteColoredString(out, trim_size, "CFS_WC_REVISION", CFS_WC_REVISION);
    out << endl;

    WriteColoredString(out, trim_size, "CMAKE_BUILD_TYPE", CMAKE_BUILD_TYPE);
    WriteColoredString(out, trim_size, "C++ compiler", CFS_CXX_COMPILER_NAME, CFS_CXX_COMPILER_VER);
    WriteColoredString(out, trim_size, "Fortran compiler", CFS_FORTRAN_COMPILER_NAME, CFS_FORTRAN_COMPILER_VER);
    WriteColoredString(out, trim_size, "CFS_DISTRO", CFS_DISTRO, CFS_DISTRO_VER);
    WriteColoredString(out, trim_size, "CFS_ARCH", CFS_ARCH);

    string cxxFlags = CMAKE_CXX_FLAGS;
    string ldFlags  = CFS_LINKER_FLAGS;
#ifndef NDEBUG
    cxxFlags += " " + string(CMAKE_CXX_FLAGS_DEBUG);
    ldFlags  += " " + string(CMAKE_EXE_LINKER_FLAGS_DEBUG);
#else
    cxxFlags += " " + string(CMAKE_CXX_FLAGS_RELEASE);
    ldFlags  += " " + string(CMAKE_EXE_LINKER_FLAGS_RELEASE);
#endif
    WriteColoredString(out, trim_size, "COMPILE_FLAGS", cxxFlags);
    WriteColoredString(out, trim_size, "LINK_FLAGS", ldFlags);
    out << endl;

    // deps.Dump();

    WriteColoredString(out, trim_size, "external licenses",
        string("distributable=") + (deps.IsDistributable() ? "YES": "NO"),
        string("commercial=") + (deps.HasLicense(Dependencies::COMMERCIAL) ? "YES": "NO"),
        string("closed=") + (deps.HasLicense(Dependencies::CLOSED) ? "YES": "NO"),
        string("gpl=") + (deps.HasLicense(Dependencies::GPL) ? "YES": "NO"));
    out << endl;

    // the actual licenses are not displayed as one needs to be accurate then and
    // actually almost any lib has a own licencse. For the standard licenses one
    // would have to give the flavours (number, new-bsd, ...)
    for(auto& dep: deps.data)
      if(dep.IsSwitchable())
      {
        // build artefacts are not relevant for the cfs binary itself for user output
        if(dep.IsBuildOption())
          continue;
        if(dep.active)
          WriteColoredString(out, trim_size, dep.cmake, dep.version != "" ? dep.version : "unknown version", dep.date, dep.comment);
        else
          WriteColoredString(out, trim_size, dep.cmake, "NO");
      }
    out << endl;

    for(auto& dep: deps.data)
      if(!dep.IsSwitchable())
      {
        assert(dep.active);
        WriteColoredString(out, trim_size, dep.name, dep.version != "" ? dep.version : "YES", dep.date, dep.comment, dep.cmake);
      }
    out << endl;

    WriteColoredString(out, trim_size, "sizeof(int)", sizeof(int));
    WriteColoredString(out, trim_size, "sizeof(Integer)", sizeof(Integer));
    WriteColoredString(out, trim_size, "sizeof(size_t)", sizeof(size_t));
    WriteColoredString(out, trim_size, "sizeof(void*)", sizeof(void*));
    
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
        << "  The optimization group is finally back on the FE-space trunk. Not everyone yet, " << endl
        << "  not all features yet, but the the boys are back in town! :)" << endl
        << endl
        << "15.11, Back To The Future" << endl
        << "  Precompiled CFSDEPS are back and eamc080 is a new mirror server for CFSDEPS." << endl
        << "  Tests now are able to compare info.xml files." << endl
        << endl
        << "16.1, Concurrent Monorail" << endl
        << "  Starting point of making classes thread safe in preparation to parallelize assembly loop." << endl
        << "  Introducing CFSDat program for lightweight, pipeline based data processing." << endl
        << endl
        << "18.8, AMGme" << endl
        << "  Contains Klaus' algebraic multigrid." << endl
        << endl
        << "20.9, Frohes Fensterln" << endl
        << "  Brings native Windows support and several features from the fork (thanks to Hermann, Jens and Simon)." << endl
        << "  We are on the track of going open source as openCFS! :)" << endl;
  }

  void ProgramOptions::GetHeaderString(std::ostream & out)
  {
    string omp = getenv("OMP_NUM_THREADS") != NULL ? string(getenv("OMP_NUM_THREADS")) : "-";
    string mkl = getenv("MKL_NUM_THREADS") != NULL ? string(getenv("MKL_NUM_THREADS")) : "-";
    if(IsQuiet())
    {
      out << ">> CFS++ '" << CFS_VERSION << " " << CFS_NAME << "'"
          << " Compiled: '" << __DATE__ << "'"
          << " Build: '" << CMAKE_BUILD_TYPE << "'"
          << " (" << (deps_.IsDistributable() ? "" : "NOT ") << "distributable)" << endl;
      #ifdef USE_OPENMP
        out << ">> OpenMP *_NUM_THREADS: CFS=" << CFS_NUM_THREADS  << ", OMP=" << omp << ", MKL=" << mkl << endl;
      #endif
    }
    else
    {
      // CFS_VERSION and CFS_NAME are to be set in source/CMakeLists.txt
      out << endl
          << "=======================================================================" << endl;
      out << " CFS++ - Coupled Field Simulation" << endl << endl
          << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'" << " " << CFS_WC_REVISION << endl
          << " compiled " << __DATE__  << " as " << CMAKE_BUILD_TYPE
          << " (" << (deps_.IsDistributable() ? "" : "NOT ") << "distributable)" << endl
          << " CFS++ routines use " << CFS_NUM_THREADS << " threads for this run."
          << " OMP/ MKL threads: " << omp << "/ " << mkl << endl;
      out << "=======================================================================" << endl << endl;
    }
  }
}
