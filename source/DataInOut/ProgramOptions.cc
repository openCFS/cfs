// include general defines
#include <def_cfs_stats.hh>

#include <def_cfs_fortran_interface.hh>

#include "ProgramOptions.hh"
#include "ColoredConsole.hh"

#include <sstream>
#include <cstdlib>

#include <boost/bind/bind.hpp>
#include <boost/predef/os.h>
#include <filesystem>
#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "Utils/ToolsFull.hh"

// https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
#if defined(WIN32)
  #include <stdlib.h>
#endif
#if (BOOST_OS_LINUX) // note that this NOT a define (boost sucks :()
  #include <unistd.h>
  #include <limits.h>
#endif
#if (BOOST_OS_MACOS)
  #include <mach-o/dyld.h>
#endif

DEFINE_LOG(progopts, "progOpts")

using std::string;
using std::endl;
using std::cout;
using namespace boost::placeholders;

namespace CoupledField {

  // external definition
  ProgramOptions* progOpts = NULL;

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

  std::filesystem::path ProgramOptions::ObtainCFSRootFromSystem()
  {
    // obtaining the path of the current executable is a non-trivial and non-portable task

    // argv[0] from the main() call has this information for
    // - Windows
    // - Linux and macOs when the executable was called absolutely or relatively.
    //   In the link case only the link name is given.
    // see https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe

    // some systems/shells have the executable in "_"
    // https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
    const char* underline = getenv("_");
    LOG_DBG(progopts) << "OCR: env('_')=" << underline;

    // most Linux systems /proc/self/exe. Plain old Ubunut (14) have this for root only
    // https://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
#if (BOOST_OS_WINDOWS)
    char *exePath;
    if(_get_pgmptr(&exePath) != 0)
      exePath = const_cast<char*>("");
    LOG_DBG(progopts) << "OCR: win get_pgmptr=" << exePath;
#endif
#if (BOOST_OS_LINUX)
	assert(PATH_MAX > 10);
    char exePath[PATH_MAX];
    ssize_t l_len = ::readlink("/proc/self/exe", exePath, sizeof(exePath));
    if(l_len == -1 || l_len == sizeof(exePath))
      l_len = 0;
    exePath[l_len] = '\0';
    LOG_DBG(progopts) << "OCR: /proc/self/exe=" << exePath;
#endif
#if (BOOST_OS_MACOS)
	assert(PATH_MAX > 10);
    char exePath[PATH_MAX];
    // in case of a link _NSGetExecutablePath does not resolve the link
    uint32_t len = sizeof(exePath);
    if (_NSGetExecutablePath(exePath, &len) == 0)
    {
      char *canonicalPath = realpath(exePath, NULL);
      if (canonicalPath == NULL)
        exePath[0] = '\0';
      else
      {
        strncpy(exePath,canonicalPath,len);
        free(canonicalPath);
      }
    }
    else
      exePath[0] = '\0';
    LOG_DBG(progopts) << "OCR: macOS realpath=" << exePath;
#endif

    fs::path root = fs::path(exePath).parent_path().parent_path();
    LOG_DBG(progopts) << "OCR: root cand=" << root.string();
    bool good = fs::exists(root);
    LOG_DBG(progopts) << "OCR: root exits=" << good;
    if(!good)
      throw Exception("guessed exePath " + string(exePath) + " but " + root.string() + " doesn't exist");
    return root;
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
        "information about the openCFS executable" )

      ( "history,H",
        "history of revisions" )

      ( "numThreads,t", po::value<UInt>()->default_value(1),
#ifdef USE_OPENMP
        "set number of threads to use\nalternatively, set env CFS_NUM_THREADS" )
#else
        "not available, not compiled with USE_OPENMP" )
#endif

      ( "meshFile,m", po::value<string>(),
        "name of mesh file for the simulation" )

      ( "paramFile,p", po::value<string>(),
        "name of XML parameter file for the simulation" )

      ( "schemaRoot,s", po::value<string>(),
        "XML schema directory (env CFS_SCHEMA_ROOT)")

      ( "ersatz,x", po::value<string>(),
        "name of ersatz material density file")

      ( "restart,r",
        "read restart file of previous simulation run" )

      ( "detailed,d",
        "detailed output to info.xml")

      ( "printGrid,g",
        "read grid from input and write it to output file" )

      ( "exportGrid,G",
        "export the grid to .info.xml, best with -g")

      ( "equationMap,M",
        "create equation mapping data file (.map)")

      ( "dependencies,c", po::value<string>(),
        "write CMake file with dependencies to given file")

      ( "logConfFile,l", po::value<string>(),
        "logging config file (.xml), only debug build" )

      ( "forceSegFault,f",
        "force a segmentation fault at exceptions")

      ( "quiet,q",
        "more compressed console output (env CFS_QUIET)")

      ( "noColor",
        "turn off colored output")

      ( "id", po::value<string>(),
        "set cfsInfo/header/@id in .info.xml")
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

    // 4) Combine visible and invisible options to commandLine options
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

    // Initialize logging class (read parameters from file if desired)
    // such that we can debug FindSchemaPath()
    std::string confFile = progOpts->GetLogConfFileStr();
    LogConfigurator::ParseLogConfFile(confFile);

    // obtain schema path once and prior --version prints the information
    schemaPath_ = FindSchemaPath();

    // Check if colored output should be switched off
    if( varMap_.count("noColor") != 0) {
      ColoredConsole::suppressed = true;
    }

    // Check for version
    if( varMap_.count("version") != 0  ) 
    {
      PrintHeader(cout);
      PrintVersion( cout, true );
      exit( EXIT_SUCCESS );
    }

    // Check for history
    if(varMap_.count("history") != 0)
    {
      PrintHeader(cout);
      PrintHistory(cout);
      exit( EXIT_SUCCESS );
    }

    // do we want to output dependency as cmake file? Exit afterwards
    if(varMap_.count("dependencies"))
    {
      PrintHeader(cout);
      deps_.WriteCMakeUSE(varMap_["dependencies"].as<string>());
      exit(EXIT_SUCCESS);
    }

    // Check for help
    if( varMap_.count("help") != 0 || varMap_.count("simName") == 0) 
    {
      PrintHeader(cout);
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
       fs::path simPath ( varMap_["simName"].as<string>());

       // return path to simulation
       fs::path parent = simPath.parent_path();
       return parent.empty() ? fs::current_path() : fs::absolute(parent);
    } else {
       return fs::current_path().string();
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
      return fs::absolute( paramPath ); //
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
      return fs::absolute( filePath);
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
    // XMLSCHEMA is set by compile time and will fail on any distributed cfs.

    LOG_DBG(progopts) << "FSP: compiled XMLSCHEMA: " << XMLSCHEMA;

    string schema;

    bool given = false;
    // command line option or CFS_SCHEMA_ROOT environment variable beats CFS_SCHEMA_ROOT
    if(varMap_.count("schemaRoot"))
    {
      schema = varMap_["schemaRoot"].as<string>();
#if defined(WIN32)
      // watch out for leading and closing " in schema string
      boost::trim_if(schema, boost::is_any_of(" \t\"\'"));
#endif      
      LOG_DBG(progopts) << "FSP: arg given: " << schema;
      given = true;
    }
    else
    {
      // from CMakeLists.txt with ${CFS_SOURCE_DIR}/share/xml as default
      // is invalid if the binaries are copied to another location (e.g. binary distribution)
      schema = XMLSCHEMA;
    }

    // make a normalized schemaPath - works also in the non-existing case
    fs::path schemaPath = fs::absolute(schema); // if it did not start with root inserts current working directory, which is clearly nonsense but does not throw an error
    schemaPath = schemaPath.lexically_normal(); // resolves stuff like bla/../fasel. Is depreciated and should work without

    if(fs::exists(schemaPath))
      return schemaPath;

    // now either CFS_SCHEMA_ROOT and command line are not given or wrong and also XMLSCHEMA is invalid
    // generate from cfs root but make a smart message
    if(given)
      std::cout << "Warning: given xml schema path '" + schemaPath.string() + "' is invalid. We construct from executable location\n"; // see .info.xml or --version

    fs::path root = ObtainCFSRootFromSystem();
    root = root.lexically_normal(); // shall be save to remove when the depreciaton hurts
    string exe_schema_root = root.string() + "/share/xml"; // shall work also on Windows

    if(fs::exists(exe_schema_root))
      return exe_schema_root;
    else
      EXCEPTION( "Schema path " << exe_schema_root << " does not exist. Use -s or environment variable CFS_SCHEMA_ROOT." );
  }


  fs::path ProgramOptions::GetMeshFile() const
  {

    if( varMap_.count( "meshFile") != 0 ) {
      fs::path meshPath( varMap_["meshFile"].as<string>() );
      return fs::absolute( meshPath );
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
    if( varMap_.count( "numThreads") != 0 )
      return varMap_["numThreads"].as<UInt>();
    else
      return 1;
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
    in_cfs->Get("compiled")->SetValue(__DATE__);
    in_cfs->Get("exe")->SetValue(exe_);

    // openmp information
    PtrParamNode in_omp = in->Get("openmp");
    #ifdef USE_OPENMP
      in_omp->Get("CFS_NUM_THREADS")->SetValue(CFS_NUM_THREADS);
      if(HasBlasThreadsEnvVariable()) {  // Dependencies.cc for comments
        const char* otherenv = GetBlasThreadsEnvVariable(); // MKL_NUM_THREADS or VECLIB_MAXIMUM_THREADS
        in_omp->Get(otherenv)->SetValue(getenv(otherenv) != NULL ? getenv(otherenv) : "-");
      }
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
  
  void ProgramOptions::PrintVersion(std::ostream & out, bool colorise)
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
      fs::path fn = fs::absolute(progOpts->exe_);
      fn = fn.lexically_normal();
      WriteColoredString(out, trim_size, "executable", fn.string());
      WriteColoredString(out, trim_size, "xml schema", progOpts->GetSchemaPathStr());
      WriteColoredString(out, trim_size, "compiled XMLSCHEMA", XMLSCHEMA);
    }
    else
      WriteColoredString(out, trim_size, "XMLSCHEMA", XMLSCHEMA);
    
    WriteColoredString(out, trim_size, "CFS_BUILD_HOST", CFS_BUILD_HOST);
    WriteColoredString(out, trim_size, "CFS_BUILD_DISTRO", CFS_BUILD_DISTRO);
    WriteColoredString(out, trim_size, "CFS_GIT_BRANCH", CFS_GIT_BRANCH);
    WriteColoredString(out, trim_size, "CFS_GIT_COMMIT", CFS_GIT_COMMIT);
    out << endl;

    WriteColoredString(out, trim_size, "CMAKE_BUILD_TYPE", CMAKE_BUILD_TYPE);
    WriteColoredString(out, trim_size, "C++ compiler", CMAKE_CXX_COMPILER_ID, CMAKE_CXX_COMPILER_VERSION);
    WriteColoredString(out, trim_size, "Fortran compiler", CMAKE_Fortran_COMPILER_ID, CMAKE_Fortran_COMPILER_VERSION);
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
    WriteColoredString(out, trim_size, "sizeof(long)", sizeof(long));
    WriteColoredString(out, trim_size, "sizeof(size_t)", sizeof(size_t));
    WriteColoredString(out, trim_size, "sizeof(void*)", sizeof(void*));
    
    ColoredConsole::colorise = colTmp;
  }

  void ProgramOptions::PrintHistory( std::ostream & out)
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
        << "16.01, Concurrent Monorail" << endl
        << "  Starting point of making classes thread safe in preparation to parallelize assembly loop." << endl
        << "  Introducing CFSDat program for lightweight, pipeline based data processing." << endl
        << endl
        << "18.08, AMGme" << endl
        << "  Contains Klaus' algebraic multigrid." << endl
        << endl
        << "20.09, Frohes Fensterln" << endl
        << "  Brings native Windows support and several features from the fork (thanks to Hermann, Jens and Simon)." << endl
        << "  We are on the track of going open source as openCFS! :)" << endl
        << endl
        << "22.10, Welcome World" << endl
        << "  CFS++ became officially open source and got renamed to openCFS (www.opencfs.org)." << endl
        << "  Feel free to contribute via gitlab.com/openCFS :)" << endl;
  }

  void ProgramOptions::PrintHeader(std::ostream & out)
  {
    if(IsQuiet())
    {
      out << ">> openCFS '" << CFS_VERSION << " " << CFS_NAME << "'"
          << " Compiled: '" << __DATE__ << "'"
          << " Build: '" << CMAKE_BUILD_TYPE << "'"
          << " (" << (deps_.IsDistributable() ? "" : "NOT ") << "distributable)" << endl;
    }
    else
    {
      // CFS_VERSION and CFS_NAME are to be set in source/CMakeLists.txt
      string git_ref(CFS_GIT_COMMIT); // "591da40ea0b41e336e169ee632618fd63313231f (modified)"
      string git_short = git_ref.substr(0,6); // first 6 characters
      if(git_ref.find(' ') != string::npos) // re-add a potential " (modified)"
        git_short += git_ref.substr(git_ref.find(' '));

      out << endl
          << "=======================================================================" << endl;
      out << " openCFS - Coupled Field Simulation" << endl << endl
          << " v. " << CFS_VERSION << " - '" << CFS_NAME << "'" << " " << git_short << endl
          << " compiled " << __DATE__  << " as " << CMAKE_BUILD_TYPE
          << " (" << (deps_.IsDistributable() ? "" : "NOT ") << "distributable)" << endl << endl;

    }
  }

  void ProgramOptions::PrintNumThreads(std::ostream& out, bool quiet)
  {
    string omp = getenv("OMP_NUM_THREADS") != NULL ? string(getenv("OMP_NUM_THREADS")) : "-";
    const char* otherenv = GetBlasThreadsEnvVariable(true); // MKL_NUM_THREADS or VECLIB_MAXIMUM_THREADS
    string other = getenv(otherenv) != NULL ? string(getenv(otherenv)) : "-";
#ifdef USE_OPENMP
    if(quiet)
    {
      out << ">> OpenMP *_NUM_THREADS: CFS=" << CFS_NUM_THREADS  << ", OMP=" << omp;
      if(HasBlasThreadsEnvVariable())
        out << ", " << GetBlasThreadsEnvVariable(false) << "=" << other; // MKL or VECLIB
      out << endl;
    }
    else
    {
      out << " openCFS is using " << CFS_NUM_THREADS << " OpenMP threads";
      string cfs = std::to_string(CFS_NUM_THREADS);
      if((cfs != omp) || (cfs != other))
        out << " and OMP_NUM_THREADS=" << omp << " and " << otherenv << "=" << other;
      out << endl;
      out << "=======================================================================" << endl << endl;
    }
#endif
  }
}
