// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_PROGRAM_OPTIONS_HH
#define FILE_CFS_PROGRAM_OPTIONS_HH

#include <boost/program_options.hpp>
#include <filesystem>

// Required for the CFS own data types
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Dependencies.hh"

namespace po = boost::program_options;
namespace fs = std::filesystem;

// Helper function to safely call fs::absolute() with empty path handling
// std::filesystem::absolute() throws on empty paths, unlike boost::filesystem
inline fs::path safe_absolute(const fs::path& p) {
  return p.empty() ? fs::current_path() : fs::absolute(p);
}

namespace CoupledField
{
  class ProgramOptions;

  //! Define global instance of this class
  extern ProgramOptions* progOpts;

  class ProgramOptions
  {

  public:

    //! Standard constructor
    ProgramOptions( Integer argc, const char **argv );

    ~ProgramOptions();

    /** Gather information from command line and environment
     * Logging is initialized based on -l and schemaPath is obtained */
    void ParseData();

    /** static general purpose helper function which tries to obtain the root of the
     * cfs installation.
     * @return absolute file path to root of bin/share/lib or exception */
    static std::filesystem::path ObtainCFSRootFromSystem();

    // =======================================================================
    // QUERY METHODS FOR PARAMETERS
    // =======================================================================

    //@{
    //! \name Query methods for command line parameters
    //! For each command line parameter there exists an associated query
    //! method. If no command line value was specified than the query method
    //! will return the corresponding default, if it exists.

    //! Return base name of simulation run (without path information)

    //! This method can be used to query the name of the current simulation
    //! run, which may be something like e.g. plate3D. This name is used by
    //! CFS internally as basename for the generation of several (default)
    //! names for input and output files.
    std::string GetSimName() const;

    //! Return path to simulation files

    //! This method can be used to query the directory of the current,
    //! i.e. where the output files will be created.
    fs::path GetSimPath() const;
    std::string GetSimPathStr() const;

    //! Return name of XML parameter file (including path)

    //! This method can be used to query the name of the parameter file in
    //! XML format that contains the steering parameters for the simulation.
    fs::path GetParamFile() const;
    std::string GetParamFileStr() const;

    //! Return name of log configuration file

    //! This method returns the full name of an (optional) log configuration 
    //! file, i.e. a xml file, which contains the module names, log levels
    //! and destination, where the log stream gets logged to. 
    //!  If no filename was provided an empty string is returned.
    fs::path GetLogConfFile() const;
    std::string GetLogConfFileStr() const;

    /** Return the optional ersatz  material density file
     * @return "" if nothing given. */
    std::string GetErsatzMaterialStr() const;

    //! Return path to XML schema file

    //! This method can be used to query the path to the XML schema file
    //! used by validating XML parsers to verify the formal correctness
    //! of the XML parameter file.
    //! \note
    //! - There is currently no way to specify the name of the schema
    //!   file itself. This must be called CFS.xsd!
    //! - This path is also used to locate the default XML-file that is
    //!   currently still needed by the XMLParamHandler.
    fs::path GetSchemaPath() const { return schemaPath_; }
    std::string GetSchemaPathStr() const { return schemaPath_.string(); }

    //! Return name of mesh file (including path)

    //! This method can be used to query the name of the mesh file containing
    //! the description of the FEM mesh for the simulation.
    fs::path GetMeshFile() const;
    std::string GetMeshFileStr() const;

    /** the value provided by --id to be written to the cfsInfo/header/@id in the info.xml.
     * This way we can easily postproc a parameter study where only the mesh or material file content is changed */
    std::string GetId() const;


    //! Return printGrid flag

    //! This method can be used to query the status of the printGrid flag.
    //! By specifying this flag one instructs the executable to do not
    //! perform an actual simulation, but to only import the grid and
    //! re-export it to an output file in the format specified in the XML
    //! parameter file.
    bool GetPrintGrid() const;

    /** exports the grid to the info.xml file.
     * Might get really big!! */
    bool DoExportGrid() const;

    /** shall a .map be created using StdPDE::CreateEquationMapFile() */
    bool DoEquationMapping() const;

    //! Return Restart flag

    //! This method can be used to query the status of the restart flag.
    //! If this flag is true the simulation restarts from an previous state.
    bool GetRestart() const;

    //! Return forceSegfault flag

    //! This method can be used to query the status of the forceSegFault flag.
    //! If this flag is set now exception will be thrown, but a segmentation
    //! fault will be forced instead, which enables one to use a debugger
    //! to get a stack trace.
    bool GetForceSegFault() const;
    //@}

    /** Also more detailed info.xml output as with DoListMapping */
    bool DoDetailedInfo() const;

    /** Is cfs invoked with the quite flag to compress console output to just a minimu. */
    bool IsQuiet() const;

    /** Get number of threads for CFS supplied on the command line.
     * When not USE_OPENMP this is always 1 */
    UInt GetNumThreads() const;

    // =======================================================================
    // AUXILLIARY METHODS FOR OUTPUTTING INFORMATION
    // =======================================================================

    //@{ \name Auxilliary methods for outputting information

    //! Print help information to command line
    void PrintHelp( std::ostream &out );

    /** Write the command line options to the info.xml file */
    void ToInfo(PtrParamNode in) const;

    /** collects all available data to the string 
     *  It containts valuable information about the executable like the 
     *  distro on which it was built on, the compiler it was built with and so on. */
    static void PrintVersion(std::ostream & outstr, bool colorise);
    
    /** Major release history notes. From Dec. 08 */
    static void PrintHistory(std::ostream& out);
    
    /** This gives the head line of openCFS printed to the given string (usually cout) */
    void PrintHeader(std::ostream& out);

    /** prints a line with comprehensive CFS_/OMP_[/MKL_NUM_THREADS] or /VECLIB_... info.
     * In not USE_OPENMP case nothing is written.
     * @param quiet the CFS_QUIET (ProgramOptions::IsQuiet()) value */
    static void PrintNumThreads(std::ostream& out, bool quiet);
    // @}

  protected:


    /** obtain root path to the schema directory.
     * This is "<base>/share/xml" such that stuff like "/CFS-Simulation/CFS.xsd" can be added.
     *
     * This are the rules to obtain the path:
     * 1) use program option -s or --schemaRoot
     * 2) use the environment variable CFS_SCHEMA_ROOT
     * 3) use the compile time information XMLSCHEMA (advanced cmake option)
     * 4) if this all fails, try to obtain from runtime location (Linux, Windows, does not work on macOS)  */
    fs::path FindSchemaPath() const;

    /** schema root path obtained by FindSchemaPath */
    fs::path schemaPath_;

    //! Command line arguments as vector
    std::vector<std::string> args_;

    //! Path to executable. This might include the path (Linux, Windows) but does not for macOS
    std::string exe_;

    //! Help message as string
    std::string helpMsg_;

    //! Boost's argument map
    po::variables_map varMap_;

  private:
    //! Helper function for mapping environment variables to internal names
    std::string EnvironmentNameMapper( const std::string& envVarName );

    /** Helper for colored output */
    static void WriteColoredString(std::ostream& out, int trim_size, const string& head, const string& data, const string& opt1 = "", const string& opt2 = "", const string& opt3 = "");

    static void WriteColoredString(std::ostream& out, int trim_size, const string& head, int data);

    // static PrintVersion() has own variant
    Dependencies deps_;
  };
}
#endif
