// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_PROGRAM_OPTIONS_HH
#define FILE_CFS_PROGRAM_OPTIONS_HH

#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>


// Include defs
#include <def_build_type_options.hh>
#include <def_use_scripting.hh>


// Required for the CFS own data types
#include "General/environment.hh"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace CoupledField {

  //! Define global instance of this class
  class ProgramOptions;
  extern ProgramOptions * progOpts;


  //! Base class for handling program options (command line / environment )

  //! This class is the base class for all classes that implement the handling
  //! of parameters specified on the command line. The class specifies a set
  //! ofmethods that define the public query methods a command
  //! line handler for NACS must provide. These are also the only methods
  //! which other NACS classes should call on a command line handler. This
  //! will allow for easy exchangebility of a concrete implementation of the
  //! command line handler concept. Another part of the concept is that the
  //! command line parameters are made available to all parts of NACS by
  //! a global pointer to an instance of this object called ::commandLine .
  //!
  //!
  //! The class requires that for each allowed command line parameter (with
  //! the exception of --help) a corresponding query method must be provided
  //! in the derived class.
  //!
  //! The following table gives an overview of the
  //! allowed parameters, their types and the associated query method.
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="5" align="center">
  //!         <b>Overview of NACS command line parameters</b>
  //!       </td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center"><b>command line parameter</b></td>
  //!       <td align="center"><b>parameter type</b></td>
  //!       <td align="center"><b>description</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>query method</b></td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">\<simName\></td>
  //!       <td align="center">string</td>
  //!       <td><em>name of the simulation run</em><br>
  //!       <td align="center">---</td>
  //!       <td align="center">%GetSimName()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-m / --meshFile = \<meshFileName\></td>
  //!       <td align="center">string</td>
  //!       <td><em>name of mesh file</em></td>
  //!       <td align="center">\<simName\>.mesh</td>
  //!       <td align="center">%GetMeshFile()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-s / --schemaPath = \<schemaPath\></td>
  //!       <td align="center">string</td>
  //!       <td><em>path to the XML schema file</em></td>
  //!       <td align="center">(assigned during compilation)</td>
  //!       <td align="center">%GetSchemaPath()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-p / --paramFile = \<xmlFileName\></td>
  //!       <td align="center">string</td>
  //!       <td><em>name of XML parameter file</em></td>
  //!       <td align="center">\<simName\>.xml</td>
  //!       <td align="center">%GetParamFile()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-w / --writeSkeleton</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>write skeleton of XML file for subsequent
  //!           simulation</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetWriteSkeleton()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-t / --traceDepth = \<traceDepth\></td>
  //!       <td align="center">non-negative integer</td>
  //!       <td><em>depth of function tracing</em></td>
  //!       <td align="center">0</td>
  //!       <td align="center">%GetTraceDepth()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-g / --printGrid</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>only read grid from mesh-file and write it to
  //!           output file</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetPrintGrid()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-e / --showEqnMap</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>print the dof to equation number mapping to OLAS report
  //!           file</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetShowEqnMap()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-d / --doProfile</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>do use system calls for generating profiling
  //!           information</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetDoProfile()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-h / --help</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>print usage information</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">---</td>
  //!     </tr>
  //!   </table>
  //! </center>
  class ProgramOptions
  {

  public:

    //! Standard constructor
    ProgramOptions( Integer argc,
                    const char **argv );

    //! Virtual destructor

    //! The destructor is made virtual
    //! for inheritance purposes.
    virtual ~ProgramOptions();
      
    //! Gather information from commandline and environment

    //! This method triggers the reading of information from the command line
    //! and the environment.
    void ParseData();
    
    // =======================================================================
    // QUERY METHODS FOR PARAMETERS
    // =======================================================================

    //@{
    //! \name Query methods for command line parameters
    //! For each command line parameter there exists an associated query
    //! method. If no command line value was specified than the query method
    //! will return the corresponding default, if it exists.

    //! Return flag for print help
    virtual bool GetPrintHelp();

    //! Return base name of simulation run (without path information)

    //! This method can be used to query the name of the current simulation
    //! run, which may be something like e.g. plate3D. This name is used by
    //! NACS internally as basename for the generation of several (default)
    //! names for input and output files.
    virtual std::string GetSimName() const;
    
    //! Return path to simulation files 
    
    //! This method can be used to query the directory of the current,
    //! i.e. where the output files will be created.
    virtual fs::path GetSimPath() const;
    virtual std::string GetSimPathStr() const;

    //! Return name of XML parameter file (including path)

    //! This method can be used to query the name of the parameter file in
    //! XML format that contains the steering parameters for the simulation.
    virtual fs::path GetParamFile() const;
    virtual std::string GetParamFileStr() const;

#ifdef USE_SCRIPTING
    //! Return (optional) name of scripting file

    //! This method returns the full name of an (optional) scripting file,
    //! e.g. a TCL-script that contains routines provided by the user to
    //! interact dynamically with NACS. If no filename was provided an
    //! empty string is returned.
    virtual fs::path GetScriptFile() const;
    virtual std::string GetScriptFileStr() const;
#endif    

    //! Return path to XML schema file

    //! This method can be used to query the path to the XML schema file
    //! used by validating XML parsers to verify the formal correctness
    //! of the XML parameter file.
    //! \note
    //! - There is currently no way to specify the name of the schema
    //!   file itself. This must be called NACS.xsd!
    //! - This path is also used to locate the default XML-file that is
    //!   currently still needed by the XMLParamHandler.
    virtual fs::path GetSchemaPath() const;
    virtual std::string GetSchemaPathStr() const;

    //! Return name of mesh file (including path)

    //! This method can be used to query the name of the mesh file containing
    //! the description of the FEM mesh for the simulation.
    virtual fs::path GetMeshFile() const;
    virtual std::string GetMeshFileStr() const;

    //! Return printGrid flag

    //! This method can be used to query the status of the printGrid flag.
    //! By specifying this flag one instructs the executable to do not
    //! perform an actual simulation, but to only import the grid and
    //! re-export it to an output file in the format specified in the XML
    //! parameter file.
    virtual bool GetPrintGrid() const;

    //! Return doProfile flag

    //! This method can be used to query the status of the doProfile flag.
    //! By specifying this flag one instructs the executable to generate
    //! profiling information. By default system calls are used (under Linux)
    //! to obtain information on things like e.g. memory footprint of the
    //! simulation run.
    //! \note The flag is only of use in the case that profiling was enabled
    //!       during compilation by defining the PROFILING macro.
    virtual bool GetDoProfile() const;


    //! Return Restart flag

    //! This method can be used to query the status of the restart flag.
    //! If this flag is true the simulation restarts from an previous state.
    virtual bool GetRestart() const;

    //! Return writeSkeleton flag

    //! This method can be used to query the status of the writeSkeleton flag.
    //! As a convenience for the NACS user it is possible to let the
    //! executable write a skeleton XML parameter file that must then be
    //! filled out by the user for a subsequent simulation run.
    virtual bool GetWriteSkeleton() const;

    //! Return depth of function tracing

    //! This method can be used to query the depth desired for generating
    //! function trace information.
    //! \note Function tracing is a time-consuming process that can easily
    //!       lead to the generation of a tremendously large output file
    //!       containing the trace information. Thus, tracing must explicitely
    //!       be compiled into the executable by turning on the respective
    //!       option in <em>Makefile.option</em>. Specifying a trace depth of
    //!       zero avoids generation of the trace file, but does not remove
    //!       the work associated with function tracing.
    virtual UInt GetTraceDepth() const;

    //! Return getDumpStats flag

    //! This method can be used to query the status of the getDumpStats flag.
    //! Starting NACS with the dump stats flags lets NACS print out
    //! valuable information about the executable like the distro on which
    //! it was built on, the compiler it was built with and so on.
    //! The execution stops after the printout.
    virtual bool GetDumpStats() const;
    
    //! Returns license path
    
#ifdef DEBUG
    //! Return forceSegfault flag

    //! This method can be used to query the status of the forceSegFault flag.
    //! If this flag is set now exception will be thrown, but a segmentation 
    //! fault will be forced instead, which enables one to use a debugger
    //! to get a stack trace.
    virtual bool GetForceSegFault() const;
#endif
    //@}
    

    // =======================================================================
    // AUXILLIARY METHODS FOR OUTPUTTING INFORMATION
    // =======================================================================

    //@{ \name Auxilliary methods for outputting information

    //! Print help information to command line
    void PrintHelp( std::ostream &out );
    
    //! Print a summary of the values of the command line parameters
    
    //! Calling this method will print a summary of the values of the command
    //! line parameters to the specified output stream. The values include
    //! the parameters supplied on the command line and the default values for
    //! those parameters where no value was supplied.
    void PrintParams( std::ostream &out, bool colorise = true );
    
    //! Print information about NACSexecutable

    //! Calling this method will print a number of valuable infos about
    //! about the executable to the standard output.
    void DumpStats( std::ostream & out, bool colorise = true );

    //! Helper method for DumpStats()
    void GetDumpString( std::ostream & outstr, bool colorise);
    // @}
    
 private:
    // =======================================================================
    // INTERNAL HELPER METHODS
    // =======================================================================

    //@{ \name Additional internal helper methods
    
    //! Helper function for mapping environment variables to internal names
    std::string EnvironmentNameMapper( std::string envVarName );
    
    //@}

    // =======================================================================
    // INTERNAL DATA
    // =======================================================================

    //! Command line arguments as vector
    std::vector<std::string> args_;

    //! Help message as string
    std::string helpMsg_;

    //! Boost's argument map
    po::variables_map varMap_;

  };

}
#endif
