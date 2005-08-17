#ifndef FILE_BASE_COMMANDLINE_HANDLER
#define FILE_BASE_COMMANDLINE_HANDLER


// Required for the CFS own data types
#include "General/environment.hh"
#include "Utils/StdVector.hh"


namespace CoupledField {


  //! Base class for handling parameters specified on the command line

  //! This class is the base class for all classes that implement the handling
  //! of parameters specified on the command line. The class specifies a set
  //! of pure virtual methods that define the public query methods a command
  //! line handler for CFS++ must provide. These are also the only methods
  //! which other CFS classes should call on a command line handler. This
  //! will allow for easy exchangebility of a concrete implementation of the
  //! command line handler concept. Another part of the concept is that the
  //! command line parameters are made available to all parts of CFS++ by
  //! a global pointer to an instance of this object called ::commandLine .
  //!
  //! Besides this, the class defines for each parameter a short and a long
  //! version for setting it on the command line (e.g. -h and --help). These
  //! are denoted as markers and stored as strings (e.g. markerHelp_ and
  //! markerLongHelp_). Derived classes should only work using theses strings,
  //! so that the definition of the marker is restricted to this base class.
  //! In addition the class defines for for each parameter a short description
  //! in string format (e.g. helpHelp_) that is used for generating the usage
  //! information and may also be employed in message generation in a derived
  //! class.
  //!
  //! The class requires that for each allowed command line parameter (with
  //! the exception of --help) a corresponding query method must be provided
  //! in the derived class.
  //! The assignment of default values for parameters that are not explicitely
  //! assigned on the command line is also done in this base class via
  //! methods like e.g. DefaultMeshFile(). Derived classes should again make
  //! use of these methods for assigning default values.
  //!
  //! The following table gives an overview of the
  //! allowed parameters, their types and the associated query method.
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="5" align="center">
  //!         <b>Overview of CFS++ command line parameters</b>
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
  //!       <td align="center">%GetTrace()</td>
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
  //!       <td align="center">-n / --noProfile</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>do not use system call for generating profiling
  //!           information</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetNoProfile()</td>
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
  class BaseCommandLineHandler
  {

  public:

    //! Virtual destructor

    //! The destructor is made virtual
    //! for inheritance purposes.
    virtual ~BaseCommandLineHandler() {
      ENTER_FCN( "BaseCommandLineHandler::~BaseCommandLineHandler" );
    };

    // =======================================================================
    // QUERY METHODS FOR PARAMETERS
    // =======================================================================

    //@{
    //! \name Query methods for command line parameters
    //! For each command line parameter there exists an associated query
    //! method. If no command line value was specified than the query method
    //! will return the corresponding default, if it exists.

    //! Return base name of simulation run

    //! This method can be used to query the name of the current simulation
    //! run, which may be something like e.g. plate3D. This name is used by
    //! CFS++ internally as basename for the generation of several (default)
    //! names for input and output files.
    virtual std::string GetSimName() const = 0;

    //! Return name of XML parameter file

    //! This method can be used to query the name of the parameter file in
    //! XML format that contains the steering parameters for the simulation.
    virtual std::string GetParamFile() const = 0;

    //! Return path to XML schema file

    //! This method can be used to query the path to the XML schema file
    //! used by validating XML parsers to verify the formal correctness
    //! of the XML parameter file.
    //! \note
    //! - There is currently no way to specify the name of the schema
    //!   file itself. This must be called CFS.xsd!
    //! - This path is also used to locate the default XML-file that is
    //!   currently still needed by the XMLParamHandler.
    virtual std::string GetSchemaPath() const = 0;

    //! Return name of mesh file

    //! This method can be used to query the name of the mesh file containing
    //! the description of the FEM mesh for the simulation.
    virtual std::string GetMeshFile() const = 0;

    //! Return printGrid flag

    //! This method can be used to query the status of the printGrid flag.
    //! By specifying this flag one instructs the executable to do not
    //! perform an actual simulation, but to only import the grid and
    //! re-export it to an output file in the format specified in the XML
    //! parameter file.
    virtual Boolean GetPrintGrid() const = 0;

    //! Return showEqnMap flag

    //! This method can be used to query the status of the showEqnMap flag.
    //! By specifying this flag one instructs the executable to report the
    //! maps that relates global to region-local node numbers and node number
    //! to equation number in the algebraic system. The maps are written to
    //! the standard OLAS report file.
    virtual Boolean GetShowEqnMap() const = 0;

    //! Return noProfile flag

    //! This method can be used to query the status of the noProfile flag.
    //! By specifying this flag one instructs the executable to not generate
    //! profiling information. By default system call are used (under Linux)
    //! to obtain information on things like e.g. memory footprint of the
    //! simulation run.
    //! \note The flag is only of use in the case that profiling was enabled
    //!       during compilation by defining the PROFILING macro.
    virtual Boolean GetNoProfile() const = 0;

    //! Return writeSkeleton flag

    //! This method can be used to query the status of the writeSkeleton flag.
    //! As a convenience for the CFS++ user it is possible to let the
    //! executable write a skeleton XML parameter file that must then be
    //! filled out by the user for a subsequent simulation run.
    virtual Boolean GetWriteSkeleton() const = 0;

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
    virtual UInt GetTraceDepth() const = 0;

    //@}


    // =======================================================================
    // STATIC STRING MEMBERS
    // =======================================================================

    //! \name Strings containing explanation of command line parameters
    //! These strings contain for each command line parameter a brief
    //! description of its meaning. They are employed e.g. in the PrintUsage()
    //! method.
    //@{
    const static std::string helpSimName_;
    const static std::string helpParamFile_;
    const static std::string helpMeshFile_;
    const static std::string helpTraceDepth_;
    const static std::string helpWriteSkeleton_;
    const static std::string helpPrintGrid_;
    const static std::string helpHelp_;
    const static std::string helpSchemaPath_;
    const static std::string helpShowEqnMap_;
    const static std::string helpNoProfile_;
    //@}

    //! \name Strings containing short markers for command line parameters
    //! These strings specify the command line parameters by defining their
    //! short form textual representation
    //@{
    const static std::string markerParamFile_;
    const static std::string markerMeshFile_;
    const static std::string markerTraceDepth_;
    const static std::string markerWriteSkeleton_;
    const static std::string markerPrintGrid_;
    const static std::string markerHelp_;
    const static std::string markerSchemaPath_;
    const static std::string markerShowEqnMap_;
    const static std::string markerNoProfile_;
    //@}

    //! \name Strings containing long markers for command line parameters
    //! These strings specify the command line parameters by defining their
    //! long form textual representation
    //@{
    const static std::string markerLongParamFile_;
    const static std::string markerLongMeshFile_;
    const static std::string markerLongTraceDepth_;
    const static std::string markerLongWriteSkeleton_;
    const static std::string markerLongPrintGrid_;
    const static std::string markerLongHelp_;
    const static std::string markerLongSchemaPath_;
    const static std::string markerLongShowEqnMap_;
    const static std::string markerLongNoProfile_;
    //@}


    // =======================================================================
    // ASSIGNMENT OF DEFAULT VALUES
    // =======================================================================

    //@{
    //! \name Setting default values for command line parameters
    //!       These methods provide default values for the allowed command
    //!       line parameters. Derived classes should use calls to these
    //!       methods for handling of defaults in their implementation of the
    //!       query methods instead of providing their own defaults.

    //! Returns default value for --traceDepth parameter

    //! This method returns default value for --traceDepth parameter. The
    //! current default is 0, which indicates no function tracing.
    //! \return 0
    UInt DefaultTraceDepth() const {
      return 0;
    }

    //! Returns default value for --paramFile parameter

    //! This method returns default value for --paramFile parameter. The
    //! current default is to compose the name of the parameter file, by
    //! adding a <em>.xml</em> postfix to the name of the current simulation
    //! run
    //! \return \<simName\>.xml
    std::string DefaultParamFile() const {
      return GetSimName() + ".xml";
    }

    //! Returns default value for --meshFile parameter

    //! This method returns default value for --meshFile parameter. The
    //! current default is to compose the name of the parameter file, by
    //! adding a <em>.mesh</em> postfix to the name of the current simulation
    //! run
    //! \return \<simName\>.mesh
    std::string DefaultMeshFile() const {
      return GetSimName() + ".mesh";
    }

    //! Returns default value for --printGrid parameter

    //! This method returns default value for --printGrid parameter. The
    //! current default is to perform a full simulation and to not only
    //! print the grid.
    //! \return FALSE
    Boolean DefaultPrintGrid() const {
      return FALSE;
    }

    //! Returns default value for --showEqnMap parameter

    //! This method returns the default value for the --showEqnMap parameter.
    //! The current default is not to report the node number <-> equation
    //! number mapping.
    //! \return FALSE
    Boolean DefaultShowEqnMap() const {
      return FALSE;
    }

    //! Returns default value for --writeSkeleton parameter

    //! This method returns default value for --writeSkeleton parameter. The
    //! current default is to perform a full simulation and to not only
    //! write a skeleton XML-file.
    //! \return FALSE
    Boolean DefaultWriteSkeleton() const {
      return FALSE;
    }

    //! Returns default value for --schemaPath parameter

    //! This method returns default value for --schemaPath parameter. The
    //! default name and location of the XML schema file is hard-coded
    //! during compile time making use of the \<XMLSCHEMA\> compile macro.
    //! \return \<XMLSCHEMA\>
    std::string DefaultSchemaPath() const {
      return XMLSCHEMA;
    }

    //! Returns default value for --noProfile parameter

    //! This method returns default value for --noProfile parameter. The
    //! default is to generate profile information, if profiling was enabled
    //! during compile time.
    //! \return - FALSE if PROFILING was enabled during compilation
    //!         - TRUE otherwise
    Boolean DefaultNoProfile() const {
#ifdef PROFILING
      return FALSE;
#else
      return TRUE;
#endif
    }

    //@}


    // =======================================================================
    // AUXILLIARY METHODS FOR OUTPUTTING INFORMATION
    // =======================================================================

    //@{ \name Auxilliary methods for outputting information

    //! Print summary of possible command line parameters

    //! Calling this method will print a summary of the possible command line
    //! parameters to the standard output.
    void PrintUsage();

    //! Print a summary of the values of the command line parameters

    //! Calling this method will print a summary of the values of the command
    //! line parameters to the specified output stream. The values include
    //! the parameters supplied on the command line and the default values for
    //! those parameters where no value was supplied.
    void PrintParams( std::ostream &out, Boolean colorise = TRUE );

    //@}
  };

}

#endif
