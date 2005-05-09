#ifndef FILE_BASE_COMMANDLINE_HANDLER
#define FILE_BASE_COMMANDLINE_HANDLER


// Required for the CFS own data types
#include "General/environment.hh"
#include "Utils/StdVector.hh"



namespace CoupledField
{


  //! Base class for handling parameters specified on the command line

  //! This class is the base class for all classes that implement the handling
  //! of parameters specified on the command line. It is a pure virtual class
  //! which defines the public query methods a command line handler for CFS++
  //! must provide. These are also the only methods which other CFS classes
  //! should call on a command line handler.
  //!
  //!
  //! The class requires that for each allowed command line parameter a
  //! corresponding query method must be provided. The following table
  //! gives an overview of the allowed parameters and the associated
  //! query method.
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
  //!       <td align="center">-p / --paramFile = \<xmlFileName\></td>
  //!       <td align="center">string</td>
  //!       <td><em>name of XML parameter file</em></td>
  //!       <td align="center">\<simName\>.xml</td>
  //!       <td align="center">%GetParamFile()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-s / --writeSkeleton</td>
  //!       <td align="center">flag / boolean</td>
  //!       <td><em>write skeleton of XML file for subsequent
  //!           simulation</em></td>
  //!       <td align="center">false = not set</td>
  //!       <td align="center">%GetWriteSkeleton()</td>
  //!     </tr>
  //!
  //!     <tr>
  //!       <td align="center">-t / --trace = \<traceDepth\></td>
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
    virtual Integer GetTraceDepth() const = 0;

    //@}


    // =======================================================================
    // STATIC STRING MEMBERS
    // =======================================================================

    //@{
    //! \name Strings containing textual explanation of command line
    //!       parameters
    const static std::string helpSimName_;
    const static std::string helpParamFile_;
    const static std::string helpMeshFile_;
    const static std::string helpTraceDepth_;
    const static std::string helpWriteSkeleton_;
    const static std::string helpPrintGrid_;
    //@}

    //@{
    //! \name Strings containing textual markers for command line parameters
    const static std::string markerParamFile_;
    const static std::string markerMeshFile_;
    const static std::string markerTraceDepth_;
    const static std::string markerWriteSkeleton_;
    const static std::string markerPrintGrid_;

    const static std::string markerLongParamFile_;
    const static std::string markerLongMeshFile_;
    const static std::string markerLongTraceDepth_;
    const static std::string markerLongWriteSkeleton_;
    const static std::string markerLongPrintGrid_;
    //@}

  };


  // ------------------------------------------------
  //  Definition of static data members of the class
  // ------------------------------------------------

  // Help strings
  const std::string BaseCommandLineHandler::helpSimName_       =
  "name of simulation run";

  const std::string BaseCommandLineHandler::helpParamFile_     =
  "name of XML parameter file for the simulation";

  const std::string BaseCommandLineHandler::helpMeshFile_      =
  "name of mesh file for the simulation";  

  const std::string BaseCommandLineHandler::helpTraceDepth_    =
  "depth of function tracing";

  const std::string BaseCommandLineHandler::helpWriteSkeleton_ =
  "write skeleton of XML file for subsequent simulation";

  const std::string BaseCommandLineHandler::helpPrintGrid_     =
  "only read grid from mesh-file and write it to output file";

  // Short forms of markers
  const std::string BaseCommandLineHandler::markerParamFile_     = "-p";
  const std::string BaseCommandLineHandler::markerMeshFile_      = "-m";
  const std::string BaseCommandLineHandler::markerTraceDepth_    = "-t";
  const std::string BaseCommandLineHandler::markerWriteSkeleton_ = "-s";
  const std::string BaseCommandLineHandler::markerPrintGrid_     = "-g";

  // Long forms of markers
  const std::string BaseCommandLineHandler::markerLongParamFile_     =
  "--paramFile";
  const std::string BaseCommandLineHandler::markerLongMeshFile_      =
  "--meshFile";
  const std::string BaseCommandLineHandler::markerLongTraceDepth_    =
  "--traceDepth";
  const std::string BaseCommandLineHandler::markerLongWriteSkeleton_ =
  "--writeSkeleton";
  const std::string BaseCommandLineHandler::markerLongPrintGrid_     =
  "--PrintGrid";

}

#endif
