#ifndef FILE_COMMANDLINE_HANDLER_SETTING_HH
#define FILE_COMMANDLINE_HANDLER_SETTING_HH


#include "Settings.hh"
#include "BaseCommandLineHandler.hh"


namespace CoupledField {


  //! Command line handler built on Settings library by Uwe Fabricius
  class CommandLineHandlerSetting : public BaseCommandLineHandler {


  public:

    // =======================================================================
    // CONSTRUCTIION and DESTRUCTION
    // =======================================================================

    //@{
    //! \name Methods for construction, destruction and setup

    //! The constructor requires passing of the command line arguments
    CommandLineHandlerSetting( Integer argc, const char **argv );

    //! Destructor
    ~CommandLineHandlerSetting() {
      ENTER_FCN( "CommandLineHandlerSetting::~CommandLineHandlerSetting" );
    };

    //@}


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
    std::string GetSimName() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetSimName" );
      return simName_;
    };

    //! Return name of XML parameter file

    //! This method can be used to query the name of the parameter file in
    //! XML format that contains the steering parameters for the simulation.
    std::string GetParamFile() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetParamFile" );
      std::string paramFile = DefaultParamFile();
      Setting *aux = commandLine_.getSetting( markerParamFile_.c_str() );
      if ( aux != NULL ) {
        paramFile = aux->getString();
      }
      return paramFile;
    };

    //! Return path to XML schema file

    //! This method can be used to query the path to the XML schema file
    //! used by validating XML parsers to verify the formal correctness
    //! of the XML parameter file.
    //! \note
    //! - There is currently no way to specify the name of the schema
    //!   file itself. This must be called CFS.xsd!
    //! - This path is also used to locate the default XML-file that is
    //!   currently still needed by the XMLParamHandler.
    std::string GetSchemaPath() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetSchemaPath" );
      std::string schemaPath = DefaultSchemaPath();
      Setting *aux = commandLine_.getSetting( markerSchemaPath_.c_str() );
      if ( aux != NULL ) {
        schemaPath = aux->getString();
      }
      return schemaPath;
    };

    //! Return name of mesh file

    //! This method can be used to query the name of the mesh file containing
    //! the description of the FEM mesh for the simulation.
    std::string GetMeshFile() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetMeshFile" );
      std::string meshFile = DefaultMeshFile();
      Setting *aux = commandLine_.getSetting( markerMeshFile_.c_str() );
      if ( aux != NULL ) {
        meshFile = aux->getString();
      }
      return meshFile;
    };

    //! Return printGrid flag

    //! This method can be used to query the status of the printGrid flag.
    //! By specifying this flag one instructs the executable to do not
    //! perform an actual simulation, but to only import the grid and
    //! re-export it to an output file in the format specified in the XML
    //! parameter file.
    Boolean GetPrintGrid() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetPrintGrid" );
      Boolean retVal = DefaultPrintGrid();
      Setting *aux = commandLine_.getSetting( markerPrintGrid_.c_str() );
      if ( aux != NULL ) {
        retVal = TRUE;
      }
      return retVal;
    };

    //! Return showEqnMap flag

    //! This method can be used to query the status of the shwoEqnMap flag.
    //! By specifying this flag one instructs the executable to report the
    //! maps that relates global to region-local node numbers and node
    //! numbers to equation numbers in the algebraic system. The maps are
    //! written to the standard OLAS report file.
    Boolean GetShowEqnMap() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetShowEqnMap" );
      Boolean retVal = DefaultShowEqnMap();
      Setting *aux = commandLine_.getSetting( markerShowEqnMap_.c_str() );
      if ( aux != NULL ) {
        retVal = TRUE;
      }
      return retVal;
    };

    //! Return noProfile flag

    //! This method can be used to query the status of the noProfile flag.
    //! By specifying this flag one instructs the executable to not generate
    //! profiling information. By default system call are used (under Linux)
    //! to obtain information on things like e.g. memory footprint of the
    //! simulation run.
    //! \note The flag is only of use in the case that profiling was enabled
    //!       during compilation by defining the PROFILING macro.
    Boolean GetNoProfile() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetNoProfile" );
      Boolean retVal = DefaultNoProfile();
      Setting *aux = commandLine_.getSetting( markerNoProfile_.c_str() );
      if ( aux != NULL ) {
        retVal = TRUE;
      }
      return retVal;
    };

    //! Return writeSkeleton flag

    //! This method can be used to query the status of the writeSkeleton flag.
    //! As a convenience for the CFS++ user it is possible to let the
    //! executable write a skeleton XML parameter file that must then be
    //! filled out by the user for a subsequent simulation run.
    Boolean GetWriteSkeleton() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetWriteSkeleton" );
      Boolean retVal = DefaultWriteSkeleton();
      Setting *aux = commandLine_.getSetting( markerWriteSkeleton_.c_str() );
      if ( aux != NULL ) {
        retVal = TRUE;
      }
      return retVal;
    };

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
    UInt GetTraceDepth() const {
      ENTER_FCN( "CommandLineHandlerSetting::GetTraceDepth" );
      Setting *aux = commandLine_.getSetting( markerTraceDepth_.c_str() );
      UInt retVal = DefaultTraceDepth();
      if ( aux != NULL ) {
        retVal = (UInt)(aux->getInt());
      }
      return retVal;
    };

    //@}

  private:

    //! The SettingDataBase object that does the actual command line handling

    //! The %CommandLineHandlerSettings class does not directly handle the
    //! command line. Instead this work is delegated to an instance of the
    //! SettingDataBase class from the Settings library by Uwe Fabricius.
    SettingDataBase commandLine_;

    //! This string stores the name of the current simulation run

    //! This attribute stores as string the name of the current simulation
    //! run as it was specified on the command line.
    std::string simName_;

  };

}

#endif
