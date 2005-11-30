#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

// Be a little bit "colorful"
#ifndef SUPPRESS_COLORED_OUTPUT
// #define COLOR_INIT "\033[1m"
#define COLOR_INIT "\033[34m"
#define COLOR_STOP "\033[0m"
#else
#define COLOR_INIT ""
#define COLOR_STOP ""
#endif

namespace CoupledField {

  // ------------------------------------------------
  //  Definition of static data members of the class
  // ------------------------------------------------

  // Help strings
  const std::string BaseCommandLineHandler::helpSimName_       =
  "name of simulation run";

  const std::string BaseCommandLineHandler::helpParamFile_     =
  "name of XML parameter file for the simulation";

  const std::string BaseCommandLineHandler::helpSchemaPath_    =
  "path to directory containing the XML schema file";

  const std::string BaseCommandLineHandler::helpMeshFile_      =
  "name of mesh file for the simulation";  

  const std::string BaseCommandLineHandler::helpScriptFileName_ =
  "(optional) name of script file to be evaluated";  

  const std::string BaseCommandLineHandler::helpTraceDepth_    =
#ifdef TRACE
  "depth of function tracing";
#else
  "depth of function tracing \033[1m(de-activated in this version)\033[0m";
#endif

  const std::string BaseCommandLineHandler::helpWriteSkeleton_ =
  "write skeleton of XML file for subsequent simulation";

  const std::string BaseCommandLineHandler::helpPrintGrid_     =
  "only read grid from mesh-file and write it to output file";

  const std::string BaseCommandLineHandler::helpShowEqnMap_    =
  "print node to equation number mapping to OLAS report file";

  const std::string BaseCommandLineHandler::helpDoProfile_     =
  "turns on generation of profiling information";

  const std::string BaseCommandLineHandler::helpRestart_     =
  "read restart file";

  const std::string BaseCommandLineHandler::helpHelp_          =
  "print this usage information";

  // Short forms of markers
  const std::string BaseCommandLineHandler::markerParamFile_       = "-p";
  const std::string BaseCommandLineHandler::markerMeshFile_        = "-m";
  const std::string BaseCommandLineHandler::markerScriptFileName_  = "-e";
  const std::string BaseCommandLineHandler::markerTraceDepth_      = "-t";
  const std::string BaseCommandLineHandler::markerWriteSkeleton_   = "-w";
  const std::string BaseCommandLineHandler::markerPrintGrid_       = "-g";
  const std::string BaseCommandLineHandler::markerHelp_            = "-h";
  const std::string BaseCommandLineHandler::markerSchemaPath_      = "-s";
  const std::string BaseCommandLineHandler::markerDoProfile_       = "-d";
  const std::string BaseCommandLineHandler::markerRestart_         = "-r";

  // Long forms of markers
  const std::string BaseCommandLineHandler::markerLongParamFile_     =
  "--paramFile";
  const std::string BaseCommandLineHandler::markerLongMeshFile_      =
  "--meshFile";
  const std::string BaseCommandLineHandler::markerLongScriptFileName_ =
  "--scriptFile";
  const std::string BaseCommandLineHandler::markerLongTraceDepth_    =
  "--traceDepth";
  const std::string BaseCommandLineHandler::markerLongWriteSkeleton_ =
  "--writeSkeleton";
  const std::string BaseCommandLineHandler::markerLongPrintGrid_     =
  "--printGrid";
  const std::string BaseCommandLineHandler::markerLongHelp_          =
  "--help";
  const std::string BaseCommandLineHandler::markerLongSchemaPath_    =
  "--schemaPath";
  const std::string BaseCommandLineHandler::markerLongShowEqnMap_    =
  "--showEqnMap";
  const std::string BaseCommandLineHandler::markerLongDoProfile_     =
  "--doProfile";
  const std::string BaseCommandLineHandler::markerLongRestart_     =
  "--restart";


  // ---------------------------------------------------------------
  //  Implementation of some class methods too large for the header
  // ---------------------------------------------------------------


  // **************
  //   PrintUsage
  // **************
  void BaseCommandLineHandler::PrintUsage() {

    ENTER_FCN( "BaseCommandLineHandler::PrintUsage" );

    std::ostream &os = std::cout;

    os << COLOR_INIT << "\n Usage: " << COLOR_STOP

       << "./cfs [param] <name of simulation run>\n\n"
       << " where [param] is one or more of the following parameters\n\n"

      // --paramFile
       << COLOR_INIT
       << " " << markerParamFile_ << ", " << markerLongParamFile_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpParamFile_ << "\n\n"

      // --schemaPath
       << COLOR_INIT
       << " " << markerSchemaPath_ << ", " << markerLongSchemaPath_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpSchemaPath_ << "\n\n"

      // --meshFile
       << COLOR_INIT
       << " " << markerMeshFile_ << ", " << markerLongMeshFile_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpMeshFile_ << "\n\n"

      // --scriptFile
       << COLOR_INIT
       << " " << markerScriptFileName_ << ", " << markerLongScriptFileName_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpScriptFileName_ << "\n\n"

      // --traceDepth
       << " " << COLOR_INIT
       << markerTraceDepth_ << ", " << markerLongTraceDepth_
       << COLOR_STOP
       << " = <non-negative integer>\n"
       << " " << helpTraceDepth_ << "\n\n"

      // --printGrid
       << " " << COLOR_INIT
       << markerPrintGrid_ << ", " << markerLongPrintGrid_
       << COLOR_STOP
       << '\n'
       << " " << helpPrintGrid_ << "\n\n"

      // --WriteSkeleton
       << " " << COLOR_INIT
       << markerWriteSkeleton_ << ", " << markerLongWriteSkeleton_
       << COLOR_STOP
       << '\n'
       << " " << helpWriteSkeleton_ << "\n\n"

      // --showEqnMap
       << " " << COLOR_INIT
       << " "  << ", " << markerLongShowEqnMap_
       << COLOR_STOP
       << '\n'
       << " " << helpShowEqnMap_ << "\n\n"

      // --doProfile
       << " " << COLOR_INIT
       << markerDoProfile_ << ", " << markerLongDoProfile_
       << COLOR_STOP
       << '\n'
       << " " << helpDoProfile_ << "\n\n"

      // --Restart
       << " " << COLOR_INIT
       << markerRestart_ << ", " << markerLongRestart_
       << COLOR_STOP
       << '\n'
       << " " << helpRestart_ << "\n\n"

      // --help
       << " " << COLOR_INIT
       << markerHelp_ << ", " << markerLongHelp_
       << COLOR_STOP
       << '\n'
       << " " << helpHelp_ << "\n\n";

  }


  // ***************
  //   PrintParams
  // ***************
  void BaseCommandLineHandler::PrintParams( std::ostream &out,
                                            Boolean colorise ) {

    ENTER_FCN( "BaseCommandLineHandler::PrintParams" );

    std::string colorInit = "";
    std::string colorStop = "";

    if ( colorise == TRUE ) {
      colorInit = COLOR_INIT;
      colorStop = COLOR_STOP;
    }

    out << colorInit
        << "\n\nValues of command line parameters (including defaults):\n\n"
        << colorStop

        << " name of simulation run = "
        << colorInit
        << GetSimName()
        << colorStop << '\n'

        << ' ' << markerLongParamFile_ << " = "
        << colorInit
        << GetParamFile()
        << colorStop << '\n'

        << ' ' << markerLongSchemaPath_ << " = "
        << colorInit
        << GetSchemaPath()
        << colorStop << '\n'

        << ' ' << markerLongMeshFile_ << " = "
        << colorInit
        << GetMeshFile()
        << colorStop << '\n'

        << ' ' << markerLongScriptFileName_ << " = "
        << colorInit
        << GetScriptFileName()
        << colorStop << '\n'

        << ' ' << markerLongTraceDepth_ << " = "
        << colorInit
        << GetTraceDepth()
        << colorStop << '\n'

        << ' ' << markerLongPrintGrid_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetPrintGrid() == TRUE )
        << colorStop << '\n'

        << ' ' << markerLongShowEqnMap_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetShowEqnMap() == TRUE )
        << colorStop << '\n'

        << ' ' << markerLongDoProfile_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetDoProfile() == TRUE )
        << colorStop << '\n'

        << ' ' << markerLongRestart_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetRestart() == TRUE )
        << colorStop << '\n'

        << ' ' << markerLongWriteSkeleton_ << " = "
        << colorInit
        << std::boolalpha
        << ( GetWriteSkeleton() == TRUE )
        << colorStop << "\n\n";
  }

}
