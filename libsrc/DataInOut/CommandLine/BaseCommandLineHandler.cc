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

  const std::string BaseCommandLineHandler::helpMeshFile_      =
  "name of mesh file for the simulation";  

  const std::string BaseCommandLineHandler::helpTraceDepth_    =
  "depth of function tracing";

  const std::string BaseCommandLineHandler::helpWriteSkeleton_ =
  "write skeleton of XML file for subsequent simulation";

  const std::string BaseCommandLineHandler::helpPrintGrid_     =
  "only read grid from mesh-file and write it to output file";

  const std::string BaseCommandLineHandler::helpHelp_          =
  "print this usage information";

  // Short forms of markers
  const std::string BaseCommandLineHandler::markerParamFile_     = "-p";
  const std::string BaseCommandLineHandler::markerMeshFile_      = "-m";
  const std::string BaseCommandLineHandler::markerTraceDepth_    = "-t";
  const std::string BaseCommandLineHandler::markerWriteSkeleton_ = "-s";
  const std::string BaseCommandLineHandler::markerPrintGrid_     = "-g";
  const std::string BaseCommandLineHandler::markerHelp_          = "-h";

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
  "--printGrid";
  const std::string BaseCommandLineHandler::markerLongHelp_          =
  "--help";


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

      // --meshFile
       << COLOR_INIT
       << " " << markerMeshFile_ << ", " << markerLongMeshFile_
       << COLOR_STOP
       << " = <string>\n"
       << " " << helpMeshFile_ << "\n\n"

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

      // --help
       << " " << COLOR_INIT
       << markerHelp_ << ", " << markerLongHelp_
       << COLOR_STOP
       << '\n'
       << " " << helpHelp_ << "\n\n";

  };


  // ***************
  //   PrintParams
  // ***************
  void BaseCommandLineHandler::PrintParams( std::ostream &out ) {

    ENTER_FCN( "BaseCommandLineHandler::PrintParams" );

    out << COLOR_INIT
        << "\n Status of command line parameters (including defaults):\n\n"
        << COLOR_STOP

        << " name of simulation run = "
        << COLOR_INIT
        << GetSimName()
        << COLOR_STOP << '\n'

        << ' ' << markerLongParamFile_ << " = "
        << COLOR_INIT
        << GetParamFile()
        << COLOR_STOP << '\n'

        << ' ' << markerLongMeshFile_ << " = "
        << COLOR_INIT
        << GetMeshFile()
        << COLOR_STOP << '\n'

        << ' ' << markerLongTraceDepth_ << " = "
        << COLOR_INIT
        << GetTraceDepth()
        << COLOR_STOP << '\n'

        << ' ' << markerLongPrintGrid_ << " = "
        << COLOR_INIT
        << GetPrintGrid()
        << COLOR_STOP << '\n'

        << ' ' << markerLongWriteSkeleton_ << " = "
        << COLOR_INIT
        << GetWriteSkeleton()
        << COLOR_STOP << "\n\n";
  };

}
