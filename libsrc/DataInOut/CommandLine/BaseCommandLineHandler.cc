#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

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
