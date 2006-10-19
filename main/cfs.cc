#include <iostream>
#include <iomanip>
#include "Utils/myclock.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/CommandLine/CommandLineHandlerSetting.hh"
#include "Driver/driver_header.hh"
#include "Domain/domain.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "Utils/tracing.hh"
#include "Utils/profiler.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef MpCCI
#include <cci.h>
#endif

#ifdef GRIDLIB
#include "Domain/Gridlib/interface_gridlib.hh"
#endif


#ifdef NETGEN
#include "Domain/NetGen/interface_netgen.hh"
#endif


#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif


using namespace CoupledField;

#ifdef PARALLEL
#define STDOUT if (!commrank) std::cout
#else
#define STDOUT std::cout
#endif

int main( int argc, const char **argv ) {


  // =========================================================================
  // INITIALISATION OF MPI
  // =========================================================================

#ifdef PARALLEL //initialize MPI
  int commrank,commsize;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc,(char***)(&argv));
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);
#endif //parallel

  
  // =========================================================================
  // TIMING
  // =========================================================================
  MyClock oClockTotal;
  oClockTotal.Start();


  // =========================================================================
  // CREATE INFORMATION OBJECT
  // =========================================================================

  // Create object for logging information
  Info = new WriteInfo();

  // Log program startup
  Info->PrintHeader();


  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================
  commandLine = new CommandLineHandlerSetting( argc, argv );

  // =========================================================================
  // GENERATE OBJECT FOR HANDLING FILE-IO
  // =========================================================================
  Info->StartProgress( "Generating object for handling file-IO" );
  DefineInOutFiles FileHandler( commandLine->GetSimName().c_str() );
  Info->FinishProgress();
  

  // =========================================================================
  // GENERATE CENTRAL MESSENGER OBJECT
  // =========================================================================
#ifdef USE_SCRIPTING

  // Try to determine (optional) script file name
  std::string scriptFileName = commandLine->GetScriptFileName();
  
  messenger = FileHandler.CreateScriptMessenger( scriptFileName );
  
  // Check if script file was provided
  if ( scriptFileName != "" ) {
    std::stringstream msg;
    msg << "Activating Scripting using '" << scriptFileName << "'";
    Info->StartProgress( msg.str() );

    // Create new central messenger object (up to now only tcl available)
    //messenger = new TCL_CFSMessenger();
    //messenger = new PY_CFSMessenger();
    messenger->ReadScriptFile(scriptFileName );
    Info->FinishProgress();    
  }

#endif


  // =========================================================================
  // ACTIVATE DEBUGGING STUFF
  // =========================================================================

  // Activate function tracing
#ifdef TRACE
  Integer traceDepth = commandLine->GetTraceDepth();
  OutInfo::FcnTraceHandler::SetMaxTraceDepth( traceDepth );
  if ( traceDepth == 0 ) {
    Info->StartProgress( "De-activating function tracing since depth = 0" );
    Info->FinishProgress();
  }
  else {
    std::stringstream msg;
    msg << "Activating function tracing with depth = " << traceDepth;
    Info->StartProgress( msg.str() );
    FileHandler.OpenFile( TRACE_FILE );
    Info->FinishProgress();
  }
#endif

  // Open file for debugging ouput
#ifdef DEBUG
  Info->StartProgress( "Opening file for writing debugging output" );
  FileHandler.OpenFile( DEBUG_FILE );
  Info->FinishProgress();
#endif

#ifdef MEMTRACE
  Info->StartProgress( "Opening file for tracing memory allocation" );
  FileHandler.OpenFile( MEMTRACE_FILE );
  Info->FinishProgress();
#endif

#ifdef PROFILING
  if ( commandLine->GetDoProfile() == true ) {
    Info->StartProgress( "Opening file for profiling output" );
  }
  else {
    Info->StartProgress( "Skipping generation of profiling output" );
  }
  profiler = new Profiler();
  SETPROFILE( "Begin of program" );
  Info->FinishProgress();
#endif

  // =========================================================================
  // WRITE SKELETON XML-FILE
  // =========================================================================

  // This block is used for writing a skeleton XML-File to make life easier
  // for the user. This also will import some information from the mesh-file.
  // Since normally all opening of files is handled in definefiles
  // the debug and trace file are opened separately here. Definefiles cannot
  // be used, since it would try to open other files also, that need not be
  // present????
  if ( commandLine->GetWriteSkeleton() == true ) {
#ifdef TRACE
    FileHandler.OpenFile( TRACE_FILE );
#endif
 
#ifdef DEBUG
    FileHandler.OpenFile( DEBUG_FILE );
#endif

    // class writing log-information
    Info = new WriteInfo();
    Info->PrintHeader();
    SkeletonConf *ptskel = new SkeletonConf();
    ptskel->WriteConf();
    delete ptskel;

    return 0;
  }

  // =========================================================================
  // SETUP OF MPCCI
  // =========================================================================

#ifdef MpCCI
  Info->StartProgress( "Setting up MpCCI interface" );
  CCI_Init( &argc, const_cast<char***>(&argv) );  
  Info->FinishProgress();
#endif


  // =========================================================================
  // HANDLE XML-FILE
  // =========================================================================

  // Generate parameter handler and pass address to global pointer
  std::string xmlFile = commandLine->GetParamFile();

  // Write information to command line
  std::stringstream msg;
  msg << "Reading parameters from file '" << xmlFile << "'";
  Info->StartProgress( msg.str() );
  
#ifdef USE_XERCES
  std::string cfsSchema = commandLine->GetSchemaPath();
  cfsSchema += "/CFS-Simulation/CFS.xsd";
  
  std::string defaults = commandLine->GetSchemaPath();
  defaults += "/CFS-Simulation/Defaults/CFS++Defaults.xml";
  params = new XMLParamHandler( xmlFile.c_str(), cfsSchema.c_str(), 
                                defaults.c_str() );
#else
  params = new PlainXMLParamHandler( xmlFile );
#endif

  Info->FinishProgress();


  // =========================================================================
  // SETUP OF IO-STUFF
  // =========================================================================

  // read type of mesh-libraray
  std::string libmesh;
  params->Get( "meshLibrary", libmesh, "input" );

  FileType *ptInputfile;
  if ( libmesh == "structGrid" ) {
    //for struct grid we do not need a mesh-input-file
    ptInputfile = NULL;
  }
  else {
    // Generate mesh reader
    Info->StartProgress( "Generating mesh reader" );
    ptInputfile = FileHandler.CreateMeshFileHandler();
    Info->FinishProgress();
  }

  // generate material handler
  MaterialHandler * ptMatHandler;
  Info->StartProgress( "Generating material reader");
  ptMatHandler = FileHandler.CreateMaterialHandler();
  Info->FinishProgress();
  
  Info->StartProgress( "Generating remaining output files" );

  // Open file for status reports by CFS++
  Info->CreateFile();

  // Log command line parameters
  std::ostream *myInfo = Info->GetInfoStreamPointer();
  commandLine->PrintParams( *myInfo, false );

  // Open file for status reports by OLAS
  FileHandler.OpenFile( OLAS_FILE );

  // Generate Write results object
  WriteResults *ptOut = FileHandler.Create_ptWriteResults();

  Info->FinishProgress();


  // =========================================================================
  // GENERATION OF DOMAIN OBJECT
  // =========================================================================
  TimeFunc myTimeFunc;

  SETPROFILE("Before Creation of Domain");
  domain = new  Domain( ptInputfile, ptOut, &myTimeFunc, ptMatHandler );
  SETPROFILE("After Creation of Domain");


  // =========================================================================
  // Only output of grid
  // =========================================================================
  if ( commandLine->GetPrintGrid() == true ) {
    STDOUT << "Printing grid to file " << myEndl << myEndl;
    PrintGridOnly = true;
    domain->PrintGrid();
    delete domain;
    delete Info;
    delete commandLine; 
    return 0;
  }


  // =========================================================================
  // SELECTION OF DRIVER
  // =========================================================================

#ifdef ADAPTGRID
  flags = new Flags();
#endif

  BaseDriver   *ptdriver = NULL;
  AnalysisType analysisType;
  std::string  analysis;
  bool      adaptspace;

  // Determine type of analysis and, if adaptivity is to be used
  params->Get( "type", analysis, "analysis" );
  adaptspace = params->IsSet( "adaptspace" );
  String2Enum( analysis, analysisType );


  // Generate log message
  Info->StartProgress( "Generating driver" );

  // Generate driver
  switch( analysisType ) {


  case STATIC:

    // with adaptivity
    if ( adaptspace ) {
#ifdef ADAPTGRID
      ptdriver = new StaticAdaptSpaceDriver(*domain);
#else
      std::string errmsg;
      errmsg  = "Your version of cfs does not support adaptivity! ";
      errmsg += "Recompile with Adaptivity = yes";
      Info->Error( errmsg, __FILE__, __LINE__ );
#endif
    }

    // without adaptivity
    else {
      ptdriver = new StaticDriver( domain );
    }
    break;

  case TRANSIENT:
    ptdriver = new TransientDriver( domain );
    break;

  case HARMONIC:
    // calls Driver for parameter identification, using harmonic analysis
    if ( analysis == "paramIdent" ) {
      ptdriver = new piezoParamIdent( domain );
    }
    else
      ptdriver = new HarmonicDriver( domain );
    break;

  case TRANSIENTHARMONIC:
    ptdriver = new TransientHarmonicDriver( domain );
    break;

  case EIGENFREQUENCY:
    ptdriver = new EigenFrequencyDriver( domain );
    break;

    //case MULTIHARMONIC:
    //ptdriver = new MultiHarmonicDriver( domain );
    //break;

  case MULTI_SEQUENCE:
    ptdriver = new MultiSequenceDriver( domain );
    break;

  default:
    (*error) << "Driver '" << analysis << "' not supported!";
    Error( __FILE__, __LINE__ );
  }

  // Pass driver to domain object
  domain->SetDriver( ptdriver );

  Info->FinishProgress();

  // =========================================================================
  // SOLUTION PHASE
  // =========================================================================

  ptdriver->SolveProblem();
  std::cout << std::endl;
  Info->StartProgress( "Finished solving the problem" );
  Info->FinishProgress();

  std::cout << std::endl;
  std::cout << std::endl;

#ifdef USE_SCRIPTING

    // Call intialization procedure
    StdVector<std::string> context;
    context.Push_back( commandLine->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Finish, context);
#endif

  Double wTime, cTime;
  oClockTotal.GetTime(wTime, cTime);
  std::cout << std::setw(70) << std::setfill('=') << "" << std::endl;
  std::cout << "TOTAL TIME:\t\t Wall clock: " << wTime
            << " s\t CPU: " << cTime << " s" << std::endl;
  std::cout << std::setw(70) << std::setfill('=') << "" 
            << std::endl << std::endl;


  // =========================================================================
  // CLEANUP PHASE
  // =========================================================================

#ifdef MpCCI
  CCI_Finalize();
#endif

  // Delete objects
  delete ptdriver;
  delete domain;
  delete Info;
  delete params;
  delete commandLine;

#ifdef USE_SCRIPTING
  delete messenger;
  messenger = NULL;
#endif

#ifdef PROFILING
  delete profiler;
  profiler = NULL;
#endif

  // As the last thing we close the trace file (if exists)
#ifdef TRACE
  delete trace;
  trace = NULL;
#endif

  // Cleanup parallel stuff if it exists
#ifdef PARALLEL
  MPI_Finalize();
#endif

  // Delete global string streams
  delete error;
  delete warning;

  // Seems that everything went fine
  return 0;
}
