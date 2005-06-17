#include "Utils/myclock.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/CommandLine/CommandLineHandlerSetting.hh"
#include "Driver/driver_header.hh"
#include "Domain/domain.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "Utils/tracing.hh"


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
  oClockTotal.ClockCount(MyClock::beg);


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


  // =========================================================================
  // WRITE SKELETON XML-FILE
  // =========================================================================

  // This block is used for writing a skeleton XML-File to make life easier
  // for the user. This also will import some information from the mesh-file.
  // Since normally all opening of files is handled in definefiles
  // the debug and trace file are opened separately here. Definefiles cannot
  // be used, since it would try to open other files also, that need not be
  // present????
  if ( commandLine->GetWriteSkeleton() == TRUE ) {

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
  CCI_Init( &argc, &argv );  
  Info->FinishProgress();
#endif


  // =========================================================================
  // HANDLE XML-FILE
  // =========================================================================

  // Generate parameter handler and pass address to global pointer
  std::string xmlFile = commandLine->GetParamFile();
#ifdef USE_XERCES
  params = new XMLParamHandler( xmlFile.c_str() );
#else
  params = new PlainXMLParamHandler( xmlFile.c_str() );
#endif


  // =========================================================================
  // SETUP OF IO-STUFF
  // =========================================================================

  // Generate mesh reader
  Info->StartProgress( "Generating mesh reader" );
  FileType *ptInputfile = FileHandler.CreateMeshFileHandler();
  Info->FinishProgress();


  Info->StartProgress( "Generating remaining output files" );

  // Open file for status reports by CFS++
  Info->CreateFile();

  // Open file for status reports by OLAS
  FileHandler.OpenFile( OLAS_FILE );

  // Generate Write results object
  WriteResults *ptOut = FileHandler.Create_ptWriteResults();

  Info->FinishProgress();


  // =========================================================================
  // GENERATION OF DOMAIN OBJECT
  // =========================================================================
  TimeFunc myTimeFunc;
  Domain domain( ptInputfile, ptOut, &myTimeFunc );


  // =========================================================================
  // Only output of grid
  // =========================================================================
  if ( commandLine->GetPrintGrid() == TRUE ) {
    STDOUT << "Printing grid to file " <<  commandLine->GetSimName()
           << ".unv" << myEndl << myEndl;
    domain.PrintGrid();
    exit(0);
  }


  // =========================================================================
  // SELECTION OF DRIVER
  // =========================================================================

#ifdef ADAPTGRID
  flags = new Flags();
#endif

  BaseDriver   *ptdriver;  
  AnalysisType analysisType;
  std::string  analysis;
  Boolean      adaptspace;

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
      ptdriver = new StaticAdaptSpaceDriver(domain);
#else
      std::string errmsg;
      errmsg  = "Your version of cfs does not support adaptivity! ";
      errmsg += "Recompile with Adaptivity = yes";
      Info->Error( errmsg, __FILE__, __LINE__ );
#endif
    }

    // without adaptivity
    else {
      ptdriver = new StaticDriver( &domain );
    }
    break;

  case TRANSIENT:
    ptdriver = new TransientDriver( &domain );
    break;

  case HARMONIC:
    // calls Driver for parameter identification, using harmonic analysis
    if ( analysis == "paramIdent" ) {
      ptdriver = new piezoParamIdent( &domain );
    }
    else if ( analysis == "multiHarmonic" )
      ptdriver = new MultiHarmonicDriver( &domain );
    else
      ptdriver = new HarmonicDriver( &domain );
    break;

  case MULTI_SEQUENCE:
    ptdriver = new MultiSequenceDriver( &domain );
    break;

  case BUBBLEDYNAMIC:
    ptdriver = new BubbleDriver( &domain );
    break;

  default:
    (*error) << "Driver '" << analysis << "' not supported!";
    Error( __FILE__, __LINE__ );
  }

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
  oClockTotal.ClockCount(MyClock::end,"Total time");


  // =========================================================================
  // CLEANUP PHASE
  // =========================================================================

#ifdef MpCCI
  CCI_Finalize();
#endif

  // Delete objects
  delete ptdriver;
  // delete domain;
  // delete ptTimeFunc;
  delete Info;
  delete params;

  // As the last thing we close the trace file (if exists)
#ifdef TRACE
  delete trace;
  trace = NULL;
#endif

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;

}
