#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <vector>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Utils/myclock.hh"
#include "Utils/tracing.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/timefunc.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "ParamIdent/piezoParamIdent.hh"

#ifdef MpCCI
#include <cci.h>
#endif

#ifdef GRIDLIB
#include "Domain/Gridlib/interface_gridlib.hh"
#endif

#include "Driver/driver_header.hh"
#include "Domain/domain.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "Domain/GridCFS/interface_gridcfs.hh"
#include "General/environment.hh"

#ifdef NETGEN
#include "Domain/NetGen/interface_netgen.hh"
#endif

using namespace CoupledField;

#ifdef PARALLEL
#define STDOUT if (!commrank) std::cout
#else
#define STDOUT std::cout
#endif

int main(int argc, char *argv[]) {

  Integer numargs=2;
  Boolean SkeletonPrint=FALSE;

  // ==========================================================================
  // HANDLING OF COMMAND LINE PARAMETERS
  // ==========================================================================

  if (argc > 1) {
    for ( Integer i = 1; i < argc; i++ ) {
      if ( !strcmp( "-skel", argv[i]) ) {
	SkeletonPrint = TRUE;
	numargs++;
      }
      else if ( !strcmp( "-grid", argv[i]) ) {
	PrintGridOnly = TRUE;
	numargs++;
      }

#ifdef TRACE
      else if ( !strcmp("-trace", argv[i]) ) {
	OutInfo::FcnTraceHandler::SetMaxTraceDepth(TRACE);
	int dummy = TRACE;
	std::cerr << "Set tracing depth to" << dummy << std::endl;
	numargs++;
      }
#endif

      //---------------------------------------------
      //  Choice of bubbledynamical model
      else if ( !strcmp("-bubbletyp", argv[i]) ) {
	numargs += 2;
	std::string aux( argv[i+1] );
	//	String2Enum( aux, bubbleDyn );
      }
    }
    
  }

  if (argc < numargs) {
    STDOUT << std::endl;
    STDOUT << " \033[36mUsage\033[0m : cfs [-options] name \n"
	   << "\t\033[36m name    \033[0m:"
	   << " name of input file without extension\n"
	   << "\t\033[36m options \033[0m:"
	   << " -skel for writing a skeleton of a config-file\n"
	   << "\t           -grid for writing just the grid to result-file\n"
#ifdef DEBUG
	   << "\t           -trace write a trace file\n"
#endif
	   << "\n \033[31mERROR:\033[0m"
	   << " Invalid parameter specification. See usage details above.\n\n";
    exit(1);
  }

  // Obtain basename for CFS++ IO-files
  char *name = argv[argc-1];


  // ==========================================================================
  // INITIALISATION OF MPI
  // ==========================================================================

#ifdef PARALLEL //initialize MPI
  int commrank,commsize;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(comm, &commsize);
  MPI_Comm_rank(comm, &commrank);
#endif //parallel


  // ==========================================================================
  // SKELETON OF CONFIG-FILE
  // ==========================================================================

  // This block is used for writing a skeleton of a .conf or .xml parameter
  // file. Since normally all opening of files is handled in definefiles
  // the debug and trace file are opened separately here. Definefiles cannot
  // be used, since it would try to open other files also, that need not be
  // present????

  Char *filename = new char[100];

  // for writing a skeleton of a config file by using the information
  // from the mesh-file
  if (SkeletonPrint==TRUE)
    {
#ifdef TRACE
      strcpy(filename, name);
      trace = new std::ofstream(strcat(filename,".trace"));
      if (!trace) Error("Can't open trace-file");
#endif
 
#ifdef DEBUG
      strcpy(filename, name);
      debug=new std::ofstream(strcat(filename,".deb"));
      if (!debug) Error("Can't open debug-file");
#endif

      // class writing log-information
      Info = new WriteInfo();
      Info->PrintHeader();
      SkeletonConf *ptskel = new SkeletonConf(name);
      ptskel->WriteConf();
      delete ptskel;

      return 0;
    }


  // ==========================================================================
  // SETUP OF MPCCI
  // ==========================================================================

#ifdef MpCCI
  CCI_Init( &argc, &argv );  
#endif


  // ==========================================================================
  // SETUP OF IO-STUFF
  // ==========================================================================

  // Create object for logging information
  Info = new WriteInfo();

  // Log program startup
  Info->PrintHeader();

  // Open files for input/output
  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);
  Info->CreateFile(name);

  if (PrintGridOnly)
    STDOUT << "Printing grid to file " << name << ".unv" << myEndl << myEndl;

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);
  
  FileType *ptInputfile = ptDefineFiles->Create_ptFileType();

  WriteResults *ptOut = ptDefineFiles->Create_ptWriteResults(ptInputfile);
  
  
  // TimeFunc * ptTimeFunc=NULL;
  TimeFunc *ptTimeFunc = new TimeFunc(ptInputfile);

  Domain *domain = new Domain(ptInputfile, ptOut, ptTimeFunc);


  // ==========================================================================
  // SELECTION OF DRIVER
  // ==========================================================================

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
	ptdriver = new StaticDriver(domain);
      }
      break;

    case TRANSIENT:
      ptdriver = new TransientDriver(domain);
      break;

    case HARMONIC:
      // calls Driver for parameter identification, using harmonic analysis
      if (analysis == "paramIdent"){ 
	ptdriver = new piezoParamIdent(domain);
      }
      else if (analysis =="multiHarmonic")
  	ptdriver = new MultiHarmonicDriver(domain);
      else
	ptdriver = new HarmonicDriver(domain);
      break;

    case MULTI_SEQUENCE:
      ptdriver = new MultiSequenceDriver(domain);
      break;

    case BUBBLEDYNAMIC:
      ptdriver = new BubbleDriver(domain, name);
      break;

    default:
      {
	std::string errMsg = "Driver '";
	errMsg += analysis;
	errMsg += "' not supported!";
	Info->Error( errMsg.c_str(),  __FILE__, __LINE__ );
      }
    }

  Info->FinishProgress();


  // ==========================================================================
  // SOLUTION PHASE
  // ==========================================================================

  ptdriver->SolveProblem();
  std::cout << std::endl;
  Info->StartProgress( "Finished solving the problem" );
  Info->FinishProgress();

  std::cout << std::endl;
  std::cout << std::endl;
  oClockTotal.ClockCount(MyClock::end,"Total time");


  // ==========================================================================
  // CLEANUP PHASE
  // ==========================================================================

#ifdef MpCCI
  CCI_Finalize();
#endif

  // Delete objects
  delete ptdriver;
  delete domain;
  delete ptTimeFunc;
  delete Info;

  // ptDefineFiles must be destroyed as the last of all objects!
  // The reason is that if TRACE is defined all methods want to
  // log into (*trace) which is deleted by the destructor of the
  // DefineInOutFiles class!
  delete ptDefineFiles;

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;

}
