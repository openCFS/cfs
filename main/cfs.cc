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

int main(int argc, char *argv[])
{

  Integer numargs=2;  
  Boolean SkeletonPrint=FALSE;

  if (argc > 1)
    for (Integer i=1; i<argc; i++)
      {
	if (!strcmp("-skel", argv[i])) 
	  {
	    SkeletonPrint=TRUE;
	    numargs++;
	  }
	else if (!strcmp("-grid", argv[i])) 
	  {
	    PrintGridOnly = TRUE;
	    numargs++;
	  }
#ifdef TRACE
	else if (!strcmp("-notrace", argv[i])) 
	  {

	    OutInfo::FcnTraceHandler::SetMaxTraceDepth(0);
	    numargs++;
	  }
#endif
	
      }
  
#ifdef PARALLEL //initialize MPI
  int commrank,commsize;
  	
  MPI_Comm comm = MPI_COMM_WORLD;

	MPI_Init(&argc,&argv);
          
    MPI_Comm_size(comm, &commsize);
    MPI_Comm_rank(comm, &commrank);
#endif //parallel

//  OutInfo::trace = new std::ofstream("trace.data");
  
  if (argc < numargs) 
    {
      STDOUT << std::endl;
      STDOUT << " \033[36mUsage\033[0m : cfs [-options] name "<< std::endl 
	     << "\t \033[36m name    \033[0m: name of input file without extension" 
	     << std::endl 
	     << "\t \033[36m options \033[0m: -skel for writing a skeleton of a config-file" << std::endl 
	     << "\t            -grid for writing just the grid to result-file" << std::endl 
#ifdef DEBUG
	     << "\t            -notrace prevent writing a trace file" << std::endl 
#endif
	     << std::endl ;
      std::cerr << std::endl << "\033[31mERROR:\033[0m " << "Invalid Running of CFS" << std::endl;
      exit(1);
    }
  
  char *name = argv[argc-1];
  
  Char * filename=new char[100];

  // for writing a skeleton of a config file by using the information
  // from the mesh-file
  if (SkeletonPrint==TRUE)
    {
#ifdef TRACE
      strcpy(filename, name);
      trace=new std::ofstream(strcat(filename,".trace"));
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


#ifdef MpCCI
  CCI_Init(&argc, &argv);  
#endif  

  
  // class writing log-information
  Info = new WriteInfo();
  Info->PrintHeader();

  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);
  Info->CreateFile(name);
  
   if (PrintGridOnly)
    STDOUT << "Printing grid to file " << name << ".unv" << myEndl << myEndl;

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  
  FileType * ptInputfile=ptDefineFiles->Create_ptFileType();

  WriteResults * ptOut=ptDefineFiles->Create_ptWriteResults(ptInputfile);
  
  
  // TimeFunc * ptTimeFunc=NULL;
  TimeFunc * ptTimeFunc = new TimeFunc(ptInputfile);

  Domain * domain=new Domain(ptInputfile, ptOut, ptTimeFunc);


  // choose your driver
  BaseDriver       * ptdriver;  
  std::string      analysis, errMsg;
  Boolean          adaptspace;
  AnalysisType analysisType;

#ifndef XMLPARAMS
  conf->get("analysis", analysis);
  adaptspace=conf->get_option("adaptspace");
#else
  params->Get( "type", analysis, "analysis" );
  adaptspace = params->IsSet( "adaptspace" );
#endif

  Info->StartProgress("Creating driver");
  String2Enum(analysis, analysisType);

  switch(analysisType)
    {
    case STATIC:
      if (adaptspace)
	{
#ifdef ADAPTGRID
	  ptdriver = new StaticAdaptSpaceDriver(domain);
#else
	  std::string errmsg;
	  errmsg  = "Your version of cfs does not support adaptivity! ";
	  errmsg += "Recompile with Adaptivity = yes";
	  Info->Error( errmsg, __FILE__, __LINE__ );
#endif
	}
      else
	ptdriver = new StaticDriver(domain);
      break;
    case TRANSIENT:
      ptdriver = new TransientDriver(domain);
      break;
    case HARMONIC:
       if (analysis == "paramIdent"){ // calls Driver for parameter identification, using harmonic analysis
          ptdriver = new piezoParamIdent(domain);
	  std::cout<<"calls driver for parameter identification"<< std::endl;
       }
       else
          ptdriver = new HarmonicDriver(domain);
      break;
    case MULTI_SEQUENCE:
      ptdriver = new MultiSequenceDriver(domain);
      break;
    default:
      errMsg = "Driver '";
      errMsg += analysis;
      errMsg += "' not supported!";
      Info->Error( errMsg.c_str(),  __FILE__, __LINE__ );
    }
  Info->FinishProgress();
  ptdriver->SolveProblem();
  std::cout << std::endl;
  Info->StartProgress("Finished solving the problem");
  Info->FinishProgress();

  std::cout << std::endl;
  std::cout << std::endl;
  oClockTotal.ClockCount(MyClock::end,"Total time");
 
#ifdef MpCCI
  CCI_Finalize();
#endif
 
  // Delete objects
  if (ptdriver) 
    delete ptdriver;
  if (ptTimeFunc) 
    delete ptTimeFunc;
  if (Info)
    delete Info;

  // ptDefineFiles must be destroyed as the last of all objects!
  // The reason is that if TRACE is defined all methods want to
  // log into (*trace) which is deleted by the destructor of the
  // DefineInOutFiles class!
  if (ptDefineFiles) 
    delete ptDefineFiles;

#ifdef PARALLEL
MPI_Finalize();
#endif

  return 0;

}



