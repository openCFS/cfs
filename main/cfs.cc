#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
//#include <math.h>
#include <cmath>
#include <vector>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Utils/myclock.hh"
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
#include "DataInOut/SkeletonConf.hh"
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
      SkeletonConf *ptskel = new SkeletonConf(name);
      ptskel->WriteConf();
      delete ptskel;

      return 0;
    }


#ifdef MpCCI
  CCI_Init(&argc, &argv);  
#endif  

  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);
  // class writing log-information
  Info = new WriteInfo(name);
  Info->PrintHeader();


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
  std::string      analysis;
  Boolean          adaptspace;

#ifndef XMLPARAMS
  conf->get("analysis", analysis);
  adaptspace=conf->get_option("adaptspace");
#else
  params->Get( "type", analysis, "analysis" );
  adaptspace = params->IsSet( "adaptspace" );
#endif

  if (analysis=="static") 
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

  else if (analysis=="transient") 
    ptdriver = new TransientDriver(domain);

  else if (analysis=="harmonic")
    ptdriver = new HarmonicDriver(domain);
  else
    Info->Error( "Driver not supported", __FILE__, __LINE__ );

  ptdriver->SolveProblem();
  
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



