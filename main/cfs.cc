#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <math.h>
#include <vector>

#include <Utils/myclock.hh>
#include <DataInOut/DefineFiles/definefiles.hh>
#include <DataInOut/timefunc.hh>
#include <DataInOut/WriteInfo.hh>

#ifdef MpCCI
#include <cci.h>
#endif

#ifdef GRIDLIB
#include <Domain/Gridlib/interface_gridlib.hh>
#endif

#include <AlgebraicSystem/abstractAlgSys.hh>
#include <Driver/driver_header.hh>
#include <Domain/domain.hh>
#include <DataInOut/SkeletonConf.hh>
#include <Domain/GridCFS/interface_gridcfs.hh>
#include <General/environment.hh>
//#include <Utils/storeSolution.hh>


#ifdef NETGEN
#include <Domain/NetGen/interface_netgen.hh>
#endif

using namespace CoupledField;

Integer main(int argc, char *argv[])
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
  
  if (argc < numargs) 
    {
      std::cout << std::endl;
      std::cout << " \033[36mUsage\033[0m : cfs [-options] name "<< std::endl 
		<< "\t \033[36m name    \033[0m: name of input file without extension" 
		<< std::endl 
		<< "\t \033[36m options \033[0m: -skel for writing a skeleton of a config-file" << std::endl 
		<< "\t            -grid for writing just the grid to result-file" << std::endl 
                << std::endl ;
      Error("Invalid running of cfs. See Usage above.");
    }
  
  Char * name=argv[argc-1];
  Char * filename=new Char[100];

  //for writing a skeleton of a config file by using the information from the mesh-file
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
    std::cout << "Printing grid to file " << name << ".unv" << myEndl << myEndl;

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  FileType * ptInputfile=ptDefineFiles->Create_ptFileType();

  WriteResults * ptOut=ptDefineFiles->Create_ptWriteResults(ptInputfile);
  

  //	TimeFunc * ptTimeFunc=NULL;
  TimeFunc * ptTimeFunc = new TimeFunc(ptInputfile);

  Domain * domain=new Domain(ptInputfile, ptOut, ptTimeFunc);


  //choose your driver
  BaseDriver       * ptdriver;  
  std::string      analysis;
  Boolean          adaptspace;
  conf->get("analysis", analysis);
  adaptspace=conf->get_option("adaptspace");

  if (analysis=="static") 
    if (adaptspace)   
#ifdef ADAPTGRID
      ptdriver = new StaticAdaptSpaceDriver(domain);
#else
      Error("Your version do not support adaptivity; recompile with Adaptivity = yes",__FILE__,__LINE__);
#endif
    else
      ptdriver = new StaticDriver(domain);

  else if (analysis=="transient") 
    ptdriver = new TransientDriver(domain);

  else if (analysis=="harmonic")
    ptdriver = new HarmonicDriver(domain);
  else
    Error("Driver not supported",__FILE__,__LINE__);

  ptdriver->SolveProblem();
  
  oClockTotal.ClockCount(MyClock::end,"Total time");

#ifdef MpCCI
  CCI_Finalize();
#endif

  //delete objects
  if (ptdriver) 
    delete ptdriver;
  if (ptTimeFunc) 
    delete ptTimeFunc;
  if (ptDefineFiles) 
    delete ptDefineFiles; 
  if (Info)
    delete Info;
  return 0;

}
