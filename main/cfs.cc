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

#ifdef MpCCI
#include <cci.h>
#endif

#ifdef GRIDLIB
#include <Domain/Gridlib/interface_gridlib.hh>
#endif

#include <AlgebraicSystem/abstractAlgSys.hh>
#include <Driver/transientdriver.hh>
#include <Driver/staticdriver.hh>
#include <Driver/harmonicDriver.hh>

#include <Domain/GridCFS/interface_gridcfs.hh>

#ifdef NETGEN
#include <Domain/NetGen/interface_netgen.hh>
#endif

using namespace CoupledField;

Integer main(int argc, char *argv[])
{
  std::cout << std::endl;
  std::cout << " Welcome to CFS++ session. " << std::endl << std::endl;
 
  if (argc < 2) 
    {
      std::cout << " \033[36mUsage\033[0m : cfs name "<< std::endl 
		<< "\t \033[36m name \033[0m: name of input file without extension" 
		<< std::endl << std::endl;
      Error("Invalid running of cfs. See Usage above.");
    }

  // always write info file: material parameters, data of nonlin iteration ...
  InfoPrint=TRUE;

  Char * name=argv[argc-1];

#ifdef MpCCI
  CCI_Init(&argc, &argv);  
#endif  

  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  FileType * ptInputfile=ptDefineFiles->Create_ptFileType();

  WriteResults * ptOut=ptDefineFiles->Create_ptWriteResults();

  //	TimeFunc * ptTimeFunc=NULL;
  TimeFunc * ptTimeFunc = new TimeFunc(ptInputfile);

  Domain * domain=new Domain(ptInputfile, ptOut, ptTimeFunc);

  //choose your driver
  BaseDriver * ptdriver;  
  std::string analysis;
  conf->get("analysis", analysis);

  if (analysis=="static") 
    ptdriver = new StaticDriver(domain);
  else if (analysis=="transient") 
    ptdriver = new TransientDriver(domain);
  else if (analysis=="harmonic")
    ptdriver = new HarmonicDriver(domain);
  else if (analysis=="harmonic")
    ptdriver = new HarmonicDriver(domain);
  else
    Error("Driver not supported",__FILE__,__LINE__);

  //solve your problem
  std::string adaptTimeOn, adaptSpaceOn;

  adaptSpaceOn =  "no";
  conf->ifget("adapttime",adaptTimeOn);
  adaptTimeOn = "no";
  conf->ifget("adaptspace",adaptSpaceOn);

  if (adaptTimeOn == "yes")  
    ptdriver->SolveProblemAdapt();
  else
    if 
      (adaptSpaceOn == "yes") ptdriver->SolveProblemAdaptSpace();
    else 
      ptdriver->SolveProblem();

  oClockTotal.ClockCount(MyClock::end,"Total time");

#ifdef MpCCI
    CCI_Finalize();
#endif

  //delete objects
  if (ptdriver) delete ptdriver;
  if (ptTimeFunc) delete ptTimeFunc;
  //  if (domain) delete domain;
  if (ptDefineFiles) delete ptDefineFiles; // it should be deleted the last

  return 0;
}
