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
#include <DataInOut/material.hh>
#include <DataInOut/timefunc.hh>

#ifdef GRIDLIB
#include <Domain/Gridlib/interface_gridlib.hh>
#endif

#include <AlgebraicSystem/abstractAlgSys.hh>
#include <Driver/transientdriver.hh>
#include <Driver/staticdriver.hh>

#include <Domain/GridCFS/interface_gridcfs.hh>

#ifdef NETGEN
#include <Domain/NetGen/interface_netgen.hh>
#endif

using namespace CoupledField;

Integer main(int argc, char *argv[])
{
  std::cout << " Wellcome to sample session. " << argc << std::endl ;
  std::cout << " \033[36mUsage\033[0m : cfs [-i] name "<< std::endl 
	    << "\t \033[36m i \033[0m: to create info-file " << std::endl
	    << "\t \033[36m name \033[0m: name of input file without extension" << std::endl << std::endl;

  if (argc < 2) Error("Invalid running of cfs. See Usage above.");

  if (!strcmp("-i", argv[1])) InfoPrint=TRUE;

  Char * name=argv[argc-1];
  
  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  Material * ptMaterial=NULL;
  std::string material;
  conf->get("material_file",material);
  if (material != "non") ptMaterial=new Material(material.c_str());
	
  FileType * ptInputfile=ptDefineFiles->Create_ptFileType();

  WriteResults * ptOut=ptDefineFiles->Create_ptWriteResults();

  //	TimeFunc * ptTimeFunc=NULL;
  TimeFunc * ptTimeFunc = new TimeFunc(ptInputfile);

  Domain * domain=new Domain(ptInputfile,ptOut,ptMaterial, ptTimeFunc);

  //	domain->PrintGrid(0);
  //choose your driver
  BaseDriver * ptdriver;  
  std::string analysis;
  conf->get("analysis", analysis);

  if (analysis=="static") 
    ptdriver = new StaticDriver(domain);
  else  
    {
    ptdriver = new TransientDriver(domain);
    //	ptTimeFunc = new TimeFunc(ptInputfile);
    }

  //solve your problem
  std::string adaptTimeOn, adaptSpaceOn;
  //  conf->get("adapttime",adaptTimeOn);
  conf->get("adaptspace",adaptSpaceOn);
  adaptTimeOn = "no";

  if (adaptTimeOn == "yes")  
    ptdriver->SolveProblemAdapt();
  else
    if 
      (adaptSpaceOn == "yes") ptdriver->SolveProblemAdaptSpace();
    else 
      ptdriver->SolveProblem();

  oClockTotal.ClockCount(MyClock::end,"Total time");

  /// Putzen
    if (ptdriver) delete ptdriver;
    if (ptMaterial) delete ptMaterial;
    if (ptTimeFunc) delete ptTimeFunc;
    if (domain) delete domain;
    if (ptDefineFiles) delete ptDefineFiles; // it should be deleted the last

    return 1;
}
