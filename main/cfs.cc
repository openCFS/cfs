#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <math.h>

#include "myclock.hh"
#include "definefiles.hh"
#include "material.hh"
#include "timefunc.hh"
#include "acousticPDE.hh"

#ifdef GRIDLIB
#include "interface_gridlib.hh"
#endif

#include "driver.hh"
#include <abstractAlgSys.hh>
//#include <interface_piles.hh>
#include "transientdriver.hh"
#include "staticdriver.hh"
#include "interface_gridcfs.hh"

//#include "interface_netgen.hh"

using namespace CoupledField;

void main(int argc, char *argv[])
{
  std::cout << " Wellcome to sample session. " << argc << std::endl ;
  std::cout << " \033[36mUsage\033[0m : cfs -ext [-i] name [-m materialfile]"<< std::endl 
       << "\t \033[36m ext \033[0m: format of input file( implemented: dat ) " << std::endl
       << "\t \033[36m i \033[0m: to create info-file " << std::endl
       << "\t \033[36m name \033[0m: name of input file without extension" << std::endl << std::endl;

  if (argc < 3) Error("Invalid running of scfe. See Usage above.");

  if (!strcmp("-i", argv[2])) InfoPrint=TRUE;

  Char * name=NULL;
  Material * ptMaterial=NULL;
  if (!strcmp("-m", argv[argc-2])) 
    { 
      name=argv[argc-3];
      ptMaterial=new Material(argv[argc-1]);
    }
  else name=argv[argc-1];

  DefineInOutFiles oDefFiles(name);

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  FileType * ptInputfile=oDefFiles.Create_ptFileType(argv[1]);

//  WriteResults<Point2D> * ptOut=oDefFiles.Create_ptWriteResults2d();
  WriteResults<Point3D> * ptOut=oDefFiles.Create_ptWriteResults3d();  

  TimeFunc * ptTimeFunc=new TimeFunc(ptInputfile);

//  Domain<Point2D> *domain=new Domain<Point2D>(ptInputfile,ptOut,ptMaterial, ptTimeFunc);
   Domain<Point3D> *domain=new Domain<Point3D>(ptInputfile,ptOut,ptMaterial, ptTimeFunc);

  // print grid to unverg-file
  domain->PrintGrid(0);

  //choose your driver
//  BaseDriver *ptdriver = new TransientDriver(domain);
   BaseDriver<Point3D> *ptdriver = new StaticDriver<Point3D>(domain);

  //solve your problem
  std::string adaptTimeOn;
  conf->get("adapttime",adaptTimeOn,"Acoustic");
  if (adaptTimeOn == "yes")  ptdriver->SolveProblemAdapt();
                          else ptdriver->SolveProblem();

  oClockTotal.ClockCount(MyClock::end,"Total time");

/*
   Grid<Point3D> * ptgrid=new GridInterfaceCFS<Point3D>(ptInputfile);
   ptgrid->Read();

   WriteResults<Point3D> * ptOut=new WriteResultsUnverg<Point3D>("3d");

   ptOut->Init(ptgrid);
   ptOut->WriteGrid(0);
*/

/// Putzen

  if (ptInputfile) delete ptInputfile;
  if (ptOut) delete ptOut;
//  if (ptdriver) delete ptdriver;
  if (ptMaterial) delete ptMaterial;
}
