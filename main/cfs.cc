#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>
#include <math.h>

#include <vector>

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
#include "interface_netgen.hh"

using namespace CoupledField;

void main(int argc, char *argv[])
{
  std::cout << " Wellcome to sample session. " << argc << std::endl ;
  std::cout << " \033[36mUsage\033[0m : cfs [-i] name [-m materialfile]"<< std::endl 
       << "\t \033[36m i \033[0m: to create info-file " << std::endl
       << "\t \033[36m name \033[0m: name of input file without extension" << std::endl << std::endl;

  if (argc < 2) Error("Invalid running of scfe. See Usage above.");

  if (!strcmp("-i", argv[1])) InfoPrint=TRUE;

  Char * name=NULL;
  Material * ptMaterial=NULL;
  if (!strcmp("-m", argv[argc-2])) 
    { 
      name=argv[argc-3];
      ptMaterial=new Material(argv[argc-1]);
    }
  else name=argv[argc-1];

  DefineInOutFiles * ptDefineFiles=new DefineInOutFiles(name);

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  FileType * ptInputfile=ptDefineFiles->Create_ptFileType();

  WriteResults * ptOut=ptDefineFiles->Create_ptWriteResults2d();

  TimeFunc * ptTimeFunc=new TimeFunc(ptInputfile);

  Domain * domain=new Domain(ptInputfile,ptOut,ptMaterial, ptTimeFunc);

  // print grid to unverg-file
  domain->PrintGrid(0);

  //choose your driver
  BaseDriver *ptdriver = new TransientDriver(domain);

  //solve your problem
  std::string adaptTimeOn;
  conf->get("adapttime",adaptTimeOn,"Acoustic");
  if (adaptTimeOn == "yes")  ptdriver->SolveProblemAdapt();
                          else ptdriver->SolveProblem();

  oClockTotal.ClockCount(MyClock::end,"Total time");

/*
// for testing refinement of mesh
   WriteResults<Point2D> * ptOut=oDefFiles.Create_ptWriteResults2d();

   Grid * ptgrid=new InterfaceNetGen<Point2D>(ptInputfile);

   ptgrid->Read();
 
//   ptOut->Init(ptgrid);
//   ptOut->WriteGrid(0);

   ptgrid->SubdivideUniform(0);
   ptOut->Init(ptgrid);
   ptOut->WriteGrid(0);

   if (ptgrid) delete ptgrid;
*/
/*
   WriteResults<Point3D> * ptOut=oDefFiles.Create_ptWriteResults3d();

   Grid * ptgrid=new InterfaceNetGen<Point3D>(ptInputfile);

   ptgrid->Read();

   ptOut->Init(ptgrid);
   ptOut->WriteGrid(0);

//   ptgrid->SubdivideUniform(0);
//   ptOut->Init(ptgrid);
//   ptOut->WriteGrid(0);

   if (ptgrid) delete ptgrid;
   if (ptOut) delete ptOut;
*/
/*
   Grid * ptgrid=new GridInterfaceCFS<Point3D>(ptInputfile);
   ptgrid->Read();

   WriteResults<Point3D> * ptOut=new WriteResultsUnverg<Point3D>("3d");

   ptOut->Init(ptgrid);
   ptOut->WriteGrid(0);
*/

/// Putzen
  if (ptdriver) delete ptdriver;
  if (ptMaterial) delete ptMaterial;
//  if (ptDefineFiles) delete ptDefineFiles;

}
