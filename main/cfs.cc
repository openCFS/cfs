#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>
#include <list>

#include "myclock.hh"
#include "definefiles.hh"
#include "material.hh"
#include "timefunc.hh"
#include "acousticPDE.hh"
//#include "interface_gridlib.hh"
#include "driver.hh"
#include <abstractAlgSys.hh>
//#include <interface_piles.hh>
#include "transientdriver.hh"

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
  Material * materialdata=NULL;
  if (!strcmp("-m", argv[argc-2])) 
    { 
      name=argv[argc-3];
      materialdata=new Material(argv[argc-1]);
    }
  else name=argv[argc-1];

  DefineInOutFiles oDefFiles(name);

  MyClock oClockTotal;
  oClockTotal.ClockCount(MyClock::beg);

  FileType * ptInputfile=oDefFiles.Create_ptFileType(argv[1]);
  WriteResults<Point2D> * ptOut=oDefFiles.Create_ptWriteResults2d();

  TimeFunc * ptTimeFunc=new TimeFunc(ptInputfile);

/*
  /// Test 3D mesh

  Integer data[1];
  ptInputfile-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);
  std::cout << "Number of nodes" << data[0] << std::endl;
 
  Point3D * ptCoord=new Point3D[data[0]];
  ptInputfile->ReadCoordinate(ptCoord, data[0]);
  std::cout << "We have read coordinates" << std::endl;

  Grid<Point3D> * ptGridlib=new InterfaceGridlib<Point3D>(ptInputfile);
  ptGridlib->Read();
  
  ptOut->Create(ptGridlib,0);

*/
//    Driver<Point2D> * ptDriver=new Driver<Point2D>(ptInputfile,1,materialdata);
//    ptDriver->SolveNewmarkMethod(ptOut);

    Domain<Point2D> *domain=new Domain<Point2D>(ptInputfile,ptOut,materialdata, ptTimeFunc);

  // print grid to unverg-file
  domain->PrintGrid(0);

  //choose your driver
  BaseDriver *ptdriver = new TransientDriver(domain);

  //solve your problem
  ptdriver->SolveProblemAdapt();

#ifdef TestGMRes
  Double eps=1e-15;
  LinSystem<Double, Matrix<Double> > * ptLS=new  LinSystem<Double, Matrix<Double> >(eps);
  ptLS->Set();
  if (ptLS->BiCGSTAB(100,Jacobi)) std::cout << "convergence is true";
  ptLS->printxscreen();
  ptLS->check();
  // if (ptLS->GMRes_m(100, Jacobi, 100)) std::cout << "convergence is true";
#endif

  oClockTotal.ClockCount(MyClock::end,"Total time");

/// Putzen

  if (ptInputfile) delete ptInputfile;
  if (ptOut) delete ptOut;
//  if (ptGridlib) delete ptGridlib;
//  if (ptDriver) delete ptDriver;
  if (materialdata) delete materialdata;
}
