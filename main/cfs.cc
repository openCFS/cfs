#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>

#include "myclock.hh"
#include "definefiles.hh"
#include "material.hh"
#include "timefunc.hh"
#include "acousticPDE.hh"
//#include "interface_gridlib.hh"
#include "driver.hh"
#include <abstractAlgSys.hh>
//#include <interface_piles.hh>
#include <basedriver.hh>

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

  TimeFunc * ptTimeFunc=new TimeFunc(ptInputfile);

  WriteResults<Point2D> * ptUnverg=new WriteResultsUnverg<Point2D>(name); 

/*
  Integer data[1];
  ptInputfile-> ReadGeneralAnalChoice(data, FileType::numnode, FileType::endGAnal);
  std::cout << "Number of nodes" << data[0] << std::endl;
 
  Point3D * ptCoord=new Point3D[data[0]];
  ptInputfile->ReadCoordinate(ptCoord, data[0]);
  std::cout << "We have read coordinates" << std::endl;

  Grid<Point3D> * ptGridlib=new InterfaceGridlib<Point3D>(ptInputfile);
  ptGridlib->Read();
  
  ptUnverg->Create(ptGridlib,0);
*/
    Driver<Point2D> * ptDriver=new Driver<Point2D>(ptInputfile,1,materialdata);
    ptDriver->SolveNewmarkMethod(ptUnverg);

/*
 //  //! choose your grid class
  Grid<Point2D> *grid =new GridInterfaceCFS<Point2D>(ptInputfile);
  
  // construct domain
  Domain<Point2D> *domain=new Domain<Point2D>(ptInputfile,ptUnverg,materialdata,grid);

  // print grid to unverg-file
  domain->PrintGrid(0);

  // setup the pde-information
  domain->InitPDE();

  //choose your algebraic system
  AbstractAlgebraicSys * algsys = new AlgSysPILES();
  
  // setup algebraic system
  domain->InitAlgSys(algsys);  

  //choose your driver
  BaseDriver *ptdriver = new BaseDriver(domain);

  //solve your problem
  ptdriver->SolveProblem();
*/

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
//  if (ptUnverg) delete ptUnverg;
//  if (ptGridlib) delete ptGridlib;
//  if (ptDriver) delete ptDriver;
  if (materialdata) delete materialdata;
}
