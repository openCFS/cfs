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
#include "interface_gridlib.hh"
#include "driver.hh"
#include <abstractAlgSys.hh>
#include <interface_piles.hh>
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

  OutResultUnverg * ptUnverg=new OutResultUnverg(name); 

  // DDD
  TimeFunc * ptTimeFunc=new TimeFunc(ptInputfile);
  // DDD

#ifdef TestGMRes
  Double eps=1e-15;
  LinSystem<Double, Matrix<Double> > * ptLS=new  LinSystem<Double, Matrix<Double> >(eps);
  ptLS->Set();
  if (ptLS->BiCGSTAB(100,Jacobi)) std::cout << "convergence is true";
  ptLS->printxscreen();
  ptLS->check();
  // if (ptLS->GMRes_m(100, Jacobi, 100)) std::cout << "convergence is true";
#endif
    Grid<Point2D> * ptGridlib=new InterfaceGridlib<Point2D>(ptInputfile);
    std::cout << " Test " << std::endl;
    std::cout << " max number of nodes " << std::endl;
    std::cout << ptGridlib->GetMaxnumnodes(0) << std::endl;
    std::cout << " max num of elements " << std::endl;
    std::cout << ptGridlib->GetMaxnumElem(0) << std::endl;
    std::cout << " number of nodes per element" << std::endl;
    std::cout <<  ptGridlib->GetNumNodesPerElem(0,0) << std::endl;
    std::cout << " connection of element " << std::endl;
    Integer * result=new Integer[3];

    ptGridlib->GetConnection(result, 0, 0, 3);
    for (int i=0; i < 3; i++)
     std::cout << result[i] << " ";
    std::cout << std::endl; 
    std::cout << "coordinates" << std::endl;  
    Point2D * ptCoord=new Point2D[3];
     ptGridlib->GetCoordOfNodesElem(0,0,3,ptCoord);  
    std::cout << " coordinates of element" << std::endl;
    for (int i=0; i<3; i++) std::cout << i <<" " << ptCoord[i].x << " " << ptCoord[i].y << std::endl;


    ptGridlib->GetCoordOfNodesElem(1,0,3,ptCoord);
    std::cout << " coordinates of element" << std::endl;
    for (int i=0; i<3; i++) std::cout << i <<" " << ptCoord[i].x << " " << ptCoord[i].y << std::endl;

    //    Driver * ptDriver=new Driver(ptInputfile,1,materialdata);
    //    ptDriver->SolveNewmarkMethod(ptUnverg);

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

  oClockTotal.ClockCount(MyClock::end,"Total time");

/// Putzen

  if (ptInputfile) delete ptInputfile;
  if (ptUnverg) delete ptUnverg;
//  if (ptGridlib) delete ptGridlib;
//  if (ptDriver) delete ptDriver;
  if (materialdata) delete materialdata;
}
