#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <iomanip>
#include <stdarg.h>

#include <general_head.hh> 
#include <utils_head.hh>
#include <datainout_head.hh>
#include <elements_head.hh>
#include <forms_head.hh>
#include <linalg_head.hh>   
#include <domain_head.hh>
#include <pde_head.hh> 
#include <driver_head.hh>


using namespace CoupledField;

void main(int argc, char *argv[])
{

  std::cout << " Wellcome to sample session. " << argc << std::endl ;
  std::cout << " \033[36mUsage\033[0m : cfs -ext [-i] name [-m materialfile]"<< std::endl 
       << "\t \033[36m ext \033[0m: format of input file( implemented: dat ) " << std::endl
       << "\t \033[36m i \033[0m: to create info-file " << std::endl
       << "\t \033[36m name \033[0m: name of input file without extension" << std::endl << endl;

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

  Clock oClockTotal;
  oClockTotal.ClockCount(Clock::beg);

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
#ifdef Grid
    Grid<Point2D> * ptGridlib=new InterfaceGridlib<Point2D>(ptInputfile);
    std::cout << " Test " << std::endl;
    std::cout << " max number of nodes " << std::endl;
    std::cout << ptGridlib->GetMaxnumnodes(0) << std::endl;
    std::cout << " max num of elements " << std::endl;
    std::cout << ptGridlib->GetMaxnumElem(0) << std::endl;
#endif

//  Driver * ptDriver=new Driver(ptInputfile,1,materialdata);
//  ptDriver->SolveNewmarkMethod(ptUnverg);

  oClockTotal.ClockCount(Clock::end,"Total time");

/// Putzen

  if (ptInputfile) delete ptInputfile;
  if (ptUnverg) delete ptUnverg;
//  if (ptDriver) delete ptDriver;
  if (materialdata) delete materialdata;
}
