#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <string>
#include <time.h>
#include <stdio.h>
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

  cout << " Wellcome to sample session. " << argc << endl ;
  cout << " Usage: cfs -ext [-i] name [-m materialfile]" << endl 
       << "\t ext: format of input file( implemented: dat ) " << endl
       << "\t i: to create info-file " << endl
       << "\t name: name of input file without extension" << endl << endl;


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
  if (ptLS->BiCGSTAB(100,Jacobi)) cout << "convergence is true";
  ptLS->printxscreen();
  ptLS->check();
  // if (ptLS->GMRes_m(100, Jacobi, 100)) cout << "convergence is true";
#endif
  Driver * ptDriver=new Driver(ptInputfile,1,materialdata);
  ptDriver->SolveNewmarkMethod(ptUnverg);

  oClockTotal.ClockCount(Clock::end,"Total time");

/// Putzen

  if (ptInputfile) delete ptInputfile;
  if (ptUnverg) delete ptUnverg;
  if (ptDriver) delete ptDriver;
  if (materialdata) delete materialdata;
}
