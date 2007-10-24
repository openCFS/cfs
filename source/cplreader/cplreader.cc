#include <iostream>
#include <fstream>
#include <string>
#include <stdarg.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#ifdef MpCCI
#include <mpcci.h>
#endif


#include "params.hh"
#include "settings.hh"
#include "filereader.hh"
#include "FASTEST/filereader_FASTEST.hh"
// #include "Stanford/filereader_Stanford.hh"
#include "CFX/filereader_CFX.hh"
#include "CFX/cfx_fortran_defs.h"
#include "MpCCIexch.hh"

#define REALTYPE CCI_DOUBLE
  typedef double Realtype;

 
using namespace CoupledField;

int main(int argc, char *argv[])
{
#ifdef MpCCI
  CCI_Init_with_id_string( &argc, &argv, "simulationcode1" );
#endif
  int ret = 0;
  FileReader* fileReader = NULL;

  try 
  {
      
    Settings& settings = Settings::Instance();


    ParamsInit(argc, argv);
    std::string type = settings.GetString("type");

    std::cout << "Name: " << settings.GetString("name")
              << " Type: " << settings.GetString("type")
              << " Dim: " << settings.GetInt("dim")
              << " CalcSrc: " << settings.GetInt("calcSrc")
              << " Verbose: " << settings.GetInt("verbose")
              << " numfiles: " << settings.GetInt("numSteps")
      //              << " LD_LIBRARY_PATH: " << getenv("LD_LIBRARY_PATH")
              << " PWD: " << getenv("PWD")
              <<std::endl;
    
    if(type == "FASTEST")
    {
      fileReader = new FileReader_FASTEST(settings.GetString("name"),
                                          settings.GetInt("dim"),
                                          settings.GetInt("numSteps"));
    }

 #if 0
    if(type == "Stanford")
    {
      fileReader = new FileReader_Stanford(params->GetName(),
                                           params->GetDim(),
                                           params->GetNumSteps());
    }
 #endif

    if(settings.GetString("type") == "CFX")
    {
      fileReader = new FileReader_CFX(settings.GetString("name"),
                                      settings.GetInt("dim"),
                                      settings.GetInt("numSteps"));
    }


    if(!fileReader)
    {
      std::cerr << "ERROR: Could not initialize " << settings.GetString("type")
                << " filereader." << std::endl;
      return 0;
    }
    
    //    CCI_Init(&argc, &argv);

    MpCCIexch mpCCIexch(argc, argv, fileReader);
    mpCCIexch.PutExchangeGrid2MpCCI();
    mpCCIexch.Couple();

  } catch (std::exception& ex)
  {
    std::cerr << "CAUGHT EXCEPTION:" << std::endl 
              << std::endl
              << ex.what()
              << std::endl;
    ret = 1;
  }
    
  delete fileReader;
    
  return ret;
}
      
