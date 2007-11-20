#include <iostream>
#include <fstream>
#include <string>
#include <stdarg.h>
#include <vector>
#include <sstream>
#include <stdlib.h>

#include <cplreaderdefs.hh>

#include "params.hh"
#include "settings.hh"
#include "filereader.hh"

#ifdef CPLREADER_FASTEST
#include "FASTEST/filereader_FASTEST.hh"
#endif

#ifdef CPLREADER_CFX
#include "CFX/filereader_CFX.hh"
#endif

#ifdef CPLREADER_OPENFOAM
#include "OPENFOAM/filereader_OPENFOAM.hh"
#endif

// #include "Stanford/filereader_Stanford.hh"
#include "CFX/cfx_fortran_defs.h"
#include "MpCCIexch.hh"

using namespace CoupledField;

int main(int argc, char *argv[])
{
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
#ifdef CPLREADER_FASTEST
      fileReader = new FileReader_FASTEST(settings.GetString("name"),
                                          settings.GetInt("dim"),
                                          settings.GetInt("numSteps"));
#else
      EXCEPTION("Reading of FASTEST files not supported!");
#endif
    }

 #if 0
    if(type == "Stanford")
    {
      fileReader = new FileReader_Stanford(params->GetName(),
                                           params->GetDim(),
                                           params->GetNumSteps());
    }
 #endif

    if(type == "CFX")
    {
#ifdef CPLREADER_CFX
      fileReader = new FileReader_CFX(settings.GetString("name"),
                                      settings.GetInt("dim"),
                                      settings.GetInt("numSteps"));
#else
      EXCEPTION("Reading of CFX files not supported!");
#endif
    }

    if(type == "OPENFOAM")
    {
#ifdef CPLREADER_OPENFOAM
      fileReader = new FileReader_OPENFOAM(settings.GetString("name"),
                                           settings.GetInt("dim"),
                                           settings.GetInt("numSteps"));
#else
      EXCEPTION("Reading of OPENFOAM files not supported!");
#endif
    }

    if(!fileReader)
    {
      std::cerr << "ERROR: Could not initialize " << type
                << " filereader." << std::endl;
      return 0;
    }
    
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
      
