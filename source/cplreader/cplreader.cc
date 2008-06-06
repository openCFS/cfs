#include <iostream>
#include <fstream>
#include <string>
#include <stdarg.h>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>

#include <def_use_mpcci.hh>

#if (MpCCI_RELEASE == 305)
#include <mpcci.h>
#endif

#include <cplreaderdefs.hh>

#include "params.hh"
#include "settings.hh"
#include "filereader.hh"

#ifdef CPLREADER_FASTEST
#include "FASTEST/filereader_FASTEST.hh"
#endif

#ifdef CPLREADER_CFX
#include "CFX/filereader_CFX.hh"
#include "CFXexport/FileReader_CFXexport.hh"
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

    /*std::cout << "Name: " << settings.GetString("name")
              << " Type: " << settings.GetString("type")
              << " Dim: " << settings.GetInt("dim")
              << " CalcSrc: " << settings.GetInt("calcSrc")
              << " Verbose: " << settings.GetInt("verbose")
              << " numfiles: " << settings.GetInt("numSteps")
      //              << " LD_LIBRARY_PATH: " << getenv("LD_LIBRARY_PATH")
              << " PWD: " << getenv("PWD")
              <<std::endl;*/
    
    if(type == "FASTEST")
    {
#ifdef CPLREADER_FASTEST
      fileReader = new FileReader_FASTEST(settings.GetString("name"),
                                          settings.GetInt("dim"),
                                          settings.GetInt("numSteps"),
                                          settings.GetInt("firstStep"));
#else
      EXCEPTION("Reading of FASTEST files not supported!");
#endif
    }

 #if 0
    if(type == "Stanford")
    {
      fileReader = new FileReader_Stanford(settings.GetString("name"),
                                           settings.GetInt("dim"),
                                           settings.GetInt("numSteps"));
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

    if(type == "CFX_EXPORT")
    {
#ifdef CPLREADER_CFX
      fileReader = new FileReader_CFXexport(settings.GetString("name"),
                                            settings.GetInt("dim"),
                                            settings.GetInt("numSteps"),
                                            settings.GetInt("firstStep"));
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
    

    /* check if file even exists */
    struct stat st;
    const char* const tmp_fileName = settings.GetString("name").c_str();
    if(lstat(tmp_fileName, &st) == -1 )
    {
      EXCEPTION("the file doesn't exist: " << tmp_fileName << std::endl);
    }


    MpCCIExchangeCPLR mpCCIexch(fileReader);
    mpCCIexch.Init(argc, argv);
    if(!settings.GetInt("justinit")) 
    {
      mpCCIexch.PutExchangeGrid2MpCCI();

      if(!settings.GetInt("justmesh")) 
      {
        mpCCIexch.Couple();
      }
    }
    mpCCIexch.Finish();

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
      
