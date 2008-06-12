#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>

// #include <muParser.h>

#include "params.hh"
#include "filereader.hh"
#include "settings.hh"

namespace CoupledField
{

  FileReader::FileReader(const std::string& name,
                         const UInt dim,
                         const UInt numFiles) :
    numRegions_(0), 
    dim_(dim),
    numSteps_(numFiles),
    maxNumElemNodes_(0)
  {
    name_ = name;
    baseName_ = "./";
    baseName_+= name_;
    baseName_+= "/";
    baseName_+= name_;
    baseName_+= "_";
    
    preferedOutputPath_ = "cplreader_hdf5_";
    preferedOutputPath_ += name_; 
  }

  FileReader::~FileReader()
  {
  }

  double FileReader::GetTimeStep(UInt stepNumber)
  {
    Settings& settings = Settings::Instance();
    Double timestep = settings.GetDouble("timeStep");
    
    double ts = timestep * stepNumber;

    return ts;
  }

  std::string FileReader::GetRegionName(const UInt regionIdx)
  {
    std::ostringstream sstr;
    
    sstr << "partition" << (regionIdx+1);
    return sstr.str();
  }
      
} // end of namespace
