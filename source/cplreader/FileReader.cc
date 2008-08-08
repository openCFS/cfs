#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>

#include <boost/tokenizer.hpp>

// #include <muParser.h>

#include "General/exception.hh"
#include "Settings.hh"
#include "FileReader.hh"

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
    Settings& settings = Settings::Instance();

    name_ = name;
    baseName_ = settings.GetString("basedir");
    baseName_+= "/";
    baseName_+= name_;
    baseName_+= "/";
    baseName_+= name_;
    baseName_+= "_";

    preferredOutputPath_ = "cplreader_hdf5_";
    preferredOutputPath_ += name_;

    // Initialize vector with required results
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");
    Tok t(settings.GetString("outputfields"), sep);

    if(*t.begin() == "all")
    {
      requiredResults_[NO_SOLUTION_TYPE] = true;
    }
    else
    {
      Tok::const_iterator it, end;
      SolutionType st;

      it = t.begin();
      end = t.end();

      for( ; it != end; it++)
      {
        String2Enum(*it, st);
        requiredResults_[st] = true;
      }
    }

    if(requiredResults_[ACOU_RHS_LOAD] || settings.GetInt("calcsrc"))
      requiredResults_[FLUIDMECH_VELOCITY] = true;
  }

  FileReader::~FileReader()
  {
  }

  double FileReader::GetTimeStep(UInt stepNumber)
  {
    Settings& settings = Settings::Instance();
    Double timestep = settings.GetDouble("timestep");

    double ts = timestep * stepNumber;

    return ts;
  }

  std::string FileReader::GetRegionName(const UInt regionIdx)
  {
    std::ostringstream sstr;

    sstr << "partition" << (regionIdx+1);
    return sstr.str();
  }

  void FileReader::GetNodeGroups(std::map<std::string,
                                 std::vector<UInt> >& nodeGroups)
  {
    Settings& settings = Settings::Instance();
    if(settings.GetInt("verbose"))
      std::cerr << "GetNodeGroups() not implemented!" << std::endl;
  }

  void FileReader::GetElemGroups(std::map<std::string,
                                 std::vector<UInt> >& elemGroups)
  {
    Settings& settings = Settings::Instance();
    if(settings.GetInt("verbose"))
      std::cerr << "GetElemGroups() not implemented!" << std::endl;

  }

} // end of namespace
