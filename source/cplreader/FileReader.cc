// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <iomanip>

#include "boost/tokenizer.hpp"

// #include <muParser.h>

#include "General/exception.hh"
#include "Settings.hh"
#include "FileReader.hh"

namespace CoupledField
{

  FileReader::FileReader(const std::string& name,
                         const UInt dim,
                         const UInt numFiles,
                         const UInt startIndex) :
    numRegions_(0),
    dim_(dim),
    numSteps_(numFiles),
    startIndex_(startIndex),
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

      it = t.begin();
      end = t.end();

      for( ; it != end; it++)
      {
        requiredResults_[SolutionTypeEnum.Parse(*it)] = true;
      }
    }

    //TODO> This needs refactoring!
    if(requiredResults_[ACOU_RHS_LOAD]  || settings.GetInt("calcsrc") || requiredResults_[ACOU_LAMB_RHS] )
      requiredResults_[FLUIDMECH_VELOCITY] = true;
    if((requiredResults_[ACOUMIXED_MASS_LOAD] || requiredResults_[NO_SOLUTION_TYPE]) && settings.GetInt("calcsrc"))
      requiredResults_[FLUIDMECH_PRESSURE] = true;

    if(requiredResults_[ACOU_RHS_LOAD]){
      std::string lhSrcQuant = settings.GetString("quantityForAcouRhsLoad");
      //parse to solution type
      SolutionType lhType = SolutionTypeEnum.Parse(lhSrcQuant);
      switch(lhType){
        case FLUIDMECH_VELOCITY:
          requiredResults_[FLUIDMECH_VELOCITY] = true;
          break;
        case FLUIDMECH_PRESSURE:
          requiredResults_[FLUIDMECH_PRESSURE] = true;
          break;
        case FLUIDMECH_DIV_LH_T:
          requiredResults_[ACOU_DIV_LH_TENSOR_NODAL] = true;
          break;
        case FLUIDMECH_PRESSURE_DERIV_2:
          requiredResults_[FLUIDMECH_PRESSURE_DERIV_2] = true;
          break;
        default:

          break;
      }
    }
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
