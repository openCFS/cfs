// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <algorithm>
#include <cstring>
#include <iostream>
#include <set>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <typeinfo>

#include "boost/regex.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;
#include "boost/algorithm/string/predicate.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/tokenizer.hpp"
namespace algo=boost::algorithm;

#include "General/exception.hh"
#include "cplreader/Settings.hh"
#include "FileReader_CFX.hh"
#include "cfx_fortran_defs.h"

#define CHECK_CFX_IO(NERR)                      \
  if(NERR)                                      \
  {                                             \
    std::string errStr;                         \
    IOErrorToString(NERR, errStr);              \
    EXCEPTION(errStr);                          \
  }


namespace CoupledField
{
  std::vector<int> FileReader_CFX::intvec_;
  std::vector<double> FileReader_CFX::doublevec_;
  std::vector<float> FileReader_CFX::floatvec_;
  std::vector<char> FileReader_CFX::charvec_;

  int ind = 0;

  FileReader_CFX::FileReader_CFX(const std::string& name,
                                 const UInt dim,
                                 const UInt numFiles) :
    FileReader(name, dim, numFiles, 1),
    multiLine_(false),
    timeStep_(0.0),
    determineFloatDS_(true),
    numZones_(0),
    numVolumes_(0),
    numBCPs_(0),
    numNodes_(0),
    numElems_(0),
    numVolElems_(0),
    numSurfElems_(0),
    whatfile_(0)
  {
    Settings& settings = Settings::Instance();

    name_ = name;
    baseName_ = settings.GetString("basedir");
    baseName_+= "/";
    baseName_+= name_;
    baseName_+= "/";
  }

  FileReader_CFX::~FileReader_CFX()
  {
  }

  template<typename T>
  T FileReader_CFX::ReadScalar(std::string what,
                               std::string where,
                               Integer when,
                               bool throwEx)
  {
    T data;
    int dattyp = 0;
    int length = 1;
    int nsize  = 1;
    int iopt   = __not_stop_if_failed__;
    int ioptpar = 0;
    int nerr   = 0;

    snprintf(what_, sizeof(what_), "%s", what.c_str());
    snprintf(where_, sizeof(where_), "%s", where.c_str());
    
    redsht_(&dattyp, &length, &nerr,
            what_, where_, &when,
            &nsize, &iopt, &ioptpar,
            (float*) &data, (int*) &data, NULL, NULL, (double*) &data, NULL,
            strlen(what_), strlen(where_), 0);
    
    if ( nerr == __io_ok__ ) {
      return data;
    } else {
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return T();
    }
  }
  
  std::string FileReader_CFX::ReadString(std::string what,
                                         std::string where,
                                         Integer when,
                                         bool throwEx)
  {
    int dattyp = 0;
    int length = 1;
    int nsize  = sizeof(carr_);
    int iopt   = __not_stop_if_failed__;
    int ioptpar = 0;
    int nerr   = 0;

    snprintf(what_, sizeof(what_), "%s", what.c_str());
    snprintf(where_, sizeof(where_), "%s", where.c_str());
    
    redsht_(&dattyp, &length, &nerr,
            what_, where_, &when,
            &nsize, &iopt, &ioptpar,
            NULL, NULL, carr_, NULL, NULL, sarr_,
            strlen(what_), strlen(where_), 0);
    
    if ( nerr == __io_ok__ ) {
      switch ( dattyp ) {
      case __char_data_type__:
        return boost::algorithm::trim_right_copy(
            std::string(carr_, sizeof(carr_)));
      case __string_data_type__:
        return boost::algorithm::trim_right_copy(
            std::string(sarr_, sizeof(sarr_)));
      default:
        if (throwEx) EXCEPTION("Dataset does not contain a string");
        return std::string();
      }
    } else {
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return std::string();
    }
  }

  template<>
  bool FileReader_CFX::ReadShortVector(std::vector<std::string> &out,
                                       Integer size,
                                       std::string what,
                                       std::string where,
                                       Integer when,
                                       bool throwEx)
  {
    int dattyp = 0;
    int length  = 1;
    int iopt   = __not_stop_if_failed__;
    int ioptar = 0;
    int nerr   = 0;

    if ( size < 1 ) size = 1;
    
    snprintf(what_, sizeof(what_), "%s", what.c_str());
    snprintf(where_, sizeof(where_), "%s", where.c_str());
    
    charvec_.resize(IO_BUFSIZE*size);

    redsht_(&dattyp, &length, &nerr,
            what_, where_, &when,
            &size, &iopt, &ioptar,
            NULL, NULL, NULL, NULL, NULL, charvec_.data(),
            strlen(what_), strlen(where_), 0);
    
    if ( nerr == __io_ok__ ) {
      if ( dattyp != __string_data_type__ ) {
        if (throwEx) EXCEPTION("Dataset does not contain a string array");
        return false;
      }
      
      out.clear();
      out.reserve(size);

      for ( UInt i=0; i<(UInt)size; ++i )
      {
        out.push_back(
            boost::algorithm::trim_right_copy(
                std::string(&charvec_[i*IO_BUFSIZE], IO_BUFSIZE)));
      }
      
      return true;
    } else {
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return false;
    }
  }
  
  template<typename T>
  bool FileReader_CFX::ReadShortVector(std::vector<T> &out,
                                       Integer size,
                                       std::string what,
                                       std::string where,
                                       Integer when,
                                       bool throwEx)
  {
    int dattyp = 0;
    int length = 1;
    int iopt   = __not_stop_if_failed__;
    int ioptpar = 0;
    int nerr   = 0;
    T *datarray;
    
    if (size < 1) size = 1;
    out.resize(size*length);
    datarray = out.data();
    
    snprintf(what_, sizeof(what_), "%s", what.c_str());
    snprintf(where_, sizeof(where_), "%s", where.c_str());
    

    redsht_(&dattyp, &length, &nerr,
        what_, where_, &when,
        &size, &iopt, &ioptpar,
        (float*)datarray, (int*)datarray, NULL, NULL, (double*)datarray, NULL,
        strlen(what_), strlen(where_), 0);

    if ( nerr == __io_ok__ ) {
      out.resize(size*length);
      return true;
    } else { 
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return false;
    }

  }

  template<typename T>
  bool FileReader_CFX::ReadLongVector(std::vector<T> &out,
                                      Integer size,
                                      std::string what,
                                      std::string where,
                                      Integer when,
                                      bool throwEx)
  {
    int dattyp;
    int length = 1;
    int iopt   = __not_stop_if_failed__;
    int nerr   = 0;
    T *datarray;
    
    if (size < 1) size = 1;
    out.resize(size*length);
    datarray = out.data();
    
    const std::type_info &tiOut = typeid(T);
    if ( tiOut == typeid(float) ) {
      dattyp = __real_data_type__;
    } else if ( tiOut == typeid(int) || tiOut == typeid(unsigned int) ) {
      dattyp = __int_data_type__;
    } else if ( tiOut == typeid(double) ) {
      dattyp = __double_data_type__;
    } else {
      if (throwEx) EXCEPTION("Data type not supported: " << tiOut.name());
      return false;
    }

    snprintf(what_, sizeof(what_), "%s", what.c_str());
    snprintf(where_, sizeof(where_), "%s", where.c_str());

    readlong_(&dattyp, &nerr,
        what_, where_, &when,
        &size, &iopt,
        (float*)datarray, (int*)datarray, NULL, NULL, (double*)datarray, NULL,
        strlen(what_), strlen(where_), 0);

    if ( nerr == __io_ok__ ) {
      return true;
    } else { 
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return false;
    }

  }
  
  bool FileReader_CFX::OpenCFXFile(const std::string &filename, bool throwEx) {
    int nerr;
    int whatfile = __io_open_primaryfile__;
    
    snprintf(fn_, sizeof(fn_), "%s", filename.c_str());

    openfile_(&nerr, fn_, &whatfile, strlen(fn_));
    
    if ( nerr == __io_ok__ ) {
      return true;
    } else {
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return false;
    }
  }
  
  bool FileReader_CFX::CloseCFXFile(bool throwEx) {
    int nerr;
    int whatfile = __io_close_primaryfile__;
    
    closefile_(&nerr, &whatfile);
    
    if ( nerr == __io_ok__ ) {
      return true;
    } else {
      if ( throwEx ) CHECK_CFX_IO(nerr);
      return false;
    }
  }
  
  void FileReader_CFX::Init()
  {
    Settings& settings = Settings::Instance();
    
    bool beVerbose = settings.GetInt("verbose");
    int nerr = 0;
    
    std::stringstream sstr;
    sstr << baseName_ << name_ << ".res";
    std::string resFileName = sstr.str();
    std::cout << "resFileName: " << resFileName << std::endl;

    if(fs::exists(resFileName))
    {

      //-----------------------------------------------------------------------
      //     Open RESULTS file
      //-----------------------------------------------------------------------

      if ( beVerbose ) {
        std::cout << "Trying to open CFX results file " << resFileName
                  << "." << std::endl;
      }

      OpenCFXFile(resFileName);

      //-----------------------------------------------------------------------
      //     Read the latest time step number
      //-----------------------------------------------------------------------
      
      sprintf(what_, "G/TRANSIENT");
      sprintf(where_, "ZN1");
      
      int when   = -1;
      int dattyp = 0;
      int length = 1;
      int nsize  = 3;
      int iopt   = __not_stop_if_failed__;
      int ioptpar =   0;
      intvec_.resize(length*nsize);
      floatvec_.resize(length*nsize);
      
      redsht_(&dattyp, &length, &nerr, what_, where_, &when, &nsize,
              &iopt, &ioptpar,
              floatvec_.data(), intvec_.data(), NULL, NULL, NULL, NULL,
              strlen(what_), strlen(where_), 0);
      CHECK_CFX_IO(nerr);
      
      int its = intvec_[2];

      //-----------------------------------------------------------------------
      //     Read the trninfo data set
      //-----------------------------------------------------------------------
      int ntrn;

      sprintf(what_, "G/TRN_INFO");
      sprintf(where_, "EVERY");
      
      when   = its;
      dattyp =  -3;
      length =   1;
      nsize  = its;
      iopt   = __not_stop_if_failed__;
      ioptpar =   0;
      floatvec_.resize(length*nsize);
      intvec_.resize(length*nsize);
      charvec_.resize(length*nsize*80);

      redsht_(&dattyp, &length, &nerr,
              what_, where_, &when,
              &nsize, &iopt, &ioptpar,
              &floatvec_[0], &intvec_[0], carr_, NULL, darr_, &charvec_[0],
              strlen(what_), strlen(where_), 0);
      CHECK_CFX_IO(nerr);

      ntrn = nsize;

      if ( settings.GetDouble("timestep") <= 0 )
      {
        settings.SetDouble("timestep",
            (floatvec_[1]-floatvec_[0])/(float)(intvec_[1]-intvec_[0]));
      }

      numSteps_ = 0;
      trnFilenames_.clear();
      timeStepNumbers_.clear();
      trnFilenames_.reserve(ntrn);
      timeStepNumbers_.reserve(ntrn);

      for(int i = 0; i< ntrn; i++)
      {
        int its  = intvec_[i];

        char* trnnam = &charvec_[i*80];
        int n;
        for(n = 0; n< 80; n++)
        {
          if(isblank((int)trnnam[n]))
          {
            trnnam[n] = 0;
            break;
          }
        }

        std::string fileName;
        sstr.clear();
        sstr.str("");
        sstr << baseName_ << name_ << "/" << trnnam;
        fileName = sstr.str();

        inFile_.clear();
        inFile_.open(fileName.c_str());
        if (inFile_)
        {
          inFile_.close();

          numSteps_++;
          trnFilenames_.push_back(fileName);
          timeStepNumbers_.push_back(its);
        }

      }

      //-----------------------------------------------------------------------
      //     read the G/COMMANDS string from the res file
      //     to get the name of the def file and the time step
      //-----------------------------------------------------------------------

      std::vector<std::string> commands;
      ReadShortVector(commands, 10000, "G/COMMANDS");

      CloseCFXFile();

      //-----------------------------------------------------------------------
      //     Parse the command string from the CFX results file and
      //     get infos about definition file, time unit and timestep
      //-----------------------------------------------------------------------

      GetInfosFromCommand(commands);
    }
    else
    {
      fs::path trnDir( baseName_ +  name_);
      fs::directory_iterator end_iter;
      std::set<UInt> stepNumSet;
      std::set<UInt>::const_iterator it, end;
      UInt stepNum;
      std::string fn;

      for ( fs::directory_iterator dir_itr( trnDir );
            dir_itr != end_iter;
            ++dir_itr )
      {
        if ( !fs::is_directory( *dir_itr ) )
        {
          fn = dir_itr->path().filename().string();

          if(algo::ends_with(fn, ".trn"))
          {
            sstr.clear(); sstr.str("");
            sstr << fn;
            sstr >> stepNum;
            stepNumSet.insert(stepNum);
          }
        }
      }

      it = stepNumSet.begin();
      end = stepNumSet.end();
      
      trnFilenames_.reserve(stepNumSet.size());
      timeStepNumbers_.reserve(stepNumSet.size());

      for( ; it != end; it++ )
      {
        stepNum = *it;
        sstr.clear(); sstr.str("");
        sstr << baseName_ << name_ << "/" << stepNum << ".trn";

        fn = sstr.str();

        numSteps_++;
        trnFilenames_.push_back(fn);
        timeStepNumbers_.push_back(stepNum);
      }
    }

    if(settings.GetDouble("timestep") <= 0)
      EXCEPTION("Time step could not be determined. Please specify it using --timestep X.");

    CheckTransientFiles();

    if (settings.GetInt("numsteps"))
    {
      UInt tmp = (UInt) settings.GetInt("numsteps");
      /* only take argument if tmp does not exceed the maximal number of timesteps possible */
      if (tmp < numSteps_)
      {
        numSteps_ = tmp;
      }
    }

    //-----------------------------------------------------------------------
    //     Open DEFINITION file
    //-----------------------------------------------------------------------

    std::vector< std::string > defFileNames;
    if(settings.GetString("deffile") != "")
      defFileNames.push_back(baseName_ + settings.GetString("deffile").c_str());
    if(defFile_ != "")
      defFileNames.push_back(baseName_ + defFile_);
    defFileNames.push_back(baseName_ + name_ + ".def");
    defFile_ = "";

    for(UInt i=0; i<defFileNames.size(); i++)
    {
      if(beVerbose)
      {
        std::cerr << "Trying to open deffile: " << defFileNames[i] << " ";
      }

      inFile_.clear();
      inFile_.open(defFileNames[i].c_str());
      if (inFile_)
      {
        defFile_ = defFileNames[i];
        inFile_.close();
        if(beVerbose)
        {
          std::cerr << "-> OK!" << std::endl;
        }

        break;
      }

      if(beVerbose)
      {
        std::cerr << "-> failed!" << std::endl;
      }
    }

    if(defFile_ == "")
    {
      EXCEPTION("Cannot find definition file.");
    }

    OpenCFXFile(defFile_);
    
    //-----------------------------------------------------------------------
    // number of zones
    //-----------------------------------------------------------------------
    
    numZones_ =  ReadScalar<UInt>("G/NZN");
    if ( beVerbose ) {
      printf("Number of domains: %d\n", numZones_);
    }
    
    nodeOffsetPerZone_.resize(numZones_);
    elemOffsetPerZone_.resize(numZones_);
    
    std::vector< std::string > zoneNames;
    zoneNames.resize(numZones_);
    
    //-----------------------------------------------------------------------
    // number of volume regions
    //-----------------------------------------------------------------------
    
    numVolumes_ = ReadScalar<UInt>("G/NVL");
    if ( beVerbose ) {
      printf("Number of volume regions: %d\n", numVolumes_);
    }

    volumes_.resize(numVolumes_);
    
    //-----------------------------------------------------------------------
    // number of surface regions
    //-----------------------------------------------------------------------
    
    numBCPs_ = ReadScalar<UInt>("G/NBCP");
    if ( beVerbose ) {
      printf("Number of surface regions: %d\n", numBCPs_);
    }

    faceSetsPerBCP_.resize(numBCPs_);
    numRegions_ = numVolumes_ + numBCPs_;
    numNodesPerRegion_.resize(numRegions_, 0);
    numElemsPerRegion_.resize(numRegions_, 0);
    regionNames_.resize(numRegions_);
    
    //-----------------------------------------------------------------------
    // Names of domains and boundary condition patches
    //-----------------------------------------------------------------------
    Integer nzif = ReadScalar<Integer>("G/NZIF", "EVERY", 0, false),
            nvp = ReadScalar<Integer>("G/NVP", "EVERY", 0, false);
    UInt tmp;
    std::vector< std::string > nameMap;
    std::istringstream iss;
    ReadShortVector(nameMap, 4*(numZones_+numBCPs_+nzif+nvp), "G/NAMEMAP");
    
    for ( UInt i=0, numEntries=nameMap.size()/4; i<numEntries; ++i ) {
      if ( nameMap[4*i+3] == "ZN" ) {
        nameMap[4*i+1].replace(0, 2, "");
        iss.clear();
        iss.str(nameMap[4*i+1]);
        iss >> tmp;
        if ( tmp > 0 && tmp <= numZones_ ) {
          replace(nameMap[4*i].begin(), nameMap[4*i].end(), ' ', '_');
          zoneNames[tmp-1] = nameMap[4*i];
        }
      } else if ( nameMap[4*i+3] == "BCP" ) {
        nameMap[4*i+1].replace(0, 3, "");
        iss.clear();
        iss.str(nameMap[4*i+1]);
        iss >> tmp;
        if ( tmp > 0 && tmp <= numBCPs_ ) {
          replace(nameMap[4*i].begin(), nameMap[4*i].end(), ' ', '_');
          regionNames_[tmp-1+numVolumes_] = nameMap[4*i];
        }
      }
    }

    //------------------------------------------------------------------------
    // Loop over all zones
    //------------------------------------------------------------------------
    Integer nvx, nel, nvl, nes;
    UInt numElemNodes, currVol = 0;
    ElemSet newES;
    std::ostringstream zoneStr;
    std::vector<Integer> kesvl, ipesvl, iltpes, neles;
    std::vector<std::string> volNames;

    for ( UInt iZone=1; iZone<=numZones_; ++iZone ) {
      zoneStr.str("");
      zoneStr << "ZN" << iZone;

      //----------------------------------------------------------------------
      // number of vertices: NVX
      //----------------------------------------------------------------------
      nvx = ReadScalar<Integer>("G/NVX", zoneStr.str());
      
      nodeOffsetPerZone_[iZone-1] = numNodes_;
      numNodes_ += (UInt) nvx;

      //----------------------------------------------------------------------
      // number of elements: NEL
      //----------------------------------------------------------------------
      nel = ReadScalar<Integer>("G/NEL", zoneStr.str());
      
      elemOffsetPerZone_[iZone-1] = numElems_;
      numElems_ += (UInt) nel;
      
      //----------------------------------------------------------------------
      // number of volumes: NVL
      //----------------------------------------------------------------------
      nvl = ReadScalar<Integer>("G/NVL", zoneStr.str());
      
      //----------------------------------------------------------------------
      // volume labels: CLBVL
      //----------------------------------------------------------------------
      ReadShortVector(volNames, nvl, "G/CLBVL", zoneStr.str());
      
      //----------------------------------------------------------------------
      // number of element sets: NES
      //----------------------------------------------------------------------
      nes = ReadScalar<Integer>("G/NES", zoneStr.str());
      
      //----------------------------------------------------------------------
      // number of elements per element sets: NELES
      //----------------------------------------------------------------------
      ReadShortVector(neles, nes, "G/NELES", zoneStr.str());
      
      //----------------------------------------------------------------------
      // types of element sets: ILTPES
      //----------------------------------------------------------------------
      //---   elem type = 4: tet  , 4 nodes
      //---             = 5: wedge, 6 nodes
      //---             = 6: hex,   8 nodes
      //---             = 7: pyr,   5 nodes
      ReadShortVector(iltpes, nes, "G/ILTPES", zoneStr.str());
      
      //----------------------------------------------------------------------
      // element sets per volume: KESVL, IPESVL
      //----------------------------------------------------------------------
      ReadShortVector(kesvl, nes, "G/KESVL", zoneStr.str());
      ReadShortVector(ipesvl, nvl+1, "G/IPESVL", zoneStr.str());

      for ( Integer iVol=0; iVol<nvl; ++iVol, ++currVol ) {
        volumes_[currVol].zone = iZone;
        regionNames_[currVol] = zoneNames[iZone-1] + "." + volNames[iVol];
        
        for ( Integer ies=ipesvl[iVol]; ies<ipesvl[iVol+1]; ++ies ) {
          newES.num = (UInt) kesvl[ies-1];
          newES.numElems = (UInt) neles[ies-1];
          switch (iltpes[ies-1])
          {
            case 4:
              newES.elemType = Elem::TET4;
              break;
            case 5:
              newES.elemType = Elem::WEDGE6;
              break;
            case 6:
              newES.elemType = Elem::HEXA8;
              break;
            case 7:
              newES.elemType = Elem::PYRA5;
              break;
            default:
              newES.elemType = Elem::UNDEF;
              break;
          }

          volumes_[currVol].elemSets.push_back(newES);
          
          numElemNodes = Elem::GetNumElemNodes(newES.elemType);
          if ( numElemNodes > maxNumElemNodes_ )
            maxNumElemNodes_ = numElemNodes;

          numElemsPerRegion_[currVol] += newES.numElems;
        }
        
        if (beVerbose)
        {
          std::cout << "- Volume region  '" << regionNames_[currVol] << "':\t"
              << numElemsPerRegion_[currVol] << " Elements\n";
        }
      } // loop over iVol
    } // loop over iZone

    numVolElems_ = numElems_;

    //-----------------------------------------------------------------------
    // Surfaces and face sets
    //-----------------------------------------------------------------------
    Integer nfs, nsf, sfIdx;
    std::vector<Integer> ksffs, ipsffs, nfcfs;
    std::vector<FaceSet> surfaces;
    
    nsf = ReadScalar<Integer>("G/NSF");   // total number of surfaces
    surfaces.resize(nsf);

    for ( UInt iZone=1; iZone<=numZones_; ++iZone ) {
      zoneStr.str("");
      zoneStr << "ZN" << iZone;

      nfs = ReadScalar<Integer>("G/NFS", zoneStr.str());
      nsf = ReadScalar<Integer>("G/NSF", zoneStr.str());
      if ( nfs != nsf ) {
        WARN("The number of face sets (NFS=" << nfs
             << ") is not equal to the number of surfaces (NSF=" << nsf
             << ") in domain '" << zoneStr.str()
             << "'. This case has never been tested!");
      }

      
      ReadShortVector(ksffs, nsf, "G/KSFFS", zoneStr.str());
      ReadShortVector(ipsffs, nfs+1, "G/IPSFFS", zoneStr.str());
      ReadShortVector(nfcfs, nfs, "G/NFCFS", zoneStr.str());
      
      for ( Integer ifs=0; ifs<nfs; ++ifs ) {
        for ( Integer iSurf=ipsffs[ifs]; iSurf<ipsffs[ifs+1]; ++iSurf ) {
          sfIdx = ksffs[iSurf-1] - 1;
          surfaces[sfIdx].zone = iZone;
          surfaces[sfIdx].num = ifs+1;
          surfaces[sfIdx].numElems = nfcfs[ifs];
        }
      }
    }
    
    numSurfElems_ = numElems_ - numVolElems_;
    
    //-----------------------------------------------------------------------
    // Surfaces per boundary condition patch
    //-----------------------------------------------------------------------
    UInt numFaces;
    std::vector<Integer> ksfbcp, ipsfbcp;
    ReadShortVector(ksfbcp, surfaces.size(), "G/KSFBCP");
    ReadShortVector(ipsfbcp, numBCPs_+1, "G/IPSFBCP");
    
    for ( UInt ibcp=0; ibcp<numBCPs_; ++ibcp ) {
      numFaces = 0;
      faceSetsPerBCP_[ibcp].reserve(ipsfbcp[ibcp+1]-ipsfbcp[ibcp]);
      
      for ( Integer isf=ipsfbcp[ibcp]; isf<ipsfbcp[ibcp+1]; ++isf ) {
        if ( ksfbcp[isf-1] > 0) { // ksfbcp[isf-1]==0 means surface is unused
          faceSetsPerBCP_[ibcp].push_back(surfaces[ksfbcp[isf-1]-1]);
          numFaces += surfaces[ksfbcp[isf-1]-1].numElems;
        }
      }
      
      numElemsPerRegion_[numVolumes_+ibcp] = numFaces;
      numElems_ += numFaces;
      
      if (beVerbose)
      {
        std::cout << "- Surface region '" << regionNames_[numVolumes_+ibcp]
                  << "':\t" << numFaces << " Elements\n";
      }
    }

    if ( beVerbose ) {
      printf("Total number of nodes: %d\n", numNodes_);
      printf("Total number of elements: %d\n", numElems_);
    }

    //
    //---- read CFX Release No. for UserData
    //
    userDataCFXRelease_ = ReadString("G/CRELNO", "EVERY", 0, false);
    
    //
    //---- close file
    //
    CloseCFXFile();
  }

  void FileReader_CFX::ReadNodalCoords(std::vector<Double> & nodeCoord)
  {
    Settings& settings = Settings::Instance();
    
    bool beVerbose = settings.GetInt("verbose");
    Integer size;
    std::vector<Double> buffer;
    std::ostringstream zoneStr;
    
    nodeCoord.resize(numNodes_*3);

    OpenCFXFile(defFile_);

    for ( UInt iZone=1; iZone<=numZones_; ++iZone ) {
      zoneStr.str("");
      zoneStr << "ZN" << iZone << "/VX";
      
      size = (iZone==numZones_) ? numNodes_ - nodeOffsetPerZone_[iZone-1]
             : nodeOffsetPerZone_[iZone] - nodeOffsetPerZone_[iZone-1];
      
      if (beVerbose) {
        std::cout << "Reading " << size << " nodes in domain " << iZone
                  << "... ";
        std::cout.flush();
      }

      if (!ReadLongVector(buffer, 3*size, "G/CRDVX", zoneStr.str(), 0, false))
      {
        ReadShortVector(buffer, 3*size, "G/CRDVX", zoneStr.str());
      }
      std::copy(buffer.begin(), buffer.end(),
                nodeCoord.begin() + 3*nodeOffsetPerZone_[iZone-1]);
      
      if (beVerbose) {
        std::cout << "done.\n";
      }
    }
    
    CloseCFXFile();
  }

  void FileReader_CFX::ReadTopology( std::vector<UInt> & topologyData,
                                     std::vector<UInt> & elemTypes )
  {
    Settings& settings = Settings::Instance();
    
    bool beVerbose = settings.GetInt("verbose");
    UInt i, j;
    UInt nodeOffset, numElemNodes, locIdx, globIdx, volElemIdx;
    std::vector<UInt> elConnect(maxNumElemNodes_);
    std::vector<UInt> buffer1, buffer2; 
    std::ostringstream zoneStr;

    std::cout  << "Reading element connectivity" << std::endl;

    // open .def file
    OpenCFXFile(defFile_);

    // allocate memory
    topologyData.resize(numElems_ * maxNumElemNodes_, 0);
    elemTypes.resize(numElems_, Elem::UNDEF);
    regionElems_.resize(numVolumes_+numBCPs_);
    
    // first read volume regions
    for ( UInt iVol=0; iVol<numVolumes_; ++iVol ) {
      std::vector<ElemSet>::iterator esIt = volumes_[iVol].elemSets.begin(),
                                     esEnd = volumes_[iVol].elemSets.end();
      
      regionElems_[iVol].reserve(numElemsPerRegion_[iVol]);
      
      for ( ; esIt != esEnd; ++esIt ) {
        zoneStr.str("");
        zoneStr << "ZN" << volumes_[iVol].zone << "/ES" << esIt->num;
        
        ReadLongVector(buffer1, esIt->numElems, "G/KELPE", zoneStr.str(), 0);
        
        numElemNodes = Elem::GetNumElemNodes(esIt->elemType);
        ReadLongVector(buffer2, esIt->numElems*numElemNodes,
                       "G/KVXPE", zoneStr.str());
        
        nodeOffset = nodeOffsetPerZone_[volumes_[iVol].zone-1];
        locIdx = 0;
        
        for ( i=0; i<esIt->numElems; ++i, locIdx+=numElemNodes ) {
          globIdx = buffer1[i] + elemOffsetPerZone_[volumes_[iVol].zone-1];
          regionElems_[iVol].push_back( globIdx );
          elemTypes[--globIdx] = esIt->elemType;
          
          globIdx *= maxNumElemNodes_;
          if (esIt->elemType == Elem::HEXA8) {
            topologyData[ globIdx + 0 ] = nodeOffset + buffer2[locIdx + 4];
            topologyData[ globIdx + 1 ] = nodeOffset + buffer2[locIdx + 6];
            topologyData[ globIdx + 2 ] = nodeOffset + buffer2[locIdx + 7];
            topologyData[ globIdx + 3 ] = nodeOffset + buffer2[locIdx + 5];
            topologyData[ globIdx + 4 ] = nodeOffset + buffer2[locIdx + 0];
            topologyData[ globIdx + 5 ] = nodeOffset + buffer2[locIdx + 2];
            topologyData[ globIdx + 6 ] = nodeOffset + buffer2[locIdx + 3];
            topologyData[ globIdx + 7 ] = nodeOffset + buffer2[locIdx + 1];
          } else {
            for ( j=0; j<numElemNodes; ++j )  {
              topologyData[globIdx + j] = nodeOffset + buffer2[locIdx + j];
            }
          }
          
          volumes_[iVol].nodes.insert(buffer2.begin()+locIdx,
                                      buffer2.begin()+locIdx+numElemNodes);
        }
      }
      
      numNodesPerRegion_[iVol] = volumes_[iVol].nodes.size();
      
      if (beVerbose) {
        std::cout << "Read connectivity of volume region '"
                  << regionNames_[iVol] << "'.\n";
      }
    }
    
    // start numbering of surface elements just after last volume element
    UInt surfElemIdx = numVolElems_;
    
    for ( UInt ibcp=0; ibcp<numBCPs_; ++ibcp ) {
      std::vector<FaceSet>::iterator fsIt = faceSetsPerBCP_[ibcp].begin(),
                                     fsEnd = faceSetsPerBCP_[ibcp].end();
      
      regionElems_[numVolumes_+ibcp].reserve(numElemsPerRegion_[numVolumes_+ibcp]);
      
      for ( ; fsIt != fsEnd; ++fsIt ) {
        zoneStr.str("");
        zoneStr << "ZN" << fsIt->zone << "/FS" << fsIt->num;
        
        // DO NOT USE the dataset "G/KFCPF" here, because its numbering may
        // not be contiguous!
        ReadLongVector( buffer1, 2*fsIt->numElems, "G/KELPF", zoneStr.str() );
        
        nodeOffset = nodeOffsetPerZone_[fsIt->zone-1];
        
        for ( i=0; i<fsIt->numElems; ++i, ++surfElemIdx ) {
          regionElems_[numVolumes_+ibcp].push_back( surfElemIdx+1 );
          
          volElemIdx = buffer1[2*i] - 1 + elemOffsetPerZone_[fsIt->zone-1];
          std::vector<UInt>::const_iterator connectIt
              = topologyData.begin() + volElemIdx * maxNumElemNodes_;
          
          elemTypes[surfElemIdx] = GetFaceOfElement( elemTypes[volElemIdx],
                                                     buffer1[2*i+1],
                                                     connectIt,
                                                     elConnect );
          numElemNodes = elConnect.size();
          for ( j=0; j<numElemNodes; ++j ) {
            topologyData[surfElemIdx*maxNumElemNodes_+j] = elConnect[j]; 
          }
        }
      }
      
      if (beVerbose) {
        std::cout << "Read connectivity of surface region '"
                  << regionNames_[numVolumes_+ibcp] << "'.\n";
      }
    }
    
    // close .def file
    CloseCFXFile();
  }

  void FileReader_CFX::GetRegionElements(std::vector<UInt> & regionElements,
                                              const UInt regionIdx)
  {
    regionElements = regionElems_[regionIdx];
  }

  template<typename A, typename B>
  void FileReader_CFX::MapZoneResultToVolume( std::vector<A> &volResult,
                                              const std::vector<B> &zoneResult,
                                              const std::set<UInt> &regionNodes,
                                              UInt numDOFs )
  {
    volResult.resize(regionNodes.size()*numDOFs);
    
    std::set<UInt>::iterator nodeIt = regionNodes.begin(),
                             /*beginIt = regionNodes.begin(),*/
                             endIt = regionNodes.end();
    
    for ( UInt i=0; nodeIt != endIt; ++nodeIt, ++i ) {
      for ( UInt j=0; j<numDOFs; ++j ) {
        volResult[i*numDOFs+j] = zoneResult[(*nodeIt-1)*numDOFs+j];
      }
    }    
  }
  
  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CFX::ReadNodalValues(
      std::vector<FlowDataType>& nodalFlowData,
      const std::vector<bool>& activeParts,
      const UInt timeStepIdx )
  {
    Settings& settings = Settings::Instance();
    
    bool beVerbose = settings.GetInt("verbose");
    bool floatDS, retVal;
    
    if(determineFloatDS_) 
    {
      floatDS = true;
    } 
    else {
      floatDS = settings.GetInt("floatDataset");
    }
    
    // Open input file
    if ( beVerbose ) {
      std::cout << "Opening file "<< trnFilenames_[timeStepIdx] << std::endl;
    }
    OpenCFXFile(trnFilenames_[timeStepIdx]);

    for ( UInt actPart=0; actPart < numVolumes_; ++actPart )
    {
      UInt currZone = volumes_[actPart].zone;
      Integer numZoneNodes;
      FlowDataType& fd = nodalFlowData[actPart];
      UInt numDOFs;
      std::ostringstream zonePath;

      if(!activeParts[actPart])
        continue;

      if(beVerbose)
      {
        std::cout << "Reading data on " << GetRegionName(actPart)
                  << std::endl;
      }

      numZoneNodes = ( currZone == numZones_ ) ?
                numNodes_ - nodeOffsetPerZone_[currZone-1] :
                nodeOffsetPerZone_[currZone] - nodeOffsetPerZone_[currZone-1];
      
      zonePath.clear();
      zonePath.str("");
      zonePath << "ZN" << volumes_[actPart].zone << "/VX";
      
      //-----------------------------------------------------------------------
      //     Reading velocity from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_VELOCITY] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        std::string velPath;
        if ( settings.GetInt("cfxUseStnFrame") ) {
          velPath = "G/VELS_FL1";
        } else {
          velPath = "G/VEL_FL1";
        }

        if ( floatDS ) {
          retVal = ReadLongVector( floatvec_, 3*numZoneNodes,
                                   velPath, zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        } else {
          retVal = ReadLongVector( doublevec_, 3*numZoneNodes,
                                   velPath, zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        }
        if (retVal) {
          if (determineFloatDS_) {
            settings.SetInt("floatDataset", 1);
            determineFloatDS_ = false;
          }
        } else {
          if (determineFloatDS_) {
            if ( ReadLongVector( doublevec_, 3*numZoneNodes,
                                 velPath, zonePath.str(),
                                 timeStepNumbers_[timeStepIdx], false ) )
            {
              settings.SetInt("floatDataset", 0);
              floatDS = false;
              determineFloatDS_ = false;
            } else {
              EXCEPTION("Could not determine if CFX datasets are single"
                        " or double precision");
            }
          } else if(beVerbose) {
            std::cerr << "WARNING: CFX dataset does not contain velocity!"
            << std::endl;
          }
        }
      
        FlowDataPartStruct& fdps = fd[FLUIDMECH_VELOCITY];
        fdps.isActive = true; // all partitions have results
        if(fdps.dofNames.empty())
        {

          fdps.definedOn = ResultInfo::NODE; // nodes
          fdps.dofNames.push_back("x");
          fdps.dofNames.push_back("y");
          if(dim_ == 3)
            fdps.dofNames.push_back("z");

          fdps.unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
          fdps.resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
          fdps.entryType = ResultInfo::VECTOR;
        }
        numDOFs = fdps.dofNames.size();

        if ( floatDS ) {
          MapZoneResultToVolume( fdps.data, floatvec_,
                                 volumes_[actPart].nodes, numDOFs );
        } else {
          MapZoneResultToVolume( fdps.data, doublevec_,
                                 volumes_[actPart].nodes, numDOFs );
        }
      }

      //-----------------------------------------------------------------------
      //     Reading pressure from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_PRESSURE] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        if(floatDS) {
          retVal = ReadLongVector( floatvec_, numZoneNodes,
                                   "G/PRES", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        } else {
          retVal = ReadLongVector( doublevec_, numZoneNodes,
                                   "G/PRES", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        }

        if (retVal) {
          FlowDataPartStruct& fdps = fd[FLUIDMECH_PRESSURE];
          fdps.isActive = true; // all partitions have results
          if(fdps.dofNames.empty())
          {
            fdps.definedOn = ResultInfo::NODE; // nodes
            fdps.dofNames.push_back("-");
            fdps.unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
            fdps.resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);
            fdps.entryType = ResultInfo::SCALAR;
          }
          numDOFs = fdps.dofNames.size();

          if ( floatDS ) {
            MapZoneResultToVolume( fdps.data, floatvec_,
                                   volumes_[actPart].nodes, numDOFs );
          } else {
            MapZoneResultToVolume( fdps.data, doublevec_,
                                   volumes_[actPart].nodes, numDOFs );
          }
        } else {
          if (beVerbose) {
            std::cerr << "WARNING: CFX dataset does not contain pressure!\n";
          }
        }
      }

      //-----------------------------------------------------------------------
      //     Reading density from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_DENSITY] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        if (floatDS) {
          retVal = ReadLongVector( floatvec_, numZoneNodes,
                                   "G/DENSITY_FL1", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        } else {
          retVal = ReadLongVector( doublevec_, numZoneNodes,
                                   "G/DENSITY_FL1", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        }

        if (retVal) {
          FlowDataPartStruct& fdps = fd[FLUIDMECH_DENSITY];
          fdps.isActive = true; // all partitions have results
          if(fdps.dofNames.empty())
          {
            fdps.definedOn = ResultInfo::NODE; // nodes
            fdps.dofNames.push_back("-");
            fdps.unit = MapSolTypeToUnit(FLUIDMECH_DENSITY);
            fdps.resultName = SolutionTypeEnum.ToString(FLUIDMECH_DENSITY);
            fdps.entryType = ResultInfo::SCALAR;
          }
          numDOFs = fdps.dofNames.size();

          if ( floatDS ) {
            MapZoneResultToVolume( fdps.data, floatvec_,
                                   volumes_[actPart].nodes, numDOFs );
          } else {
            MapZoneResultToVolume( fdps.data, doublevec_,
                                   volumes_[actPart].nodes, numDOFs );
          }
        } else {
          if (beVerbose) {
            std::cerr << "WARNING: CFX dataset does not contain density!\n";
          }
        }
      }

      //-----------------------------------------------------------------------
      //     Reading turbulent kinetic energy from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_TKE] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        if (floatDS) {
          retVal = ReadLongVector( floatvec_, numZoneNodes,
                                   "G/TKE_FL1", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        } else {
          retVal = ReadLongVector( doublevec_, numZoneNodes,
                                   "G/TKE_FL1", zonePath.str(),
                                   timeStepNumbers_[timeStepIdx], false );
        }

        if (retVal) {
          FlowDataPartStruct& fdps = fd[FLUIDMECH_TKE];
          fdps.isActive = true; // all partitions have results
          if(fdps.dofNames.empty())
          {
            fdps.definedOn = ResultInfo::NODE; // nodes
            fdps.dofNames.push_back("-");
            fdps.unit = MapSolTypeToUnit(FLUIDMECH_TKE);
            fdps.resultName = SolutionTypeEnum.ToString(FLUIDMECH_TKE);
            fdps.entryType = ResultInfo::SCALAR;
          }
          numDOFs = fdps.dofNames.size();

          if ( floatDS ) {
            MapZoneResultToVolume( fdps.data, floatvec_,
                                   volumes_[actPart].nodes, numDOFs );
          } else {
            MapZoneResultToVolume( fdps.data, doublevec_,
                                   volumes_[actPart].nodes, numDOFs );
          }
        } else {
          if (beVerbose) {
            std::cerr << "WARNING: CFX dataset does not contain turb. kin. energy!\n";
          }
        }
      }
      // TEMP_FL1
      // DENSITY
      // DENSITY_FL1
      // TKE_FL1
      // TED_FL1
      // VISCTRB_FL1
      // CONDUCT_FL1
      // SPHEATP_FL1
      // SPHEATV_FL1
      // PTOT
      // PTOTS
      // TTOT_FL1
      // TTOTS_FL1
      // TAUWVEC_FL1
      // YPLUS_FL1
      // QWALL_FL1
      // PRES
      // PSTAT
      // MACH_FL1
      // MACHS_FL1
      // TOTAL_ENTHALPY
      // ENTHTOT_FL1
    }

    //-----------------------------------------------------------------------
    //     Close INPUT file
    //-----------------------------------------------------------------------
    CloseCFXFile();
  }

  void FileReader_CFX::GetInfosFromCommand(std::vector<std::string> commands)
  {
    std::string cmd, attrib;
    int pos=0;
    std::ostringstream sstr;

    rootNode_.reset(new ParamNode(ParamNode::EX, ParamNode::ELEMENT ) );
    rootNode_->SetName("CFX_COMMANDS_DATASET");
    rootNode_->SetValue("Uninteresting");
    currNode_ = rootNode_;
    parents_.push(rootNode_);
    multiLine_ = false;
    
    std::vector<std::string>::iterator comIt = commands.begin(),
                                       endIt = commands.end();
    for ( ; comIt != endIt; ++comIt ) {
      ParseCCLLine(*comIt);
    }

    // std::string ds;
    // rootNode->ToXML(std::cout, 5);
    // std::cout << ds << "\n";

    ParseCommand(charvec_, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec_, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec_, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec_, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec_, pos, cmd, attrib, "", sstr);

    userDataCFX_COMMANDS_ = sstr.str();

#if 0
    if(solTimeUnit != "[s]")
    {
      EXCEPTION("Unrecognized time unit found: "
                << solTimeUnit);
    }

    std::map< std::string, std::string >::iterator it, end;
    std::string s = timeStepStr;

    if(settings.GetInt("verbose"))
    {
      std::cout << "===================================" << std::endl;


      std::cout << "Definition File: " << defFile << std::endl;
      std::cout << "Time Units: " << solTimeUnit << std::endl;
      std::cout << "Timesteps: " << timeStepStr << std::endl << std::endl;


      for(it = exprMap.begin(), end = exprMap.end();
          it != end;
          it++)
      {
        std::cout << "Expression: " << it->first << " = " << it->second << std::endl;
      }

      std::cout << std::endl;
    }

    for(it = exprMap.begin(), end = exprMap.end();
        it != end;
        it++)
    {
      bool match;// = pcrecpp::RE(it->first).GlobalReplace(it->second, &s);

      boost::regex e1(my_expression);

      if(match)
      {
        exprMap.erase(it);
        it = exprMap.begin();
        end = exprMap.end();
      }
    }

    if(settings.GetInt("verbose"))
    {
      std::cout << "===================================" << std::endl;
      std::cout << "regexp: " << regexp << std::endl;
    }


    //        if(!pcrecpp::RE(regexp).Replace("", &s))
    //        {
    //            std::cerr << "Error when_ replacing time unit in timestep string" << std::endl;
    //            exit(1);
    //        }
    /*
      mu::Parser parser;

      parser.SetExpr(params->GetTimeStep());
      try
      {
      timeStep = parser.Eval();
      }
      catch (mu::Parser::exception_type &e)
      {
      std::cout << e.GetMsg() << std::endl;
      exit(1);
      }

      if(timeStep == -1)
      {
      parser.SetExpr(s);

      try
      {
      timeStep = parser.Eval();
      }
      catch (mu::Parser::exception_type &e)
      {
      std::cout << e.GetMsg() << std::endl;
      exit(1);
      }
      }

      if(settings.GetInt("verbose"))
      {
      std::cout << "Timestep string after substitution: " << s << std::endl;
      std::cout << "Timestep: " << timeStep << std::endl;
      }
    */
#endif
  }

  void FileReader_CFX::ParseCCLLine(const std::string& line) 
  {
    typedef boost::char_separator<char> charsep;
    typedef boost::tokenizer< charsep > chartok;
    bool mlSwitch = multiLine_;

    // Check if the current line is a continuation of a previous line.
    if((*line.rbegin()) == '\\') 
    {
      multiLine_ = true;

      mlSwitch = (multiLine_ != mlSwitch);

      if(!mlSwitch) 
      {
        std::stringstream sstr;
        std::string val;
        val = latestNode_->As<std::string>();
        
        sstr << val << line;
        latestNode_->SetValue(sstr.str());

        return;
      }      
    } else 
    {
      if(multiLine_)
      {
        multiLine_ = false;

        std::stringstream sstr;
        std::string val;
        val = latestNode_->As<std::string>();
        
        // Replace continuation backslashes with nothing
        const boost::regex datExp("\\\\");
        const std::string name_format("");
    
        sstr << regex_replace(val, datExp, name_format,
                              boost::match_default | boost::format_sed);

        sstr << line;
        latestNode_->SetValue(sstr.str());
        return;
      }
    }    

    // Cut line into parts at = character. If we have more than one token on a
    // line which  contains a =, then we  just add child nodes  to the current
    // node and we do NOT introduce a new hierarchy level.

    charsep sep("=");
    chartok tok(line, sep);
    
    if(std::distance(tok.begin(), tok.end()) > 1) 
    {
      UInt i = 0;
      std::vector<std::string> keyValue(2);
      
      for(chartok::iterator beg=tok.begin();
          beg!=tok.end(), i < 2;
          ++beg, i++)
      {
        keyValue[i] = boost::algorithm::trim_copy(*beg);
      }

      latestNode_ = 
        PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      
      latestNode_->SetName(keyValue[0]);
      latestNode_->SetValue(keyValue[1]);

      currNode_->AddChildNode( latestNode_ );
    } else 
    {

      // Otherwise cut line into parts at  : character in order to introduce a
      // new hierarchy level

      charsep sep2(":");
      chartok tok2(line, sep2);
      std::string elemName = *tok2.begin();

      if(elemName != "END")
      {

        // If we have  encountered a new element name add a  new child node to
        // the current node.

        latestNode_ = PtrParamNode(new ParamNode(ParamNode::EX,
                                                ParamNode::ELEMENT));
        latestNode_->SetName(elemName);
        latestNode_->SetValue(elemName);

        currNode_->AddChildNode( latestNode_ );
        parents_.push(currNode_);
        currNode_ = latestNode_;

        // If there exists a second token  add it as an attribute to the newly
        // generated node.

        if(std::distance(tok2.begin(), tok2.end()) > 1) 
        {
          chartok::iterator it = tok2.begin();
          it++;

          latestNode_ = 
            PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
          latestNode_->SetName("name");
          latestNode_->SetValue(boost::algorithm::trim_copy(*it));
          currNode_->AddChildNode( latestNode_ );
        }        
      } else 
      {
        // Go one level up in the hierarchy.
        currNode_ = parents_.top();
        parents_.pop();
      }
    }
  }

  void FileReader_CFX::ParseCommand(std::vector<char>& cmdstr,
                                    int& pos,
                                    std::string& cmd,
                                    std::string& attrib,
                                    std::string indent,
                                    std::ostream& outFile)
  {
    int actPos = pos;
    int nextTokenPos;
    char endTest[4];

    while( isprint((int)cmdstr[actPos]) )
    {
      if (cmdstr[actPos] == ':')
      {

        if (cmdstr[actPos+2] == ' ')
        {
          cmdstr[actPos] = 0;
          cmd = &cmdstr[pos];
          attrib = "";
        }
        else
        {
          cmdstr[actPos] = 0;
          cmd = &cmdstr[pos];

          actPos+=2;
          pos = actPos;
          while( !(isspace((int)cmdstr[actPos]) && isspace((int)cmdstr[actPos+1])) )
          {
            actPos++;
          }
          cmdstr[actPos] = 0;
          attrib = &cmdstr[pos];
        }

        break;
      }
      actPos++;
    }

    outFile << indent << cmd << " " << attrib << std::endl;


    do
    {
      actPos++;
    }
    while( isspace((int) cmdstr[actPos]) );


    nextTokenPos = actPos;

    while( isprint((int)cmdstr[actPos]) )
    {
      if (cmdstr[actPos] == ':')
      {
        ind++;
        std::string attrib2;
        ParseCommand(cmdstr, nextTokenPos, cmd, attrib2, indent+"  ", outFile);
        actPos = nextTokenPos;
      }
      else if (cmdstr[actPos] == '=')
      {
        pos = nextTokenPos;
        std::string option, value;
        ParseOption(cmdstr, nextTokenPos, option, value, indent+"    ", outFile);
        actPos = nextTokenPos;

        if(cmd == "RUN DEFINITION" && option == "Definition File")
        {
          defFile_ = value;
        }
        if(cmd == "SOLUTION UNITS" && option == "Time Units")
        {
          solTimeUnit_ = value;
        }
        if(cmd == "TRANSIENT RESULTS" && option == "Time Interval")
        {
          timeStepStr_ = value;
        }

        if(cmd == "EXPRESSIONS")
        {
          exprMap_[option] = value;
        }

      }
      else
      {
        memcpy(endTest, &cmdstr[actPos], 3);
        endTest[3] = 0;

        if(std::string(endTest) == "END")
        {
          actPos +=3;

          outFile << indent << cmd << " " << attrib << std::endl;

          while( isspace((int)cmdstr[actPos]) )
          {
            actPos++;
          }

          pos = actPos;
          return;
        }

        actPos++;
      }

    }

    pos = actPos;
  }


  void FileReader_CFX::ParseOption(std::vector<char>& cmdstr,
                                   int& pos,
                                   std::string& option,
                                   std::string& value,
                                   std::string indent,
                                   std::ostream& outFile)
  {
    Settings& settings = Settings::Instance();
    int actPos = pos;

    while( isprint((int)cmdstr[actPos]) )
    {
      if (cmdstr[actPos] == '=')
      {
        cmdstr[actPos-1] = 0;
        option = &cmdstr[pos];
        pos = actPos+2;
        break;
      }
      actPos++;
    }

    while(isprint((int)cmdstr[actPos]))
    {
      if((cmdstr[actPos] == '\\') && isspace((int)cmdstr[actPos+1]))
      {
        actPos++;
        do
        {
          actPos++;
        }
        while (isspace((int)cmdstr[actPos]));
      }
      else if (isspace((int)cmdstr[actPos]) && isspace((int)cmdstr[actPos+1]))
      {
        break;
      }

      actPos++;
    }

    cmdstr[actPos] = 0;
    value = &cmdstr[pos];

    if(settings.GetInt("verbose"))
    {
      std::cout << indent << option << " = " << value << std::endl;
    }
    outFile << indent << option << " = " << value << std::endl;

    do
    {
      actPos++;
    }
    while( isspace((int)cmdstr[actPos]) );

    pos = actPos;
  }

  void FileReader_CFX::GetUserData(std::map<std::string, std::string>& userData)
  {
    if(userDataCFX_COMMANDS_ == "")
    {
      userDataCFX_COMMANDS_ = "ATTENTION: The COMMANDS dataset has not been read.\n"
                             "           This may be due to a missing .res file!";
    }
    userData["CFX_COMMANDS"] = userDataCFX_COMMANDS_;

    std::ostringstream sstr;

    for(UInt i=0, n=timeStepNumbers_.size(); i<n; i++)
    {
      sstr << timeStepNumbers_[i] << ".trn -> step " << (i+1) << std::endl;
    }

    userData["TRN_TO_STEP_MAP"] = sstr.str();
    
    if (userDataCFXRelease_.length() > 0)
    {
      userData["CFX_Release"] = userDataCFXRelease_;
    }
  }

  void FileReader_CFX::IOErrorToString(int ioerr, std::string& errStr)
  {
    switch(ioerr)
    {
    case __io_ok__:
      errStr = "I/O Operation was OK.";
      break;

    case __io_open_err__:
      errStr = "Failed to open file.";
      break;

    case __io_file_not_open__:
      errStr = "File not open.";
      break;

    case __io_file_not_found__:
      errStr = "File not found.";
      break;

    case __io_disk_full__:
      errStr = "Disk full.";
      break;

    case __io_end_of_file__:
      errStr = "End of file.";
      break;

    case __io_permission_denied__:
      errStr = "Permission denied.";
      break;

    case __io_close_err__:
      errStr = "Failed to close file.";
      break;

    case __io_seek_err__:
      errStr = "Seek error.";
      break;

    case __io_tell_err__:
      errStr = "Tell error.";
      break;

    case __io_read_err__:
      errStr = "Read error.";
      break;

    case __io_write_err__:
      errStr = "Write error.";
      break;

    case __io_parse_err__:
      errStr = "Parse error.";
      break;

    case __io_unknown_dataset__:
      errStr = "Unknown dataset.";
      break;

    case __io_memory_err__:
      errStr = "Memory error.";
      break;

    case __io_file_fmt_err__:
      errStr = "File format error.";
      break;

    case __io_found_more_ds__:
      errStr = "Found more datasets.";
      break;

    case __io_ds_not_found__:
      errStr = "Dataset not found.";
      break;

    case __io_str_err__:
      errStr = "String error.";
      break;

    case __io_mesg_err__:
      errStr = "Message error.";
      break;

    case __io_bad_file__:
      errStr = "Bad file.";
      break;

    case __io_locked__:
      errStr = "Locked file.";
      break;

    case __io_unknown_req__:
      errStr = "Unknown request.";
      break;

    case __io_compress_err__:
      errStr = "Compression error.";
      break;

    case __io_file_not_exist__:
      errStr = "File does not exist.";
      break;

    case __io_attribute_not_found__:
      errStr = "Attribute not found.";
      break;

    case __io_err__:
      errStr = "I/O Error.";
      break;

    case __io_convert_err__:
      errStr = "Convert error.";
      break;

    case __io_memalloc_err__:
      errStr = "Memory allocation error.";
      break;

    default:
      errStr = "Unknown CFX I/O error.";
      break;
    }

  }

  void FileReader_CFX::CheckTransientFiles() 
  {
    Settings& settings = Settings::Instance();
    std::vector< std::string > goodTRNs;
    std::vector< std::string > corruptTRNs;
    std::string lastGoodTRN;
    bool everyThingFine = true;
    uintmax_t meanFileSize = 0;
    uintmax_t actFileSize;
    double fsRatio;
    std::string s, sre;
    boost::regex re;

    std::vector< std::string >::const_iterator it, end;

    /// **************************timeStepNumbers_[timeStepIdx]; ************************

    it = trnFilenames_.begin();
    end = trnFilenames_.end();

    // Check a few .trn files to determine the mean file size
    for( ; it != end; it++ ) 
    {
      actFileSize = fs::file_size(*it);
      meanFileSize = meanFileSize < actFileSize ? actFileSize : meanFileSize;
    }
    
    // Regular expression for matching transient files
    sre = "[0-9]*\\.trn$";
    try
    {
      // Set up the regular expression for case-insensitivity
      re.assign(sre, boost::regex_constants::icase);
    }
    catch (boost::regex_error& e) {}

    it = trnFilenames_.begin();
    end = trnFilenames_.end();

    for( ; it != end; it++ ) 
    {
      // Check the format of the .trn file names
      s = fs::basename(*it) + fs::extension(*it);
      
      if (!boost::regex_match(s, re))
      {
        corruptTRNs.push_back(*it);
        everyThingFine = false;

        continue;
      }
      

      actFileSize = fs::file_size(*it);
      fsRatio = actFileSize / (double) meanFileSize;

      // If the file size diverges by more than 1%
      if(std::fabs(fsRatio - 1.0) > settings.GetDouble("trntol") )
      {
        corruptTRNs.push_back(*it);
        everyThingFine = false;

        // If we want to be fault tolerant we just add the last good .trn
        if(!settings.GetInt("strict"))
        {
          if(lastGoodTRN == "")
            EXCEPTION("Could not find last good .trn file for replacing a corrupt one!");

          goodTRNs.push_back(lastGoodTRN);
        }
      }
      else
      {
        // Check if we really have a .trn file and some .bak or other one
        goodTRNs.push_back(*it);
        lastGoodTRN = *it;
      }
    }

    if(!everyThingFine) 
    {
      std::stringstream sstr;
      
      sstr << "The following .trn files seem to be corrupt:" << std::endl;
      it = corruptTRNs.begin();
      end = corruptTRNs.end();

      for( ; it != end; it++ ) 
      {
        sstr << (*it) << std::endl;
      }

      if(settings.GetInt("strict"))
      {
        EXCEPTION(sstr.str());
      }
      else
      {
        std::cerr << sstr.str() << std::endl;
      }
    }

    trnFilenames_ = goodTRNs;

    it = trnFilenames_.begin();
    end = trnFilenames_.end();

    timeStepNumbers_.clear();
    timeStepNumbers_.reserve(trnFilenames_.size());
    
    for( ; it != end; it++ ) 
    {
      std::stringstream sstr;
      UInt stepnum;
      
      sstr << fs::basename(*it);
      sstr >> stepnum;
      
      timeStepNumbers_.push_back(stepnum);
    }

    numSteps_ = timeStepNumbers_.size();
    
  }
  
  std::string FileReader_CFX::GetRegionName(const UInt regionIdx)
  {
    if ( regionIdx >= numRegions_ ) {
      EXCEPTION("Region index too large.");
    }
    return regionNames_[regionIdx];
  }
  
  UInt FileReader_CFX::GetFaceOfElement( UInt elemType, UInt face,
                           std::vector<UInt>::const_iterator &elemConnect,
                           std::vector<UInt> &faceConnect )
  {
    /* This function fills faceConnect with the nodes of the request face.
     * The element's connectivity must be given by elemConnect.
     * 
     * Face indexes are taken from "ANSYS CFX Reference Guide Rel 13.0,
     * Sec. 3.4.7.1. cfxExportFaceNodes", except for hexahedron (which has to
     * be renumbered for compatibility with CFS++).
     */
    switch ( elemType )
    {
      case Elem::TET4:
        faceConnect.resize(3, 0);
        switch ( face )
        {
          case 1:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+1);
            faceConnect[2] = *(elemConnect+2);
            break;
          case 2:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+3);
            faceConnect[2] = *(elemConnect+1);
            break;
          case 3:
            faceConnect[0] = *(elemConnect+1);
            faceConnect[1] = *(elemConnect+3);
            faceConnect[2] = *(elemConnect+2);
            break;
          case 4:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+2);
            faceConnect[2] = *(elemConnect+3);
            break;
          default:
            EXCEPTION("Invalid face index: " << face);
            break;
        }
        return Elem::TRIA3;
        break;
        
      case Elem::PYRA5:
        faceConnect.resize(3, 0);
        switch ( face )
        {
          case 1:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+3);
            faceConnect[2] = *(elemConnect+4);
            break;
          case 2:
            faceConnect[0] = *(elemConnect+1);
            faceConnect[1] = *(elemConnect+4);
            faceConnect[2] = *(elemConnect+2);
            break;
          case 3:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+4);
            faceConnect[2] = *(elemConnect+1);
            break;
          case 4:
            faceConnect[0] = *(elemConnect+2);
            faceConnect[1] = *(elemConnect+4);
            faceConnect[2] = *(elemConnect+3);
            break;
          case 5:
            faceConnect.resize(4, 0);
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+1);
            faceConnect[2] = *(elemConnect+2);
            faceConnect[3] = *(elemConnect+3);
            return Elem::QUAD4;
            break;
          default:
            EXCEPTION("Invalid face index: " << face);
            break;
        }
        return Elem::TRIA3;
        break;
        
      case Elem::WEDGE6:
        faceConnect.resize(4, 0);
        switch ( face )
        {
          case 1:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+2);
            faceConnect[2] = *(elemConnect+5);
            faceConnect[3] = *(elemConnect+3);
            break;
          case 2:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+3);
            faceConnect[2] = *(elemConnect+4);
            faceConnect[3] = *(elemConnect+1);
            break;
          case 3:
            faceConnect[0] = *(elemConnect+1);
            faceConnect[1] = *(elemConnect+4);
            faceConnect[2] = *(elemConnect+5);
            faceConnect[3] = *(elemConnect+2);
            break;
          case 4:
            faceConnect.resize(3, 0);
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+1);
            faceConnect[2] = *(elemConnect+2);
            return Elem::TRIA3;
            break;
          case 5:
            faceConnect.resize(3, 0);
            faceConnect[0] = *(elemConnect+3);
            faceConnect[1] = *(elemConnect+5);
            faceConnect[2] = *(elemConnect+4);
            return Elem::TRIA3;
            break;
          default:
            EXCEPTION("Invalid face index: " << face);
            break;
        }
        return Elem::QUAD4;
        break;
        
      case Elem::HEXA8:
        faceConnect.resize(4, 0);
        switch ( face )
        {
          /* These indexes are different from those in the CFX Reference,
           * because the connectivity of hexahedra was renumbered before
           * (in function ReadTopology).
           */
          case 1:
            faceConnect[0] = *(elemConnect+4);
            faceConnect[1] = *(elemConnect+5);
            faceConnect[2] = *(elemConnect+1);
            faceConnect[3] = *(elemConnect+0);
            break;
          case 2:
            faceConnect[0] = *(elemConnect+7);
            faceConnect[1] = *(elemConnect+3);
            faceConnect[2] = *(elemConnect+2);
            faceConnect[3] = *(elemConnect+6);
            break;
          case 3:
            faceConnect[0] = *(elemConnect+4);
            faceConnect[1] = *(elemConnect+0);
            faceConnect[2] = *(elemConnect+3);
            faceConnect[3] = *(elemConnect+7);
            break;
          case 4:
            faceConnect[0] = *(elemConnect+5);
            faceConnect[1] = *(elemConnect+6);
            faceConnect[2] = *(elemConnect+2);
            faceConnect[3] = *(elemConnect+1);
            break;
          case 5:
            faceConnect[0] = *(elemConnect+4);
            faceConnect[1] = *(elemConnect+7);
            faceConnect[2] = *(elemConnect+6);
            faceConnect[3] = *(elemConnect+5);
            break;
          case 6:
            faceConnect[0] = *(elemConnect+0);
            faceConnect[1] = *(elemConnect+1);
            faceConnect[2] = *(elemConnect+2);
            faceConnect[3] = *(elemConnect+3);
            break;
          default:
            EXCEPTION("Invalid face index: " << face);
            break;
        }
        return Elem::QUAD4;
        break;
      default:
        EXCEPTION("Element type " << elemType << " not supported");
        break;
    }
    return Elem::UNDEF;
  }

} // end of namespace CoupledField
