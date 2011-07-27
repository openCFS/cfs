// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <map>
#include <set>
#include <fstream>
#include <stdio.h>
#include <iomanip>

#include <boost/regex.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/tokenizer.hpp>
namespace algo=boost::algorithm;

// #include <pcrecpp.h>
// #include <muParser.h>

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
  std::vector<int> FileReader_CFX::intvec;
  std::vector<double> FileReader_CFX::doublevec;
  std::vector<float> FileReader_CFX::floatvec;
  std::vector<char> FileReader_CFX::charvec;

  int ind = 0;

  FileReader_CFX::FileReader_CFX(const std::string& name,
                                 const UInt dim,
                                 const UInt numFiles) :
    FileReader(name, dim, numFiles, 1),
    determineFloatDS_(true),
    numElems_(0)
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

  void FileReader_CFX::Init()
  {
    Settings& settings = Settings::Instance();
    
    std::stringstream sstr;
    sstr << baseName_ << name_ << ".res";
    std::string resFileName = sstr.str();
    std::cout << "resFileName: " << resFileName << std::endl;

    if(fs::exists(resFileName))
    {

      //-----------------------------------------------------------------------
      //     Open RESULTS file
      //-----------------------------------------------------------------------

      snprintf(fn, sizeof(fn), "%s", resFileName.c_str());
      whatfile = __io_open_primaryfile__;
      if(settings.GetInt("verbose"))
      {
        std::cout << "Trying to open CFX results file " << fn
                  << "." << std::endl;
      }

      openfile_(&nerr, fn, &whatfile, strlen(fn));
      CHECK_CFX_IO(nerr);

      //-----------------------------------------------------------------------
      //     Read the latest time step number
      //-----------------------------------------------------------------------

      sprintf(what, "G/TRANSIENT");
      sprintf(where, "ZN1");

      when  = -1;
      length= 3;
      nsize = 1;
      iopt   = __not_stop_if_failed__;
      ioptar = 0;
      intvec.resize(length*nsize);
      floatvec.resize(length*nsize);
      int its;

      redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
              &iopt, &ioptar,
              &floatvec[0],&intvec[0],carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);

      its = intvec[2];

      //-----------------------------------------------------------------------
      //     Read the trninfo data set
      //-----------------------------------------------------------------------
      sprintf(what, "G/TRN_INFO");
      sprintf(where, "EVERY");
      when  = its;

      dattyp    = -3;
      length    = 1;
      nsize     = its;
      iopt   = __stop_if_failed__;
      ioptar = 0;
      int ntrn, dummy;
      floatvec.resize(length*nsize);
      intvec.resize(length*nsize);
      charvec.resize(length*nsize*80);

      redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
              &iopt, &ioptar,
              &floatvec[0],&intvec[0],carr,&dummy,darr,&charvec[0],
              strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);

      ntrn = nsize;

      if ( settings.GetDouble("timestep") <= 0 )
      {
        settings.SetDouble("timestep",
            (floatvec[1]-floatvec[0])/(float)(intvec[1]-intvec[0]));
      }

      numSteps_ = 0;
      transientFNs_.clear();
      timeStepNumbers_.clear();

      for(int i = 0; i< ntrn; i++)
      {
        int its  = intvec[i];

        char* trnnam = &charvec[i*80];
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
          transientFNs_.push_back(fileName);
          timeStepNumbers_.push_back(its);
        }

      }

      //-----------------------------------------------------------------------
      //     read the G/COMMANDS string from the res file
      //     to get the name of the def file and the time step
      //-----------------------------------------------------------------------

      sprintf(what, "G/COMMANDS");
      sprintf(where, "EVERY");

      dattyp = __string_data_type__;
      when  = its;
      length= 1000000;
      nsize = 1;
      iopt   = __stop_if_failed__;
      ioptar = 0;
      charvec.resize(length*nsize);
      std::fill(charvec.begin(), charvec.end(), 0);
      
      redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
              &iopt, &ioptar,
              rarr,iarr,carr,larr,darr,&charvec[0],
              strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);

      whatfile = __io_close_primaryfile__;
      closefile_(&nerr, &whatfile);
      CHECK_CFX_IO(nerr);

      //-----------------------------------------------------------------------
      //     Parse the command string from the CFX results file and
      //     get infos about definition file, time unit and timestep
      //-----------------------------------------------------------------------

      GetInfosFromCommand(nsize);
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
          fn = dir_itr->leaf();

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
      for( ; it != end; it++ )
      {
        stepNum = *it;
        sstr.clear(); sstr.str("");
        sstr << baseName_ << name_ << "/" << stepNum << ".trn";

        fn = sstr.str();

        numSteps_++;
        transientFNs_.push_back(fn);
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
    if(defFile != "")
      defFileNames.push_back(baseName_ + defFile);
    defFileNames.push_back(baseName_ + name_ + ".def");
    defFile = "";

    for(UInt i=0; i<defFileNames.size(); i++)
    {
      if(settings.GetInt("verbose"))
      {
        std::cerr << "Trying to open deffile: " << defFileNames[i] << " ";
      }

      inFile_.clear();
      inFile_.open(defFileNames[i].c_str());
      if (inFile_)
      {
        defFile = defFileNames[i];
        inFile_.close();
        if(settings.GetInt("verbose"))
        {
          std::cerr << "-> OK!" << std::endl;
        }

        break;
      }

      if(settings.GetInt("verbose"))
      {
        std::cerr << "-> failed!" << std::endl;
      }
    }

    if(defFile == "")
    {
      EXCEPTION("Can not find definition file.");
    }

    snprintf(fn, sizeof(fn),"%s", defFile.c_str());
    whatfile = __io_open_primaryfile__;
    openfile_(&nerr,
              fn,
              &whatfile,
              strlen(fn));
    CHECK_CFX_IO(nerr);

    //-----------------------------------------------------------------------
    //     Reading number of vertices: NVX
    //-----------------------------------------------------------------------

    sprintf(what, "G/NVX");
    sprintf(where, "ZN1");
    when  = 0;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = 1;
    iopt   = __stop_if_failed__;
    ioptar = 0;

    nvx = 0;

    redsht_(&dattyp, &length, &nerr,
            what, where, &when,
            &nsize,
            &iopt, &ioptar,
            rarr, &nvx, carr, larr, darr, sarr,
            strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    if(settings.GetInt("verbose"))
    {
      printf("Number of vertices: %d\n", nvx);
    }


    //-----------------------------------------------------------------------
    //     Reading element connectivity
    //-----------------------------------------------------------------------
    //
    //---- Number of element sets
    //
    sprintf(what, "G/NES");
    sprintf(where, "ZN1");
    when  = 0;
    int nes;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = 1;
    iopt   = __stop_if_failed__;
    ioptar = 0;

    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,&nes,carr,larr,darr,sarr, strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    if(settings.GetInt("verbose"))
    {
      printf("Number of element sets, NES= %d\n", nes);
    }

    numVolRegions_ = nes;

    //
    //---- Number of surfaces
    //
    sprintf(what, "G/NSF");
    sprintf(where, "ZN1");
    when  = 0;
    int nsf;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = 1;
    iopt   = __stop_if_failed__;
    ioptar = 0;

    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,&nsf,carr,larr,darr,sarr, strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    if(settings.GetInt("verbose"))
    {
      printf("Number of surface sets, NSF= %d\n", nsf);
    }

    numRegions_ = numVolRegions_ + nsf;
    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);

    //
    //---- Number of elements per set
    //
    sprintf(what, "G/NELES");
    sprintf(where, "ZN1");
    when  = 0;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = numVolRegions_;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    intvec.resize(length*nsize);

    redsht_( &dattyp, &length, &nerr,
             what, where, &when, &nsize,
             &iopt, &ioptar,
             rarr, &intvec[0], carr, larr, darr, sarr,
             strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    for ( UInt i=0; i<numVolRegions_; ++i )
    {
      numElemsPerRegion_[i] = intvec[i];
      numElems_ += numElemsPerRegion_[i];
      numNodesPerRegion_[i] = nvx;
    }
    
    //
    //---- number of surface elements per surface region
    //
    sprintf(what, "G/NFCFS");
    sprintf(where, "ZN1");
    when = 0;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = numRegions_ - numVolRegions_;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    intvec.resize(length*nsize);

    redsht_( &dattyp, &length, &nerr,
             what, where, &when, &nsize,
             &iopt, &ioptar,
             rarr, &intvec[0], carr, larr, darr, sarr,
             strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    for ( UInt i=numVolRegions_; i<numRegions_; ++i )
    {
      numElemsPerRegion_[i] = intvec[i-numVolRegions_];
      numElems_ += numElemsPerRegion_[i];
      numNodesPerRegion_[i] = nvx;
    }
    
    //---- Element type per element set
    //---   elem type = 4: tet  , 4 nodes
    //---             = 5: wedge, 6 nodes
    //---             = 6: hex,   8 nodes
    //---             = 7: pyr,   5 nodes

    sprintf(what, "G/ILTPES");
    sprintf(where, "ZN1");
    when  = 0;

    dattyp = __int_data_type__;
    length = 1;
    nsize  = numVolRegions_;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    intvec.resize(length*nsize);

    redsht_( &dattyp, &length, &nerr,
             what, where, &when, &nsize,
             &iopt, &ioptar,
             rarr, &intvec[0], carr, larr, darr, sarr,
             strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    //
    //---- loop over all element sets
    //
    for(int ies = 1; ies <= (int)numVolRegions_; ies++)
    {
      //
      //---- number of nodes per element
      //
      switch (intvec[ies-1])
      {
        case 4:
          regionElemTypes_.push_back(Elem::TET4);
          break;
        case 5:
          regionElemTypes_.push_back(Elem::WEDGE6);
          break;
        case 6:
          regionElemTypes_.push_back(Elem::HEXA8);
          break;
        case 7:
          regionElemTypes_.push_back(Elem::PYRA5);
          break;
      }

      UInt nENod = Elem::GetNumElemNodes(*regionElemTypes_.rbegin());
      if ( nENod > maxNumElemNodes_ )
        maxNumElemNodes_ = nENod;

      if(settings.GetInt("verbose"))
      {
        printf("ES: %d\nNumber of Elements: %d\nElem type: %d\n",
               ies, nsize, intvec[ies-1]);


        std::cout << "Partition " << (ies)
                  << " nodes: " << numNodesPerRegion_[ies-1]
                  << " elems: " << numElemsPerRegion_[ies-1]
                  << std::endl;
      }
    }

    
    //
    //---- read regions names
    //
    sprintf(what, "G/CLBVL");
    sprintf(where, "ZN1");

    dattyp = __string_data_type__;
    when   = 0;
    length = numVolRegions_;
    nsize  = 80;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    charvec.resize(length*nsize);
    
    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,iarr,carr,larr,darr,&charvec[0],
            strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    regionNames_.clear();
    std::ostringstream oss;
    
    for ( UInt iRegion=0; iRegion<numVolRegions_; ++iRegion )
    {
      char *regionName = &charvec[iRegion*80];
      
      for ( UInt n=79; n>0; --n )
      {
        if ( isblank(regionName[n]) )
          regionName[n] = 0;
        else
          break;
      }
      
      oss.str("");
      oss << regionName;
      regionNames_.push_back(oss.str());
      
      std::replace( regionNames_[iRegion].begin(),
                    regionNames_[iRegion].end(),
                    ' ', '_' );
    }

    
    //
    //---- read surface regions names
    //
    UInt numSurfRegions = numRegions_ - numVolRegions_;
    
    sprintf(what, "G/CLBSF");
    sprintf(where, "ZN1");

    dattyp = __string_data_type__;
    when   = 0;
    length = numSurfRegions;
    nsize  = 80;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    charvec.resize(length*nsize);
    
    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,iarr,carr,larr,darr,&charvec[0],
            strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    
    for ( UInt iRegion=0, regionIdx=numVolRegions_;
          iRegion<numSurfRegions;
          ++iRegion, ++regionIdx )
    {
      char *regionName = &charvec[iRegion*80];
      
      for ( UInt n=79; n>0; --n )
      {
        if ( isblank(regionName[n]) )
          regionName[n] = 0;
        else
          break;
      }
      
      oss.str("");
      oss << regionName;
      regionNames_.push_back(oss.str());
      
      std::replace( regionNames_[regionIdx].begin(),
                    regionNames_[regionIdx].end(),
                    ' ', '_' );
    }
    
    
    //
    //---- read CFX Release No. for UserData
    //
    sprintf(what, "G/CRELNO");
    sprintf(where, "EVERY");

    dattyp = __string_data_type__;
    when  = 0;
    length= 1;
    nsize = 101;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    charvec.resize(length*nsize);
    
    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,iarr,&charvec[0],larr,darr,carr,
            strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    oss.str("");
    oss << &charvec[0];
    userDataCFXRelease = oss.str();

    
    //
    //---- close file
    //
    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
    
    charvec.clear();
    doublevec.clear();
    floatvec.clear();
    intvec.clear();
  }

  void FileReader_CFX::ReadNodalCoords(std::vector<Double> & NODECOORD)
  {
    Settings& settings = Settings::Instance();
    NODECOORD.resize(numNodesPerRegion_[0]*3);

    snprintf(fn, sizeof(fn),"%s", defFile.c_str());
    whatfile = __io_open_primaryfile__;
    openfile_(&nerr, fn, &whatfile, strlen(fn));
    CHECK_CFX_IO(nerr);

    if(settings.GetInt("verbose"))
    {
      printf("Successfully opened %s\n", defFile.c_str());
    }

    //-----------------------------------------------------------------------
    //     Reading grid coordinates CRDVX (double precision in DEF file)
    //-----------------------------------------------------------------------

    sprintf(what, "G/CRDVX");
    sprintf(where, "ZN1/VX");
    when  = 0;

    dattyp = __double_data_type__;
    length = 3;
    nsize  = numNodesPerRegion_[0];
    iopt = __stop_if_failed__;
    doublevec.resize(numNodesPerRegion_[0]*3, 0.0);

    readlong_(&dattyp, &nerr, what,where,&when,&nsize,&iopt,
              rarr,iarr,carr,larr,&doublevec[0],sarr,
              strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    NODECOORD = doublevec;

    if(settings.GetInt("verbose"))
    {
      printf("Coordinate size = %d\n", nsize);
    }

    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
    
    doublevec.clear();
  }

  void FileReader_CFX::ReadTopology( std::vector<UInt> & TOPOLOGYDATA,
                                         std::vector<UInt> & elemTypes)
  {
    Settings& settings = Settings::Instance();
    UInt elem = 0;
    UInt elemType = Elem::UNDEF;
    UInt numRegionElems=0;
    UInt numElemNodes;
    std::vector<UInt> elConnect(maxNumElemNodes_);

    // open .def file
    snprintf(fn, sizeof(fn),"%s", defFile.c_str());
    whatfile = __io_open_primaryfile__;
    openfile_(&nerr, fn, &whatfile, strlen(fn));
    CHECK_CFX_IO(nerr);

    // allocate memory
    TOPOLOGYDATA.resize(numElems_ * maxNumElemNodes_, 0);
    elemTypes.resize(numElems_, Elem::UNDEF);

    // first read volume regions
    for ( UInt actRegion=0; actRegion<numVolRegions_; ++actRegion )
    {
      elemType = regionElemTypes_[actRegion];
      numRegionElems = numElemsPerRegion_[actRegion];
      numElemNodes = Elem::GetNumElemNodes((Elem::FEType)elemType);

      // read element numbers
      sprintf(what, "G/KELPE");
      sprintf(where, "ZN1/ES%d", actRegion+1);
      when = 0;

      dattyp = __int_data_type__;
      length = numRegionElems;
      nsize  = 1;
      iopt   = __stop_if_failed__;
      intvec.resize(length*nsize);

      readlong_( &dattyp, &nerr, what, where, &when, &nsize, &iopt,
                 rarr, &intvec[0], carr, larr, darr, sarr,
                 strlen(what), strlen(where), 0 );
      CHECK_CFX_IO(nerr);

      regionElems_[actRegion].resize(numRegionElems, 0);
      std::copy( intvec.begin(), intvec.end(),
                 regionElems_[actRegion].begin() );
      
      std::vector<int>::const_iterator it = intvec.begin(),
          itEnd = intvec.end();
      for ( ; it != itEnd; ++it ) {
        elemTypes[(*it)-1] = elemType;
      }
      
      //---- reading connectivity for each element set
      //     outer loop:  i_element
      //     inner loop:  i_vx per element (1 ... nelvx)
      sprintf(what,"G/KVXPE");

      //
      //---- where = ZN1/ESn where n is integer from 1 to nes
      //
      sprintf(where, "ZN1/ES%d", actRegion+1);
      when  = 0;
      dattyp = __int_data_type__;
      length = numElemNodes;
      nsize  = numRegionElems;
      iopt   = __stop_if_failed__;
      intvec.resize(length*nsize);

      readlong_( &dattyp, &nerr, what, where, &when, &nsize, &iopt,
                 rarr, &intvec[0], carr, larr, darr, sarr,
                 strlen(what), strlen(where), 0 );
      CHECK_CFX_IO(nerr);

      std::cout  << "Read connectivity from .def file" << std::endl;

      if ( settings.GetInt("verbose") )
      {
        printf("Length of Connectivity array: %d\n", nsize);
      }

      UInt baseIdx=0;
      for ( UInt i=0; i<numRegionElems; ++i, baseIdx += numElemNodes )
      {
        std::fill(elConnect.begin(), elConnect.end(), 0);

        if ( elemType == Elem::HEXA8 )
        {
          elConnect[0] = intvec[baseIdx + 4];
          elConnect[1] = intvec[baseIdx + 6];
          elConnect[2] = intvec[baseIdx + 7];
          elConnect[3] = intvec[baseIdx + 5];
          elConnect[4] = intvec[baseIdx + 0];
          elConnect[5] = intvec[baseIdx + 2];
          elConnect[6] = intvec[baseIdx + 3];
          elConnect[7] = intvec[baseIdx + 1];
        }
        else
        {
          std::copy( &intvec[baseIdx],
                     &intvec[baseIdx+numElemNodes],
                     &elConnect[0] );
        }

        std::copy(elConnect.begin(), elConnect.end(),
                  TOPOLOGYDATA.begin() + 
                  (regionElems_[actRegion][i]-1)*maxNumElemNodes_);
        ++elem;
      }

    }

    UInt numSurfRegions = numRegions_ - numVolRegions_;
    
    for ( UInt iRegion=0; iRegion<numSurfRegions; ++iRegion )
    {
      numRegionElems = numElemsPerRegion_[iRegion+numVolRegions_];
      
      sprintf(what,"G/KELPF");
      sprintf(where, "ZN1/FS%d", iRegion+1);
      when  = 0;
      
      dattyp = __int_data_type__;
      length = numRegionElems;
      nsize  = 2;
      iopt   = __stop_if_failed__;
      intvec.resize(length*nsize);

      readlong_( &dattyp, &nerr, what, where, &when, &nsize, &iopt,
                 rarr, &intvec[0], carr, larr, darr, sarr,
                 strlen(what), strlen(where), 0 );
      CHECK_CFX_IO(nerr);
      
      for ( UInt i=0; i<numRegionElems; ++i )
      {
        std::vector<UInt>::const_iterator connectIt
            = TOPOLOGYDATA.begin() + (intvec[2*i]-1)*maxNumElemNodes_;
        
        elemTypes[elem] = GetFaceOfElement( (UInt)elemTypes[intvec[2*i]-1],
                                            (UInt)intvec[2*i+1],
                                            connectIt,
                                            elConnect);
        
        std::copy(elConnect.begin(), elConnect.end(),
                  TOPOLOGYDATA.begin() + elem*maxNumElemNodes_);

        regionElems_[iRegion+numVolRegions_].push_back(++elem);

      }
    }

    // close .def file
    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
    
    intvec.clear();
  }

  void FileReader_CFX::GetRegionElements(std::vector<UInt> & regionElements,
                                              const UInt regionIdx)
  {
    regionElements = regionElems_[regionIdx];
  }

  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CFX::ReadNodalValues(
      std::vector<FlowDataType>& nodalFlowData,
      const std::vector<bool>& activeParts,
      const UInt timeStepIdx )
  {
    Settings& settings = Settings::Instance();
    bool floatDS;
    
    if(determineFloatDS_) 
    {
      floatDS = true;
    } 
    else {
      floatDS = settings.GetInt("floatDataset");
    }
    
    // Open input file
    snprintf(fn, sizeof(fn), "%s",
             transientFNs_[timeStepIdx].c_str());
    whatfile = __io_open_primaryfile__;

    if(settings.GetInt("verbose"))
    {
      std::cout << "Opening file "<< fn << std::endl;
    }

    openfile_(&nerr, fn, &whatfile, strlen(fn));
    CHECK_CFX_IO(nerr);

    for ( UInt actPart=0; actPart < numVolRegions_; ++actPart )
    {
      int nvx = numNodesPerRegion_[actPart];
      FlowDataType& fd = nodalFlowData[actPart];
      UInt numDOFs;

      if(!activeParts[actPart])
        continue;

      if(settings.GetInt("verbose"))
      {
        std::cout << "Reading data on " << GetRegionName(actPart)
                  << std::endl;
      }

      //-----------------------------------------------------------------------
      //     Reading velocity from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_VELOCITY] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        sprintf(what, "G/VEL_FL1");
        sprintf(where, "ZN1/VX");
        when  = timeStepNumbers_[timeStepIdx];

        if (floatDS)
          dattyp = __real_data_type__;
        else
          dattyp = __double_data_type__;
        
        length = 3;
        nsize  = nvx;
        doublevec.resize(length*nsize);
        floatvec.resize(length*nsize);

        iopt   = __stop_if_failed__;

        if(floatDS)
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              &floatvec[0],iarr,carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
        else
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              rarr,iarr,carr,larr,&doublevec[0],sarr,
              strlen(what), strlen(where), 0);

        if(nerr)
        {
          if(determineFloatDS_) 
          {
            dattyp = __double_data_type__;

            readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                      &floatvec[0],iarr,carr,larr,&doublevec[0],sarr,
                      strlen(what), strlen(where), 0);

            if(nerr)
            {
              EXCEPTION("Could not determine if CFX datasets are single"
                        " or double precision");
            }

            settings.SetInt("floatDataset", 0);
            floatDS = false;
            determineFloatDS_ = false;
          }
          else if(settings.GetInt("verbose")) {
            std::cerr << "WARNING: CFX dataset does not contain velocity!"
            << std::endl;
          }
        }
        else 
        {
          if(determineFloatDS_) 
          {
            settings.SetInt("floatDataset", 1);
            determineFloatDS_ = false;
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
        fdps.data.resize(numDOFs * nvx);

        if(floatDS)
          std::copy(floatvec.begin(),
                    floatvec.begin() + (numDOFs * nvx),
                    fdps.data.begin());
        else
          std::copy(doublevec.begin(),
                    doublevec.begin() + (numDOFs * nvx),
                    fdps.data.begin());
      }

      //-----------------------------------------------------------------------
      //     Reading pressure from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_PRESSURE] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        sprintf(what, "G/PRES");
        sprintf(where, "ZN1/VX");
        when  = timeStepNumbers_[timeStepIdx];

        if(floatDS)
          dattyp = __real_data_type__;
        else
          dattyp = __double_data_type__;

        length = 1;
        nsize  = nvx;
        doublevec.resize(length*nsize);
        floatvec.resize(length*nsize);

        iopt   = __stop_if_failed__;

        if(floatDS)
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              &floatvec[0],iarr,carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
        else
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              rarr,iarr,carr,larr,&doublevec[0],sarr,
              strlen(what), strlen(where), 0);

        if(nerr)
        {
          if(settings.GetInt("verbose"))
            std::cerr << "WARNING: CFX dataset does not contain pressure!"
            << std::endl;
        }
        else
        {
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
          fdps.data.resize(numDOFs * nvx);

          if(floatDS)
            std::copy(floatvec.begin(),
                floatvec.begin() + (numDOFs * nvx),
                fdps.data.begin());
          else
            std::copy(doublevec.begin(),
                doublevec.begin() + (numDOFs * nvx),
                fdps.data.begin());
        }
      }

      //-----------------------------------------------------------------------
      //     Reading density from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_DENSITY] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        sprintf(what, "G/DENSITY_FL1");
        sprintf(where, "ZN1/VX");
        when  = timeStepNumbers_[timeStepIdx];

        if(floatDS)
          dattyp = __real_data_type__;
        else
          dattyp = __double_data_type__;

        length = 1;
        nsize  = nvx;
        doublevec.resize(length*nsize);
        floatvec.resize(length*nsize);

        iopt   = __stop_if_failed__;

        if(floatDS)
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              &floatvec[0],iarr,carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
        else
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              rarr,iarr,carr,larr,&doublevec[0],sarr,
              strlen(what), strlen(where), 0);

        if(nerr)
        {
          if(settings.GetInt("verbose"))
            std::cerr << "WARNING: CFX dataset does not contain density!"
            << std::endl;
        }
        else
        {
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
          fdps.data.resize(numDOFs * nvx);

          if(floatDS)
            std::copy(floatvec.begin(),
                floatvec.begin() + (numDOFs * nvx),
                fdps.data.begin());
          else
            std::copy(doublevec.begin(),
                doublevec.begin() + (numDOFs * nvx),
                fdps.data.begin());
        }
      }

      //-----------------------------------------------------------------------
      //     Reading turbulent kinetic energy from input file
      //-----------------------------------------------------------------------
      if(requiredResults_[FLUIDMECH_TKE] ||
         requiredResults_[NO_SOLUTION_TYPE])
      {
        sprintf(what, "G/TKE_FL1");
        sprintf(where, "ZN1/VX");
        when  = timeStepNumbers_[timeStepIdx];

        if(floatDS)
          dattyp = __real_data_type__;
        else
          dattyp = __double_data_type__;

        length = 1;
        nsize  = nvx;
        doublevec.resize(length*nsize);
        floatvec.resize(length*nsize);

        iopt   = __stop_if_failed__;

        if(floatDS)
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              &floatvec[0],iarr,carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
        else
          readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              rarr,iarr,carr,larr,&doublevec[0],sarr,
              strlen(what), strlen(where), 0);

        if(nerr)
        {
          if(settings.GetInt("verbose"))
            std::cerr << "WARNING: CFX dataset does not contain turb. kin. energy!"
            << std::endl;
        }
        else
        {
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
          fdps.data.resize(numDOFs * nvx);

          if(floatDS)
            std::copy(floatvec.begin(),
                floatvec.begin() + (numDOFs * nvx),
                fdps.data.begin());
          else
            std::copy(doublevec.begin(),
                doublevec.begin() + (numDOFs * nvx),
                fdps.data.begin());
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
    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
  }

  void FileReader_CFX::GetInfosFromCommand(UInt numLines)
  {
    std::string cmd, attrib;
    int pos=0;
    std::ostringstream sstr;

    //    std::cout << charvec << std::endl;
    
    std::stringbuf *pbuf;
    std::stringstream ss;
    
    pbuf=ss.rdbuf();
    pbuf->sputn (&charvec[0], charvec.size());
                 //numLines*80);
    //    std::cout << pbuf->str();

    rootNode.reset(new ParamNode(ParamNode::EX, ParamNode::ELEMENT ) );
    rootNode->SetName("CFX_COMMANDS_DATASET");
    rootNode->SetValue("Uninteresting");
    currNode = rootNode;
    parents.push(rootNode);
    multiLine = false;
    
    // Split COMMANDS dataset into 80 character tokens
    std::string str = pbuf->str();
    int offsets[] = {80};
    boost::offset_separator f(offsets, offsets+1,true,false);
    boost::tokenizer<boost::offset_separator> tok(str,f);
    UInt i = 0;
    for(boost::tokenizer<boost::offset_separator>::iterator beg=tok.begin();
        beg!=tok.end(), i < numLines;
        ++beg, i++)
    {
      std::string trimmed = boost::algorithm::trim_copy(*beg);
      ParseCCLLine(trimmed);
      // std::cout << trimmed << "\n";
      // std::cout << *beg << "\n";
    }

    // std::string ds;
    // rootNode->ToXML(std::cout, 5);
    // std::cout << ds << "\n";

    ParseCommand(charvec, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec, pos, cmd, attrib, "", sstr);
    ParseCommand(charvec, pos, cmd, attrib, "", sstr);

    userDataCFX_COMMANDS = sstr.str();

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
    //            std::cerr << "Error when replacing time unit in timestep string" << std::endl;
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
    bool mlSwitch = multiLine;

    // Check if the current line is a continuation of a previous line.
    if((*line.rbegin()) == '\\') 
    {
      multiLine = true;

      mlSwitch = (multiLine != mlSwitch);

      if(!mlSwitch) 
      {
        std::stringstream sstr;
        std::string val;
        val = latestNode->As<std::string>();
        
        sstr << val << line;
        latestNode->SetValue(sstr.str());

        return;
      }      
    } else 
    {
      if(multiLine)
      {
        multiLine = false;

        std::stringstream sstr;
        std::string val;
        val = latestNode->As<std::string>();
        
        // Replace continuation backslashes with nothing
        const boost::regex datExp("\\\\");
        const std::string name_format("");
    
        sstr << regex_replace(val, datExp, name_format,
                              boost::match_default | boost::format_sed);

        sstr << line;
        latestNode->SetValue(sstr.str());
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

      latestNode = 
        PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      
      latestNode->SetName(keyValue[0]);
      latestNode->SetValue(keyValue[1]);

      currNode->AddChildNode( latestNode );
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

        latestNode = PtrParamNode(new ParamNode(ParamNode::EX,
                                                ParamNode::ELEMENT));
        latestNode->SetName(elemName);
        latestNode->SetValue(elemName);

        currNode->AddChildNode( latestNode );
        parents.push(currNode);
        currNode = latestNode;

        // If there exists a second token  add it as an attribute to the newly
        // generated node.

        if(std::distance(tok2.begin(), tok2.end()) > 1) 
        {
          chartok::iterator it = tok2.begin();
          it++;

          latestNode = 
            PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
          latestNode->SetName("name");
          latestNode->SetValue(boost::algorithm::trim_copy(*it));
          currNode->AddChildNode( latestNode );
        }        
      } else 
      {
        // Go one level up in the hierarchy.
        currNode = parents.top();
        parents.pop();
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
    Settings& settings = Settings::Instance();
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

    if(settings.GetInt("verbose"))
    {
      std::cout << indent << cmd << " " << attrib << std::endl;
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
          defFile = value;
        }
        if(cmd == "SOLUTION UNITS" && option == "Time Units")
        {
          solTimeUnit = value;
        }
        if(cmd == "TRANSIENT RESULTS" && option == "Time Interval")
        {
          timeStepStr = value;
        }

        if(cmd == "EXPRESSIONS")
        {
          exprMap[option] = value;
        }

      }
      else
      {
        memcpy(endTest, &cmdstr[actPos], 3);
        endTest[3] = 0;

        if(std::string(endTest) == "END")
        {
          actPos +=3;

          if(settings.GetInt("verbose"))
          {
            std::cout << indent << "END" << std::endl;
          }
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
    if(userDataCFX_COMMANDS == "")
    {
      userDataCFX_COMMANDS = "ATTENTION: The COMMANDS dataset has not been read.\n"
                             "           This may be due to a missing .res file!";
    }
    userData["CFX_COMMANDS"] = userDataCFX_COMMANDS;

    std::ostringstream sstr;

    for(UInt i=0, n=timeStepNumbers_.size(); i<n; i++)
    {
      sstr << timeStepNumbers_[i] << ".trn -> step " << (i+1) << std::endl;
    }

    userData["TRN_TO_STEP_MAP"] = sstr.str();
    
    if (userDataCFXRelease.length() > 0)
    {
      userData["CFX_Release"] = userDataCFXRelease;
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

    it = transientFNs_.begin();
    end = transientFNs_.end();

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

    it = transientFNs_.begin();
    end = transientFNs_.end();

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

    transientFNs_ = goodTRNs;

    it = transientFNs_.begin();
    end = transientFNs_.end();

    timeStepNumbers_.clear();
    
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
    if ( regionIdx >= numRegions_ )
      EXCEPTION("Region index too large.");
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
        }
        return Elem::QUAD4;
        break;
      default:
        EXCEPTION("Element type " << elemType << " not supported");
    }
  }

}
