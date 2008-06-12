#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <stdio.h>
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs=boost::filesystem;
#include <boost/algorithm/string/predicate.hpp>
namespace algo=boost::algorithm;

// #include <pcrecpp.h>
// #include <muParser.h>

#include "../params.hh"
#include "../settings.hh"
#include "../mpcci_defs.hh"
#include "filereader_CFX.hh"
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
  const int FileReader_CFX::BIGMEM = 10000000;
  std::vector<int> FileReader_CFX::intvec;
  std::vector<double> FileReader_CFX::doublevec;
  std::vector<float> FileReader_CFX::floatvec;
  std::vector<char> FileReader_CFX::charvec;

  int ind = 0;

  FileReader_CFX::FileReader_CFX(const std::string& name,
                                 const UInt dim,
                                 const UInt numFiles) :
    FileReader(name, dim, numFiles)
  {
    name_ = name;
    baseName_ = "./";
    baseName_+= name_;
    baseName_+= "/";
  }

  FileReader_CFX::~FileReader_CFX()
  {
  }

  void FileReader_CFX::Init()
  {
    Settings& settings = Settings::Instance();

    if(settings.GetDouble("timeStep") < 0)
      EXCEPTION("No proper time step has been specified! Use --timestep X.");
        
    std::stringstream sstr;
    sstr << baseName_ << name_ << ".res";
    std::string resFileName = sstr.str();
    std::cout << "resFileName: " << resFileName << std::endl;

    intvec.resize(BIGMEM*3);
    floatvec.resize(BIGMEM);

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
      length= 6;
      nsize = 1;
      iopt   = __not_stop_if_failed__;
      ioptar = 0;
      intvec.resize(BIGMEM*3);
      floatvec.resize(BIGMEM);
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
      length    = 3;
      nsize     = 10000;
      iopt   = __stop_if_failed__;
      ioptar = 0;
      int ntrn, dummy;
      charvec.resize(BIGMEM*3);

      redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
              &iopt, &ioptar,
              &floatvec[0],&intvec[0],carr,&dummy,&doublevec[0],&charvec[0],
              strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);

      ntrn = nsize;

      //    printf("sarrt %s, ntrn %d\n", sarrt, ntrn);
    
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
        
        snprintf(fn, sizeof(fn), "%s%s/%s", baseName_.c_str(), name_.c_str(), trnnam);
            
        inFile_.clear();
        inFile_.open(fn);
        if (inFile_)
        {
          inFile_.close();

          numSteps_++;
          transientFNs_.push_back(fn);
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
      length= BIGMEM;
      nsize = 1;
      iopt   = __stop_if_failed__;
      ioptar = 0;

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

      GetInfosFromCommand();
    }
    else 
    {
      fs::path trnDir( baseName_ +  name_);
      fs::directory_iterator end_iter;
      UInt its = 1;
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

          //          if(fn.rfind(".trn") == (fn.length() - 4)) // endswith
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

    if (settings.GetInt("numSteps"))
    {
      UInt tmp = (UInt) settings.GetInt("numSteps");
      /* only take argument if tmo does not exceed the maximal number of timesteps possible */
      if (tmp < numSteps_)
      {
        numSteps_ = tmp;
      }
    }
        
    //-----------------------------------------------------------------------
    //     Open DEFINITION file
    //-----------------------------------------------------------------------

    std::vector< std::string > defFileNames;
    if(settings.GetString("defFile") != "")
      defFileNames.push_back(baseName_ + settings.GetString("defFile").c_str());
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
    length = 3;
    nsize  = 1;
    iopt   = __stop_if_failed__;
    ioptar = 0;
    
    nvx = 0;

    redsht_(&dattyp, &length, &nerr,
            what, where, &when,
            &nsize,
            &iopt, &ioptar,
            rarr, &nvx, carr,
            larr, darr, sarr,
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
        
    numRegions_ = nes;
    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);
        
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
    nsize  = nes;
    iopt   = __stop_if_failed__;
    ioptar = 0;

    redsht_(&dattyp,&length,&nerr,what,where,&when,&nsize,
            &iopt, &ioptar,
            rarr,&intvec[0],carr,larr,darr,sarr, strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    //
    //---- loop over all element sets
    //
    std::vector<int> ielem;
    for(int ies = 0; ies < nes; ies++)
    {
      ielem.push_back(intvec[ies]);
    }

    for(int ies = 1; ies <= nes; ies++)
    {
      //
      //---- number of nodes per element
      //
      if (ielem[ies-1] == 4) regionElemTypes_.push_back(ET_TET4);
      if (ielem[ies-1] == 5) regionElemTypes_.push_back(ET_WEDGE6);
      if (ielem[ies-1] == 6) regionElemTypes_.push_back(ET_HEXA8);
      if (ielem[ies-1] == 7) regionElemTypes_.push_back(ET_PYRA5);

      UInt nENod = NUM_ELEM_NODES[*regionElemTypes_.rbegin()];
      maxNumElemNodes_ = nENod > maxNumElemNodes_ ? nENod : maxNumElemNodes_;
      
      //
      //---- reading element numbers for each element set
      //
      sprintf(what,"G/KELPE");
    
      //
      //---- where = ZN1/ESn where n is integer from 1 to nes
      //
      sprintf(where, "ZN1/ES%d", ies);
      when  = 0;

      dattyp = __int_data_type__;
      length = 1;
      nsize  = BIGMEM;
      iopt   = __stop_if_failed__;
        
      readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                rarr,&intvec[0],carr,larr,darr,sarr,
                strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);


      numElemsPerRegion_[ies-1] = nsize;
      numNodesPerRegion_[ies-1] = nvx;

      if(settings.GetInt("verbose"))
      {
        printf("ES: %d\nNumber of Elements: %d\nElem type: %d\n",
               ies, nsize, ielem[ies-1]);


        std::cout << "Partition " << (ies)
                  << " nodes: " << numNodesPerRegion_[ies-1]
                  << " elems: " << numElemsPerRegion_[ies-1]
                  << std::endl;
      }
    }

    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
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
      printf("Successfully openend %s\n", defFile.c_str());
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
    doublevec.resize(numNodesPerRegion_[0]*3);
    
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
  }
    
  void FileReader_CFX::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                        std::vector<UInt> & elemTypes)
  {
    Settings& settings = Settings::Instance();
    UInt elem=0;
    int numElems=0;
    int numElemNodes;
    int elemType = ET_UNDEF;
    std::vector<UInt> elConnect(maxNumElemNodes_);

    snprintf(fn, sizeof(fn),"%s", defFile.c_str());
    whatfile = __io_open_primaryfile__;
    openfile_(&nerr, fn, &whatfile, strlen(fn)); 
    CHECK_CFX_IO(nerr);

    // Determine total number of elements
    for(UInt actRegion=0; actRegion<numRegions_; actRegion++) {
      numElems += numElemsPerRegion_[actRegion];
    }
    
    TOPOLOGYDATA.resize(numElems * maxNumElemNodes_);

    for(UInt actRegion=0; actRegion<numRegions_; actRegion++) {
      elemType = regionElemTypes_[actRegion];
      numElems = numElemsPerRegion_[actRegion];
      numElemNodes = NUM_ELEM_NODES[elemType];

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
      nsize  = BIGMEM;
      iopt   = __stop_if_failed__;


      readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
          rarr,&intvec[0],carr,larr,darr,sarr,
          strlen(what), strlen(where), 0);
      CHECK_CFX_IO(nerr);

      std::cout  << "read connectivity from def file" << std::endl;

      if(settings.GetInt("verbose"))
      {
        printf("Length of Connectivity array: %d\n", nsize);
      }

      UInt baseIdx=0;
      for(int i=0; i<numElems; i++, baseIdx += numElemNodes) 
      {
        elemTypes.push_back(elemType);
        std::fill(elConnect.begin(), elConnect.end(), 0);

        if(elemType == ET_HEXA8)
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
          std::copy(&intvec[baseIdx],
              &intvec[baseIdx+numElemNodes],
              &elConnect[0]);
        }

        regionElems_[actRegion].push_back(elem+1);

        std::copy(elConnect.begin(), elConnect.end(),
                  TOPOLOGYDATA.begin() + elem*maxNumElemNodes_);
        elem++;
        
      }
      
    }

    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
  }
    
  void FileReader_CFX::GetRegionElements(std::vector<UInt> & regionElements,
                                              const UInt regionIdx)
  {
    regionElements = regionElems_[regionIdx];
  }

  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CFX::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                       const std::vector<bool>& activeParts,
                                       const UInt timeStepIdx)
  {
    Settings& settings = Settings::Instance();
    bool floatDS = settings.GetInt("floatDataset");
    
    // Open input file
    snprintf(fn, 
             sizeof(fn),
             "%s",
             transientFNs_[timeStepIdx].c_str());
    whatfile = __io_open_primaryfile__;

    if(settings.GetInt("verbose"))
    {
      std::cout << "Opening file "<< fn << std::endl;
    }
    
    openfile_(&nerr, fn, &whatfile, strlen(fn)); 
    CHECK_CFX_IO(nerr);

    for(UInt actPart=0; actPart < numRegions_; actPart++)
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

      sprintf(what, "G/VEL_FL1");
      sprintf(where, "ZN1/VX");
      when  = timeStepNumbers_[timeStepIdx];

      if(floatDS)
        dattyp = __real_data_type__;
      else
        dattyp = __double_data_type__;

      length = 3;
      nsize  = nvx;
      iopt   = __stop_if_failed__;
        
      if(floatDS)
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,darr,sarr,
                  strlen(what), strlen(where), 0);
      else
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,&doublevec[0],sarr,
                  strlen(what), strlen(where), 0);

      if(nerr)                                                  
      {                                                         
        if(settings.GetInt("verbose"))
          std::cerr << "WARNING: CFX dataset does not contain velocity!"
                    << std::endl;                           
      }
      else 
      {
        FlowDataPartStruct& fdps = fd[FLUIDMECH_VELOCITY];
        fdps.isActive = true; // all partitions have results
        fdps.definedOn = ResultInfo::NODE; // nodes
        fdps.dofNames.push_back("x");
        fdps.dofNames.push_back("y");
        if(dim_ == 3) 
          fdps.dofNames.push_back("z");
        fdps.unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
        Enum2String(FLUIDMECH_VELOCITY, fdps.resultName);
        numDOFs = fdps.dofNames.size();
        fdps.data.resize(numDOFs * nvx);
        fdps.entryType = ResultInfo::VECTOR;
        
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

      sprintf(what, "G/PRES");
      sprintf(where, "ZN1/VX");
      when  = timeStepNumbers_[timeStepIdx];

      if(floatDS)
        dattyp = __real_data_type__;
      else
        dattyp = __double_data_type__;

      length = 1;
      nsize  = nvx;
      iopt   = __stop_if_failed__;
        
      if(floatDS)
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,darr,sarr,
                  strlen(what), strlen(where), 0);
      else
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,&doublevec[0],sarr,
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
        fdps.definedOn = ResultInfo::NODE; // nodes
        fdps.dofNames.push_back("-");
        fdps.unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
        Enum2String(FLUIDMECH_PRESSURE, fdps.resultName);
        numDOFs = fdps.dofNames.size();
        fdps.data.resize(numDOFs * nvx);
        fdps.entryType = ResultInfo::SCALAR;
        
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
      //     Reading turbulent kinetic energy from input file
      //-----------------------------------------------------------------------

      sprintf(what, "G/TKE_FL1");
      sprintf(where, "ZN1/VX");
      when  = timeStepNumbers_[timeStepIdx];

      if(floatDS)
        dattyp = __real_data_type__;
      else
        dattyp = __double_data_type__;

      length = 1;
      nsize  = nvx;
      iopt   = __stop_if_failed__;
        
      if(floatDS)
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,darr,sarr,
                  strlen(what), strlen(where), 0);
      else
        readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
                  &floatvec[0],iarr,carr,larr,&doublevec[0],sarr,
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
        fdps.definedOn = ResultInfo::NODE; // nodes
        fdps.dofNames.push_back("-");
        fdps.unit = MapSolTypeToUnit(FLUIDMECH_TKE);
        Enum2String(FLUIDMECH_TKE, fdps.resultName);
        numDOFs = fdps.dofNames.size();
        fdps.data.resize(numDOFs * nvx);
        fdps.entryType = ResultInfo::SCALAR;
        
        if(floatDS)
          std::copy(floatvec.begin(),
                    floatvec.begin() + (numDOFs * nvx),
                    fdps.data.begin());
        else
          std::copy(doublevec.begin(), 
                    doublevec.begin() + (numDOFs * nvx),
                    fdps.data.begin());
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

  void FileReader_CFX::GetInfosFromCommand()
  {
    std::string cmd, attrib;
    int pos=0;
    Settings& settings = Settings::Instance();
    std::ostringstream sstr;
        
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
}
