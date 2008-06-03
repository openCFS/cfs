#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <stdio.h>
#include <iomanip>

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
    basename_ = "./";
    basename_+= name_;
    basename_+= "_data/";
  }

  FileReader_CFX::~FileReader_CFX()
  {
  }

  void FileReader_CFX::Init()
  {
    Settings& settings = Settings::Instance();

    if(settings.GetDouble("timeStep") < 0)
      EXCEPTION("No proper time step has been specified! Use --timestep X.");
        
    //-----------------------------------------------------------------------
    //     Open RESULTS file
    //-----------------------------------------------------------------------

    snprintf(fn, sizeof(fn), "%s%s.res", basename_.c_str(), name_.c_str());
    whatfile = __io_open_primaryfile__;
    if(settings.GetInt("verbose"))
    {
      std::cout << "Trying to open CFX results file " << fn
                << "." << std::endl; 
    }
    
    openfile_(&nerr, fn, &whatfile, strlen(fn)); 
    CHECK_CFX_IO(nerr);

    intvec.resize(BIGMEM*3);
    floatvec.resize(BIGMEM);

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
    
    numFiles_ = 0;
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
        
      snprintf(fn, sizeof(fn), "%s%s/%s", basename_.c_str(), name_.c_str(), trnnam);
            
      infile.clear();
      infile.open(fn);
      if (infile)
      {
        infile.close();

        numFiles_++;
        transientFNs_.push_back(fn);
        timeStepNumbers_.push_back(its);
      }
            
    }

    numFiles_ = numFiles_ > (UInt) settings.GetInt("numSteps") ? settings.GetInt("numSteps") : numFiles_;
        

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

        
    //-----------------------------------------------------------------------
    //     Open DEFINITION file
    //-----------------------------------------------------------------------

    std::vector< std::string > defFileNames;
    if(settings.GetString("defFile") != "")
      defFileNames.push_back(basename_ + settings.GetString("defFile").c_str());
    if(defFile != "")
      defFileNames.push_back(basename_ + defFile);
    defFileNames.push_back(basename_ + name_ + ".def");
    defFile = "";

    for(UInt i=0; i<defFileNames.size(); i++)
    {
      if(settings.GetInt("verbose"))
      {    
        std::cerr << "Trying to open deffile: " << defFileNames[i] << " ";
      }

      infile.clear();
      infile.open(defFileNames[i].c_str());
      if (infile)
      {
        defFile = defFileNames[i];
        infile.close();
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
      std::cerr << "Can not find definition file." << std::endl;
      exit(1);
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
        
    numPartitions_ = nes;
    elsize_.resize(numPartitions_);
    MpCCInodes_.resize(numPartitions_);
    MpCCIelems_.resize(numPartitions_);
        
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
      int nvxel = 0;
        
      //
      //---- number of nodes per element
      //
      if (ielem[ies-1] == 4) nvxel = 4;
      if (ielem[ies-1] == 5) nvxel = 6;
      if (ielem[ies-1] == 6) nvxel = 8;
      if (ielem[ies-1] == 7) nvxel = 5;
          
      elsize_[ies-1] = nvxel;
            
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


      MpCCIelems_[ies-1] = nsize;
      MpCCInodes_[ies-1] = nvx;

      if(settings.GetInt("verbose"))
      {
        printf("ES: %d\nNumber of Elements: %d\n", ies, nsize);

        std::cout << "Partition " << (ies)
                  << " nodes: " << MpCCInodes_[ies-1]
                  << " elems: " << MpCCIelems_[ies-1]
                  << " elsize: " << elsize_[ies-1]
                  << std::endl;
      }
    }

    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);
  }
    
  void FileReader_CFX::ReadNodalCoords(std::vector<Double> & NODECOORD,
                                       const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    NODECOORD.resize(MpCCInodes_[partitionIdx]*3);
        
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
    nsize  = MpCCInodes_[partitionIdx];
    iopt = __stop_if_failed__;
    doublevec.resize(MpCCInodes_[partitionIdx]*3);
    
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
                                    std::vector<UInt> & numNodesPerElem,
                                    std::vector<UInt> & elemTypes,
                                    const UInt partitionIdx)
  {
    Settings& settings = Settings::Instance();
    int numNodes = elsize_[partitionIdx];
    int numElems = MpCCIelems_[partitionIdx];
    int elemType = ET_UNDEF;

    snprintf(fn, sizeof(fn),"%s", defFile.c_str());
    whatfile = __io_open_primaryfile__;
    openfile_(&nerr, fn, &whatfile, strlen(fn)); 
    CHECK_CFX_IO(nerr);

    switch(numNodes)
    {
    case 4:
      elemType = ET_TET4;
      break;
        
    case 6:
      elemType = ET_WEDGE6;
      break;

    case 8:
      elemType = ET_HEXA8;
      break;

    case 5:
      elemType = ET_PYRA5;
      break;
    }
        
    //---- reading connectivity for each element set
    //     outer loop:  i_element
    //     inner loop:  i_vx per element (1 ... nelvx)
        
    sprintf(what,"G/KVXPE");
        
    //
    //---- where = ZN1/ESn where n is integer from 1 to nes
    //
        
    sprintf(where, "ZN1/ES%d", partitionIdx+1);
    when  = 0;
    
    dattyp = __int_data_type__;
    length = numNodes;
    nsize  = BIGMEM;
    iopt   = __stop_if_failed__;
        
        
    readlong_(&dattyp,&nerr,what,where,&when,&nsize,&iopt,
              rarr,&intvec[0],carr,larr,darr,sarr,
              strlen(what), strlen(where), 0);
    CHECK_CFX_IO(nerr);

    std::cout  << "read connectivity from def file" << std::endl;

    whatfile = __io_close_primaryfile__;
    closefile_(&nerr, &whatfile);
    CHECK_CFX_IO(nerr);

    if(settings.GetInt("verbose"))
    {
      printf("Length of Connectivity array: %d\n", nsize);
    }

    TOPOLOGYDATA.resize(numElems*numNodes);
    numNodesPerElem.resize(numElems);
    elemTypes.resize(numElems);

    int baseIdx;
        
    for(int i=0; i<numElems; i++) 
    {
      numNodesPerElem[i] = numNodes;
      elemTypes[i] = elemType;


      if(elemType == ET_HEXA8)
      {
        baseIdx = i*8;
        // Reihenfolge für MpCCI
        TOPOLOGYDATA[baseIdx + 0] = intvec[baseIdx + 4];
        TOPOLOGYDATA[baseIdx + 1] = intvec[baseIdx + 6];
        TOPOLOGYDATA[baseIdx + 2] = intvec[baseIdx + 7];
        TOPOLOGYDATA[baseIdx + 3] = intvec[baseIdx + 5];
        TOPOLOGYDATA[baseIdx + 4] = intvec[baseIdx + 0];
        TOPOLOGYDATA[baseIdx + 5] = intvec[baseIdx + 2];
        TOPOLOGYDATA[baseIdx + 6] = intvec[baseIdx + 3];
        TOPOLOGYDATA[baseIdx + 7] = intvec[baseIdx + 1];

        // Reihenfolge für HDF5
        /*
          TOPOLOGYDATA[baseIdx + 0] = intvec[baseIdx + 0]; // 1
          TOPOLOGYDATA[baseIdx + 1] = intvec[baseIdx + 2]; // 22
          TOPOLOGYDATA[baseIdx + 2] = intvec[baseIdx + 3]; // 23
          TOPOLOGYDATA[baseIdx + 3] = intvec[baseIdx + 1]; // 2
          TOPOLOGYDATA[baseIdx + 4] = intvec[baseIdx + 4]; // 1072
          TOPOLOGYDATA[baseIdx + 5] = intvec[baseIdx + 6]; // 1093
          TOPOLOGYDATA[baseIdx + 6] = intvec[baseIdx + 7]; // 1094
          TOPOLOGYDATA[baseIdx + 7] = intvec[baseIdx + 5]; // 1073
        */
      }
    }
  }
    
  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CFX::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
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

    for(UInt actPart=0; actPart < numPartitions_; actPart++)
    {
      int nvx = MpCCInodes_[actPart];
      FlowDataType& fd = nodalFlowData[actPart];
      UInt numDOFs;
      
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
        fdps.unit = "m s^-1";
        Enum2String(FLUIDMECH_VELOCITY, fdps.resultName);
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
        std::cerr << "WARNING: CFX dataset does not contain pressure!"
                  << std::endl;                           
      }
      else 
      {
        FlowDataPartStruct& fdps = fd[FLUIDMECH_PRESSURE];
        fdps.isActive = true; // all partitions have results
        fdps.definedOn = ResultInfo::NODE; // nodes
        fdps.dofNames.push_back("-");
        fdps.unit = "Pa";
        fdps.resultName = "fluidMechPressure";
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
      std::cerr << "Unrecognized time unit found: "
                << solTimeUnit << std::endl;
      exit(1);
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
    userData["CFX_COMMANDS"] = userDataCFX_COMMANDS;
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
