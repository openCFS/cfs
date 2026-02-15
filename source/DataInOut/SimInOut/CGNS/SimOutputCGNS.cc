// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iomanip>
#include <cmath>

#include <filesystem>
namespace fs = std::filesystem;

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"

#include "SimOutputCGNS.hh"


namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simOutputCGNS, "simOutput.CGNS")

  // ***************
  //   Constructor
  // ***************
  SimOutputCGNS::SimOutputCGNS( std::string& fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode,
                                bool isRestart )
    : SimOutput ( fileName, outputNode, infoNode, isRestart ),
      cellDim_(-1),
      outputFileOK_(false),
      writeAmplPhase_(false),
      separateFiles_(true),
      stepNumOffset_(0),
      stepValOffset_(0.0),
      writeQuadElems_(false)
  {    
    formatName_ = "cgns";
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    std::string dirString = "results_" + formatName_; 
    myParam_->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 

    myParam_->GetValue("writeQuadElems", writeQuadElems_, ParamNode::PASS );
    myParam_->GetValue("separateFiles", separateFiles_, ParamNode::PASS );

    std::string complexFormat = "realImag"; 
    myParam_->GetValue("complexFormat", complexFormat, ParamNode::PASS );
    writeAmplPhase_ = (complexFormat == "realImag") ? false : true;
  }

  // **************
  //   Destructor
  // **************
  SimOutputCGNS::~SimOutputCGNS() {
    cg_close(indexFile_[0]);

    if(separateFiles_) 
    {
      cg_close(indexFile_[1]);
    }
  }


  // *************
  //   WriteGrid
  // *************
  void SimOutputCGNS::WriteGrid() {
    if ( !outputFileOK_) {
      EXCEPTION( "File for CGNS output results is not initialized" );
    }

    if (!ptGrid_) {
      EXCEPTION("Grid pointer is not initialized" );
    }

    indexBase_.Resize(2);
    WriteNodesAndElements(0);
    WriteNodesAndElements(1);
  }

  void  SimOutputCGNS::WriteNodesAndElements(UInt baseIdx) {

    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    cellDim_ = ptGrid_->GetDim();
    cellDim_ = cellDim_ < 2 ? 2 : cellDim_;
    cellDim_ = cellDim_ > 3 ? 3 : cellDim_;
    int phys_dim = 3;

    int indexFile = indexFile_[baseIdx];

    std::string baseName = baseIdx == 0 ? "CFS_NodeSolutions" : "CFS_ElemSolutions";
    if(cg_base_write(indexFile,baseName.c_str(),cellDim_,phys_dim,&indexBase_[baseIdx]))
    {
      cg_close(indexFile);
      EXCEPTION("Cannot write base '" << baseName <<
                "':\n" << cg_get_error());
    }
    int indexBase = indexBase_[baseIdx];

    // In order to write node and element solutions afterwards, we have to
    // remember the indices of the corresponding CGNS tree nodes per region.
    UInt numRegions = ptGrid_->regionData.GetSize();
    idxZone_[baseIdx].Resize(numRegions);
    idxNodeSol_.Resize(numRegions);
    idxElemSol_.Resize(numRegions);

    // We write one zone per region, since ParaView or CFD Post do not offer a
    // way of selecting by CGNS sections.
    for(UInt i=0; i<numRegions; i++) 
    {
      RegionIdType regionId = ptGrid_->regionData[i].id;
      std::string regionName = ptGrid_->regionData[i].name;

      // Get nodes of current region.
      StdVector<UInt> nodes;
      ptGrid_->GetNodesByRegion(nodes, regionId );
      UInt numNodes = nodes.GetSize();

      // Get elements of current region.
      StdVector<Elem*> elems;
      ptGrid_->GetElems(elems, regionId);
      UInt numElems = elems.GetSize();

      #if 0
      // Try to fulfill CGNS assertion that vertices < CellDim + 1 cannot happen.
      // This may be the case when using a single LINE2 element in 2D.
      numNodes_ = ptGrid_->GetNumNodes();
      if(numNodes_ < cellDim_ + 1) 
      {
        numNodes_++;
      }
      #endif

      // Create one zone per region.
      cgsize_t isize[3][1];
      isize[0][0] = numNodes;
      isize[1][0] = numElems;
      isize[2][0] = 0;
      cg_zone_write(indexFile, indexBase, regionName.c_str(), 
                    isize[0], Unstructured, &idxZone_[baseIdx][i]);
      
      // Get node coordinates
      StdVector< Vector<double> > coords(3);
      Vector<Double> point;
      
      UInt d = 0;
      for ( ; d < 3; d++ ) {
        coords[d].Resize(numNodes);
        coords[d].Init();
      }
      UInt dim=ptGrid_->GetDim();
      for ( UInt j = 0; j < numNodes; j++ ) {
        d = 0;
        for ( ; d < dim; d++ ) {
          ptGrid_->GetNodeCoordinate(point, nodes[j]);
          coords[d][j] = point[d];
        }

        regionNodeMap_[regionId][nodes[j]] = j+1;
      }
      
      int indexCoord;
      // write grid coordinates
      cg_coord_write(indexFile,indexBase,idxZone_[baseIdx][i],
                     RealDouble,"CoordinateX",&coords[0][0],&indexCoord);
      cg_coord_write(indexFile,indexBase,idxZone_[baseIdx][i],
                     RealDouble,"CoordinateY",&coords[1][0],&indexCoord);
      cg_coord_write(indexFile,indexBase,idxZone_[baseIdx][i],
                     RealDouble,"CoordinateZ",&coords[2][0],&indexCoord);

      Elem::FEType feType = elems[0]->type;
      UInt j=1;
      
      // Check if all elems in regions are of same type.
      for( ; j<numElems; j++)
      {
        if(elems[j]->type != feType)
          break;
      }

      // If not all elements in region are of same type write a mixed section.
      if(j < numElems) {
        // Write mixed section.
        WriteMixedSection(baseIdx, elems);
      }
      else {
        // Write pure section.
        WritePureSection(baseIdx, elems);
      }


      // Write CGNS nodes for nodal and element solutions
      // First we write some information about the grid: node numbers,
      // element numbers, region ids and element types.
      int indexField;
      std::string solName;
      GridLocation_t location;

      switch(baseIdx) 
      {
      case 0:
        {
          solName = "NodeSolution";
          location = Vertex;
          
          StdVector<int> origNodeNums(numNodes);
          for(j=0; j<numNodes; j++)
          {
            origNodeNums[j] = nodes[j];
          }
          
          if(cg_sol_write(indexFile, indexBase, idxZone_[baseIdx][i],
                          solName.c_str(), location, &idxNodeSol_[i]))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write solution:\n" << cg_get_error());
          }
          
          if(cg_field_write(indexFile,indexBase,idxZone_[baseIdx][i],
                            idxNodeSol_[i],CGNSLIB_H::Integer,
                            "origNodeNums",&origNodeNums[0],&indexField))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write field:\n" << cg_get_error());
          }
        }
      break;
      
      case 1:
        {          
          StdVector<int> regionIds(numElems);
          StdVector<int> origElemNums(numElems);
          StdVector<int> elemTypes(numElems);
          for(j=0; j<numElems; j++)
          {
            Elem* ptEl = elems[j];
            regionIds[j] = ptEl->regionId;
            origElemNums[j] = ptEl->elemNum;
            elemTypes[j] = static_cast<int>(ptEl->type);
          }
          
          solName = "ElementSolution";
          location = CellCenter;
          
          if(cg_sol_write(indexFile, indexBase, idxZone_[baseIdx][i],
                          solName.c_str(), location, &idxElemSol_[i]))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write solution:\n" << cg_get_error());
          }
          
          if(cg_field_write(indexFile,indexBase,idxZone_[baseIdx][i],
                            idxElemSol_[i],CGNSLIB_H::Integer,
                            "regionId",&regionIds[0],&indexField))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write field:\n" << cg_get_error());
          }
          if(cg_field_write(indexFile,indexBase,idxZone_[baseIdx][i],
                            idxElemSol_[i],CGNSLIB_H::Integer,
                            "origElemNums",&origElemNums[0],&indexField))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write field:\n" << cg_get_error());
          }
          if(cg_field_write(indexFile,indexBase,idxZone_[baseIdx][i],
                            idxElemSol_[i],CGNSLIB_H::Integer,
                            "elemTypes",&elemTypes[0],&indexField))
          {
            cg_close(indexFile);
            EXCEPTION("Cannot write field:\n" << cg_get_error());
          }
        }
        break;
        
      default:
        break;
      }
    }
    

#if 0
    int indexFlow,indexField;

    StdVector<double> press(numNodes);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 1.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",&press[0],&indexField);

    solname = "DummyPressure1";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",&press[0],&indexField);

    for (UInt node=0; node<numNodes; node++) {
        press[node] = 2.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure2";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",&press[0],&indexField);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 3.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure3";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",&press[0],&indexField);
// create BaseIterativeData
    int nsteps=3;
    cg_biter_write(indexFile_,indexBase_,"TimeIterValues",nsteps);
// go to BaseIterativeData level and write time values
    cg_goto(indexFile_,indexBase_,"BaseIterativeData_t",1,"end");
    cgsize_t nuse=3;
    double time[3];
    time[0] = 1.0; time[1] = 2.0; time[2] = 3.0;
    cg_array_write("TimeValues",RealDouble,1,&nuse,&time);
// create ZoneIterativeData
    cg_ziter_write(indexFile_,indexBase_,indexZone_,"ZoneIterativeData");
// go to ZoneIterativeData level and give info telling which
// flow solution corresponds with which time (solname(1) corresponds
// with time(1), solname(2) with time(2), and solname(3) with time(3))
    cg_goto(indexFile_,indexBase_,"Zone_t",indexZone_,"ZoneIterativeData_t",1,"end");
    cgsize_t idata[2];
    idata[0]=32;
    idata[1]=3;
    solname="DummyPressure1      ";
    solname+="DummyPressure2      ";
    solname+="DummyPressure3      ";
    char tmp3[3][33];
    char tmp[97];
    strcpy(tmp3[0],"DummyPressure1");
    strcpy(tmp3[1],"DummyPressure2");
    strcpy(tmp3[2],"DummyPressure3");
    sprintf(tmp, "%-32s%-32s%-32s",tmp3[0],tmp3[1],tmp3[2]);
    std::cout << "tmp = " << tmp << std::endl;
    cg_array_write("FlowSolutionPointers",Character,2,idata,tmp);
// add SimulationType
    cg_simulation_type_write(indexFile_,indexBase_,TimeAccurate);
#endif
  }

  void SimOutputCGNS::WriteMixedSection(UInt baseIdx, const StdVector<Elem*>& elems) 
  {
    RegionIdType regionId = elems[0]->regionId;
    UInt numElems = elems.GetSize();
    Elem::FEType feType;
    std::string regionName = ptGrid_->regionData[regionId].name;
    
    StdVector<cgsize_t> elemData(numElems*30);

    for (UInt iElem=0, idx=0; iElem<numElems; iElem++) {
      Elem* ptEl = elems[iElem];
      feType = ptEl->type;
      UInt numElemNodes = Elem::shapes[feType].numVertices;
      if(writeQuadElems_) 
      {
        numElemNodes = Elem::shapes[feType].numNodes;
      }
      
      elemData[idx] = elemTypeMap_[feType];
      
      TranslateConnectivity(feType, &elemData[idx+1], ptEl);

      idx += numElemNodes+1;
    }

    int indexSection;
    cgsize_t nelemStart = 1, nelemEnd = numElems;
    int nbdyElem=0;
    if(cg_section_write(indexFile_[baseIdx], indexBase_[baseIdx], idxZone_[baseIdx][regionId],
                        regionName.c_str(), MIXED, nelemStart,
                        nelemEnd, nbdyElem, &elemData[0], &indexSection))
    {
      cg_close(indexFile_[baseIdx]);
      EXCEPTION("Cannot write mixed section for '" << regionName <<
                "':\n" << cg_get_error());
    }
  }
  
  void SimOutputCGNS::WritePureSection(UInt baseIdx, const StdVector<Elem*>& elems) 
  {
    RegionIdType regionId = elems[0]->regionId;
    UInt numElems = elems.GetSize();
    Elem::FEType feType = elems[0]->type;
    ElementType_t eType = elemTypeMap_[feType];
    std::string regionName = ptGrid_->regionData[regionId].name;
    UInt numElemNodes = Elem::shapes[feType].numVertices;
    if(writeQuadElems_) 
    {
      numElemNodes = Elem::shapes[feType].numNodes;
    }
    
    StdVector<cgsize_t> elemData(numElems*numElemNodes);

    for (UInt iElem=0; iElem<numElems; iElem++) {
      Elem* ptEl = elems[iElem];
      TranslateConnectivity(feType, &elemData[iElem*numElemNodes], ptEl);
    }

    int indexSection;
    cgsize_t nelemStart = 1, nelemEnd = numElems;
    int nbdyElem=0;
    if(cg_section_write(indexFile_[baseIdx], indexBase_[baseIdx], idxZone_[baseIdx][regionId],
                        regionName.c_str(), eType, nelemStart,
                        nelemEnd, nbdyElem, &elemData[0], &indexSection))
    {
      cg_close(indexFile_[baseIdx]);
      EXCEPTION("Cannot write pure section for '" << regionName <<
                "':\n" << cg_get_error());
    }
  }

  void SimOutputCGNS::InitElemTypeMap(){
    elemTypeMap_.clear();
    elemTypeMap_[Elem::ET_UNDEF]   = CGNSLIB_H::ElementTypeNull;
    elemTypeMap_[Elem::ET_POINT]   = CGNSLIB_H::NODE;
    elemTypeMap_[Elem::ET_LINE2]   = CGNSLIB_H::BAR_2;
    elemTypeMap_[Elem::ET_TRIA3]   = CGNSLIB_H::TRI_3;
    elemTypeMap_[Elem::ET_QUAD4]   = CGNSLIB_H::QUAD_4;
    elemTypeMap_[Elem::ET_TET4]    = CGNSLIB_H::TETRA_4;
    elemTypeMap_[Elem::ET_PYRA5]   = CGNSLIB_H::PYRA_5;
    elemTypeMap_[Elem::ET_WEDGE6]  = CGNSLIB_H::PENTA_6;
    elemTypeMap_[Elem::ET_HEXA8]   = CGNSLIB_H::HEXA_8;

    if(writeQuadElems_) 
    {
      elemTypeMap_[Elem::ET_LINE3]   = CGNSLIB_H::BAR_3;
      elemTypeMap_[Elem::ET_TRIA6]   = CGNSLIB_H::TRI_6;
      elemTypeMap_[Elem::ET_QUAD8]   = CGNSLIB_H::QUAD_8;
      elemTypeMap_[Elem::ET_QUAD9]   = CGNSLIB_H::QUAD_9;
      elemTypeMap_[Elem::ET_TET10]   = CGNSLIB_H::TETRA_10;
#if CGNS_VERSION > 2550
      elemTypeMap_[Elem::ET_PYRA13]  = CGNSLIB_H::PYRA_13;
#endif
      elemTypeMap_[Elem::ET_PYRA14]  = CGNSLIB_H::PYRA_14;
      elemTypeMap_[Elem::ET_WEDGE15] = CGNSLIB_H::PENTA_15;
      elemTypeMap_[Elem::ET_WEDGE18] = CGNSLIB_H::PENTA_18;
      elemTypeMap_[Elem::ET_HEXA20]  = CGNSLIB_H::HEXA_20;
      elemTypeMap_[Elem::ET_HEXA27]  = CGNSLIB_H::HEXA_27;
    }
    else 
    {
      WARN("Quadratic element types will be reduced to linear ones for CGNS!");

      elemTypeMap_[Elem::ET_LINE3]   = CGNSLIB_H::BAR_2;
      elemTypeMap_[Elem::ET_TRIA6]   = CGNSLIB_H::TRI_3;
      elemTypeMap_[Elem::ET_QUAD8]   = CGNSLIB_H::QUAD_4;
      elemTypeMap_[Elem::ET_QUAD9]   = CGNSLIB_H::QUAD_4;
      elemTypeMap_[Elem::ET_TET10]   = CGNSLIB_H::TETRA_4;
      elemTypeMap_[Elem::ET_PYRA13]  = CGNSLIB_H::PYRA_5;
      elemTypeMap_[Elem::ET_PYRA14]  = CGNSLIB_H::PYRA_5;
      elemTypeMap_[Elem::ET_WEDGE15] = CGNSLIB_H::PENTA_6;
      elemTypeMap_[Elem::ET_WEDGE18] = CGNSLIB_H::PENTA_6;
      elemTypeMap_[Elem::ET_HEXA20]  = CGNSLIB_H::HEXA_8;
      elemTypeMap_[Elem::ET_HEXA27]  = CGNSLIB_H::HEXA_8;
    }
  }

  void SimOutputCGNS::TranslateConnectivity(Elem::FEType feType,
                                           cgsize_t* cgnsConn,
                                           Elem* elem)
  {
    StdVector<UInt>& connect = elem->connect;
    RegionIdType regionId = elem->regionId;
    UInt numElemNodes = Elem::shapes[feType].numVertices;
    if(writeQuadElems_) 
    {
      numElemNodes = Elem::shapes[feType].numNodes;
    }

    static const int trDefault[27] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
      20, 21, 22, 23, 24, 25, 26
    };

    // Map from a CGNS HEXA_20 connectivity
    static const int trHEX20[20] = {
      0, 1, 2, 3,
      4, 5, 6, 7,
      8, 9, 10, 11,
      16, 17, 18, 19,
      12, 13, 14, 15
    };
    // Map from a CGNS HEXA_27 connectivity
    static const int trHEX27[27] = {
      0, 1, 2, 3,
      4, 5, 6, 7,
      8, 9, 10, 11,
      16, 17, 18, 19,
      12, 13, 14, 15,
      21, 22, 23, 24,
      20, 25,
      26
    };
    // Map from a CGNS PENTA_15 connectivity
    static const int trPRI15[15] = {
      0, 1, 2, 3, 4, 5, 6, 7, 8,
      12, 13, 14,
      9, 10, 11 
    };
    // Map from a CGNS PENTA_18 connectivity
    static const int trPRI18[18] = {
      0, 1, 2, 3, 4, 5,
      6, 7, 8, 12, 13, 14,
      9, 10, 11, 15, 16, 17
    };
    
    const int *tr;
    switch(feType) {
    case Elem::ET_HEXA20:
      tr = trHEX20;
      break;
    case Elem::ET_HEXA27:
      tr = trHEX27;
      break;
    case Elem::ET_WEDGE15:
      tr = trPRI15;
      break;
    case Elem::ET_WEDGE18:
      tr = trPRI18;
      break;
    default:
      tr = trDefault;
      break;
    }

    for(UInt n = 0; n<numElemNodes; n++ ) {
      cgnsConn[n] = regionNodeMap_[regionId][connect[tr[n]]];
    }
  }

  void  SimOutputCGNS::NodeElemDataTransient(const bool definedOnNode,
                                             RegionSolsType& regionSols,
                                             const UInt step, 
                                             const Double time)
  {
    RegionSolsType::const_iterator solIt, solEnd;
    UInt baseIdx = definedOnNode ? 0 : 1;

    // Write actual solution array.
    solIt = regionSols.begin();
    solEnd = regionSols.end();
    for( ; solIt != solEnd; solIt++ ) 
    {
      std::string solName = solIt->first;
      RegionSolsType::mapped_type::const_iterator it, end;

      it = solIt->second.begin();
      end = solIt->second.end();
      
      for( ; it != end; it++ ) 
      {
        RegionIdType regionId = it->first;
        int* indexSol = &idxNodeSol_[regionId];
    
        // Check if solution node is already present, and write it if necessary.
        if(!definedOnNode) 
        {
          indexSol = &idxElemSol_[regionId];
        }    

        const Vector<Double>& sol = it->second;
        int indexField;
    
        if(cg_field_write(indexFile_[baseIdx], indexBase_[baseIdx],
                          idxZone_[baseIdx][regionId],
                          *indexSol, RealDouble, solName.c_str(),
                          &sol[0], &indexField)) 
        {
          cg_close(indexFile_[baseIdx]);
          EXCEPTION("Cannot write field:\n" <<
                    cg_get_error());
        }
      }
    }
  }

  void SimOutputCGNS::Init(Grid * ptGrid, bool printGridOnly )
  {
    ptGrid_=ptGrid;

    InitElemTypeMap();
    
    // concatenate output file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    
    fs::path fileName (dirName_);
    
#if CGNS_VERSION > 2550
    // Set type of mid level library to use.
    int file_type = CG_FILE_ADF;
    std::string mll = "adf";
    if(myParam_->Has("mll")) 
    {
      // We set the default mid-level library type to ADF, to support also
      // older post-processors per default.

      mll = myParam_->Get("mll")->As<std::string>();

      if(mll == "hdf5") 
      {
        file_type = CG_FILE_HDF5;
      }
#if CGNS_VERSION > 3000
      if(mll == "adf2")
      {
        file_type = CG_FILE_ADF2;
      }
#endif
    }

    if(cg_set_file_type(file_type))
    {
      EXCEPTION("Cannot set file type to '" << mll <<
                "':\n" << cg_get_error());
    }
#endif

    indexFile_.Resize(2);
    if(!separateFiles_) 
    {
      fileName = fileName / (fileName_ + std::string(".cgns"));

      // Make sure a previous file with the same name gets removed.
      try {
        if(fs::exists(fileName)) { fs::remove( fileName ); }
      } catch (std::exception &ex) {
        EXCEPTION(ex.what());
      }

      // Open CGNS file for writing
      if (cg_open(fileName.string().c_str(),CG_MODE_WRITE,&indexFile_[0]))
      {
        EXCEPTION("Cannot open '" << fileName.string() <<
                  "':\n" << cg_get_error());
      }

      indexFile_[1] = indexFile_[0];
    }
    else 
    {
      // Open CGNS file for writing nodal results
      fs::path fnNode = fileName / (fileName_ + std::string("_NodeSol.cgns"));

      // Make sure a previous file with the same name gets removed.
      try {
        if(fs::exists(fnNode)) { fs::remove( fnNode ); }
      } catch (std::exception &ex) {
        EXCEPTION(ex.what());
      }

      if (cg_open(fnNode.string().c_str(),CG_MODE_WRITE,&indexFile_[0]))
      {
        EXCEPTION("Cannot open '" << fnNode.string() <<
                  "':\n" << cg_get_error());
      }

      // Open CGNS file for writing element results
      fs::path fnElem = fileName / (fileName_ + std::string("_ElemSol.cgns"));

      // Make sure a previous file with the same name gets removed.
      try {
        if(fs::exists(fnElem)) { fs::remove( fnElem ); }
      } catch (std::exception &ex) {
        EXCEPTION(ex.what());
      }

      if (cg_open(fnElem.string().c_str(),CG_MODE_WRITE,&indexFile_[1]))
      {
        EXCEPTION("Cannot open '" << fnElem.string() <<
                  "':\n" << cg_get_error());
      }
    }

    // Everything fine!
    outputFileOK_ = true;

    WriteGrid();
  }

  void SimOutputCGNS::BeginMultiSequenceStep( UInt step,
                                              BasePDE::AnalysisType type,
                                              UInt numSteps )
  {
    if(step != 1) 
    {
      EXCEPTION("Only one sequence step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }
    
    if(numSteps != 1) 
    {
      EXCEPTION("Only one time/freq step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }
  }
  
  void SimOutputCGNS::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd,
                                      bool isHistory )
  {
  }

  void SimOutputCGNS::BeginStep( UInt stepNum, Double stepVal )
  {
    if(stepNum != 1) 
    {
      EXCEPTION("Only one time/freq step supported!\n" <<
                "Note, that this is no limitation of CGNS but " <<
                "of the current implementation.");
    }
    
    actStep_ = stepNum + stepNumOffset_;
    actStepVal_ = stepVal + stepValOffset_;
    resultMap_.clear();
  }
  
  void SimOutputCGNS::AddResult( shared_ptr<BaseResult>  sol )
  {
    std::string resultName = sol->GetResultInfo()->resultName;
    resultMap_[resultName].Push_back( sol );
  }
  
  void SimOutputCGNS::FinishStep( )
  {
    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // check if result is defined on nodes or elements
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;
      if(entityType == ResultInfo::NODE ||
         entityType == ResultInfo::ELEMENT ||
         entityType == ResultInfo::SURF_ELEM ) {

        if(actInfo.resultType == NO_SOLUTION_TYPE) 
        {
          continue;
        }

        std::string title =  SolutionTypeEnum.ToString( actInfo.resultType );
        
        // if result could not be mapped, omit it
        if( title == "") {
          std::stringstream warning;
          warning  <<  "Result '" << actInfo.resultName
                   << "' could not be mappted to CGNS result type and is "
                   << "omitted!";
          WARN( warning.str().c_str() );
          continue;
        }
        
       
        // Determine type of entity the result is defined on
        bool definedOnNode = (entityType == ResultInfo::NODE);
        
        // Create map of arrays, which are then written to CGNS.
        RegionSolsType regionSols;
        FillRegionSols(regionSols, actResults, entityType );
        NodeElemDataTransient( definedOnNode, regionSols, actStep_, actStepVal_);
      } 
      else 
      {
        WARN( "CGNS can only write results on element and nodes." );
      }        
      
    } // over all result types
  }
    
  
  void SimOutputCGNS::FinishMultiSequenceStep() {
    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }

  void SimOutputCGNS::Finalize() {
    // return, if no commandLine handler or
    // global root ParaemNode are present
    if( !progOpts || !myParam_ )
      return;

    std::string paramFile = progOpts->GetParamFileStr();
    std::string matFile = myParam_->GetRoot()->Get("fileFormats")
                          ->Get("materialData")
                          ->Get("file")->As<std::string>();
    WriteXmlFiles(paramFile, matFile);
  }
  
  void SimOutputCGNS::FillRegionSols(RegionSolsType& regionSols, 
                                     const StdVector<shared_ptr<BaseResult> > & solList,
                                     ResultInfo::EntityUnknownType entityType ) {

    static const double H180DEG_OVER_PI = 180.0 / M_PI;

    BaseMatrix::EntryType entryType = solList[0]->GetEntryType();

    ResultInfo & actResultInfo = *(solList[0]->GetResultInfo());
    std::string resultName = actResultInfo.resultName;
    StdVector<std::string> dofNames = actResultInfo.dofNames;
    UInt numDofs = actResultInfo.dofNames.GetSize();
    UInt idx;

    LOG_DBG(simOutputCGNS) << "Filling vectors for result '" 
                           << resultName << "' on ";

    std::map< RegionIdType, std::map<UInt, UInt> >::const_iterator rNMIt, rNMEnd;
    rNMIt = regionNodeMap_.begin();
    rNMEnd = regionNodeMap_.end();
    
    for( ; rNMIt != rNMEnd; rNMIt++ ) {
      RegionIdType regionId = rNMIt->first;
      shared_ptr<EntityList> entityList;
      bool solDefined=false;
      shared_ptr<BaseResult> result;

      UInt i = 0;
      for( ; i < solList.GetSize(); i++ ) {
        std::string entityListName = solList[i]->GetEntityList()->GetName();
        if(regionId == ptGrid_->GetRegion().Parse(entityListName))
        {
          LOG_DBG(simOutputCGNS) << entityListName;
          entityList = solList[i]->GetEntityList();
          result = solList[i];
          solDefined = true;
          
          break;
        }
      }

      if( entityType == ResultInfo::ELEMENT ||
          entityType == ResultInfo::SURF_ELEM ) {
        
        if(!solDefined) 
        {
          ElemList* el = new ElemList(ptGrid_);
          el->SetRegion(regionId);
          entityList.reset(el);
        }
        
        UInt numEntities = entityList->GetSize();
        EntityIterator it = entityList->GetIterator();
        
        // === Element Results ===
        for(UInt iDof=0; iDof<numDofs; iDof++) 
        {
          std::string resultDofName = (numDofs == 1 ? resultName : resultName + "_" + dofNames[iDof]);

          if( entryType == BaseMatrix::DOUBLE ) {
            // real valued results
            regionSols[resultDofName][regionId].Resize( numEntities );
            regionSols[resultDofName][regionId].Init();

            Vector<Double> dummyVec;
            if(!solDefined) { dummyVec.Resize(numEntities*numDofs); dummyVec.Init(); }
            
            Vector<Double> & actSol = solDefined ? dynamic_cast<Result<Double>&>
                                      (*result).GetVector() : dummyVec;
            
            for( idx = 0, it.Begin(); !it.IsEnd(); it++, idx++ ) {
              regionSols[resultDofName][regionId][idx] = actSol[it.GetPos()*numDofs+iDof];
            }
          }
          else
          {
            // complex valued results
            std::string resNameRe = resultDofName + "_Re";
            std::string resNameIm = resultDofName + "_Im";
            std::string resNameAmpl = resultDofName + "_Ampl";
            std::string resNamePhase = resultDofName + "_Phase";
            
            if(writeAmplPhase_) 
            {
              regionSols[resNameAmpl][regionId].Resize( numEntities );
              regionSols[resNameAmpl][regionId].Init();
              regionSols[resNamePhase][regionId].Resize( numEntities );
              regionSols[resNamePhase][regionId].Init();
            }
            else
            {
              regionSols[resNameRe][regionId].Resize( numEntities );
              regionSols[resNameRe][regionId].Init();
              regionSols[resNameIm][regionId].Resize( numEntities );
              regionSols[resNameIm][regionId].Init();
            }
            
            
            Vector<Complex> dummyVec;
            if(!solDefined) { dummyVec.Resize(numEntities*numDofs); dummyVec.Init(); }
            
            Vector<Complex> & actSol = solDefined ? dynamic_cast<Result<Complex>&>
                                       (*result).GetVector() : dummyVec;
            for( idx=0, it.Begin(); !it.IsEnd(); it++, idx++ ) {
              Complex sol = actSol[it.GetPos()*numDofs+iDof];
              
              if(writeAmplPhase_) 
              {                
                regionSols[resNameAmpl][regionId][idx] = hypot(sol.real(), sol.imag());
                regionSols[resNamePhase][regionId][idx] = 
                  (std::abs(sol.imag()) > 1e-16) ?
                  std::atan2( sol.imag(), sol.real() ) * H180DEG_OVER_PI : 
                  ( sol.real() < 0.0 ) ? 180 : 0 ;
              }
              else
              {
                regionSols[resNameRe][regionId][idx] = sol.real();
                regionSols[resNameIm][regionId][idx] = sol.imag();
              }
            }
          }
        }
      } else if ( entityType == ResultInfo::NODE ) {
        std::map<UInt, UInt>& nodeMap = regionNodeMap_[regionId];
        
        if(!solDefined) 
        {
          NodeList* nl = new NodeList(ptGrid_);
          nl->SetNodesOfRegion(regionId);
          entityList.reset(nl);
        }
        
        UInt numEntities = entityList->GetSize();
        EntityIterator it = entityList->GetIterator();

        // === Nodal Results ===
        for(UInt iDof=0; iDof<numDofs; iDof++) 
        {
          std::string resultDofName = (numDofs == 1 ? resultName : resultName + "_" + dofNames[iDof]);
          
          if( entryType == BaseMatrix::DOUBLE ) {
            // real valued results
            regionSols[resultDofName][regionId].Resize( numEntities );
            regionSols[resultDofName][regionId].Init();

            Vector<Double> dummyVec;
            if(!solDefined) { dummyVec.Resize(numEntities*numDofs); dummyVec.Init(); }            

            Vector<Double> & actSol = solDefined ? dynamic_cast<Result<Double>&>
                                      (*result).GetVector() : dummyVec;
            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt nodeNum = it.GetNode();
              if(nodeMap.find(nodeNum) != nodeMap.end()) 
              {
                idx = nodeMap[nodeNum]-1;
                regionSols[resultDofName][regionId][idx] = actSol[it.GetPos()*numDofs+iDof];
              }
            }
          }
          else
          {
            // complex valued results
            std::string resNameRe = resultDofName + "_Re";
            std::string resNameIm = resultDofName + "_Im";
            std::string resNameAmpl = resultDofName + "_Ampl";
            std::string resNamePhase = resultDofName + "_Phase";
            
            if(writeAmplPhase_) 
            {
              regionSols[resNameAmpl][regionId].Resize( numEntities );
              regionSols[resNameAmpl][regionId].Init();
              regionSols[resNamePhase][regionId].Resize( numEntities );
              regionSols[resNamePhase][regionId].Init();
            }
            else
            {
              regionSols[resNameRe][regionId].Resize( numEntities );
              regionSols[resNameRe][regionId].Init();
              regionSols[resNameIm][regionId].Resize( numEntities );
              regionSols[resNameIm][regionId].Init();
            }
            
            Vector<Complex> dummyVec;
            if(!solDefined) { dummyVec.Resize(numEntities*numDofs); dummyVec.Init(); }
            
            Vector<Complex> & actSol = solDefined ? dynamic_cast<Result<Complex>&>
                                       (*result).GetVector() : dummyVec;

            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt nodeNum = it.GetNode();
              if(nodeMap.find(nodeNum) != nodeMap.end()) 
              {
                idx = nodeMap[nodeNum]-1;
                Complex sol = actSol[it.GetPos()*numDofs+iDof];
                
                if(writeAmplPhase_) 
                {
                  regionSols[resNameAmpl][regionId][idx] = hypot(sol.real(), sol.imag());
                  regionSols[resNamePhase][regionId][idx] = 
                    (std::abs(sol.imag()) > 1e-16) ?
                    std::atan2( sol.imag(), sol.real() ) * H180DEG_OVER_PI : 
                    ( sol.real() < 0.0 ) ? 180 : 0 ;
                }
                else
                {
                  regionSols[resNameRe][regionId][idx] = sol.real();
                  regionSols[resNameIm][regionId][idx] = sol.imag();
                }
              }
            }
          }
        }      
      } else {
        EXCEPTION( "Can only map nodal / element results to grid" );
      }
    }
  }

  void SimOutputCGNS::WriteXmlFiles( fs::path simFile, fs::path matFile ) {
    if( isRestart_)
          return;

    for(UInt idx = 0; idx < 2; idx++) 
    {
      cg_goto(indexFile_[idx], indexBase_[idx], "end");
      
      std::ifstream fin;
      std::ostringstream dumpStr;
      
      // open external Files
      StdVector<fs::path> filePaths(2);
      StdVector<std::string> setNames(2);
      filePaths[0] = simFile;
      setNames[0] = "ParameterFile";
      filePaths[1] = matFile;
      setNames[1] = "MaterialFile";
      
      for(UInt i=0; i<filePaths.GetSize(); i++)
      {
        fin.open( filePaths[i], std::ios::binary );
        
        if(fin.fail())
          EXCEPTION("Cannot open file '" << filePaths[i]
                    <<"' to dump into HDF5!");
        
        // seek to the end of the file
        fin.seekg (0, std::ios::end);
        UInt numBytes = fin.tellg();
        fin.seekg (0, std::ios::beg);
        
        std::string str;
        str.resize(numBytes);
        fin.read(&str[0], numBytes);
        fin.close();
        
        cg_descriptor_write(setNames[i].c_str(), str.c_str());
      }
    }
  }
} // namespace CoupledField
