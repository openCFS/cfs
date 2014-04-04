// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iomanip>
#include <cmath>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/Logging/LogConfigurator.hh"

#include "SimOutputCGNS.hh"

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(simOutputCGNS)
  DEFINE_LOG(simOutputCGNS, "simOutput.CGNS")

  // ***************
  //   Constructor
  // ***************
  SimOutputCGNS::SimOutputCGNS( std::string& fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode,
                                bool isRestart )
    : SimOutput ( fileName, outputNode, infoNode, isRestart ),
      indexFile_(-1),
      indexBase_(-1),
      indexZone_(-1),
      indexNodeSol_(-1),
      indexElemSol_(-1),
      cellDim_(-1),
      numNodes_(0),
      outputFileOK_(false),
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
  }


  // **************
  //   Destructor
  // **************
  SimOutputCGNS::~SimOutputCGNS() {
    if(indexFile_ != -1) 
    {
      cg_close(indexFile_);
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

    cellDim_ = ptGrid_->GetDim();
    cellDim_ = cellDim_ < 2 ? 2 : cellDim_;
    cellDim_ = cellDim_ > 3 ? 3 : cellDim_;
    int phys_dim = 3;

    strcpy(baseName_,"CFS_Simulation");
    cg_base_write(indexFile_,baseName_,cellDim_,phys_dim,&indexBase_);

    WriteNodesAndElements();
  }

  void  SimOutputCGNS::WriteNodesAndElements() {

    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    // Try to fulfill CGNS assertion that vertices < CellDim + 1 cannot happen.
    // This may be the case when using a single LINE2 element in 2D.
    numNodes_ = ptGrid_->GetNumNodes();
    if(numNodes_ < cellDim_ + 1) 
    {
      numNodes_++;
    }
    
    UInt numElems = 0;
    for(UInt i=0, n=ptGrid_->regionData.GetSize(); i<n; i++) 
    {
      RegionIdType regionId = ptGrid_->regionData[i].id;

      StdVector<Elem*> elems;
      
      ptGrid_->GetElems(elems, regionId);
      UInt nElems = elems.GetSize();

      numElems += nElems;
    }    

    // create zone
    strcpy(zoneName_,"CFS_Mesh");
    cgsize_t isize[3][1];
    isize[0][0] = numNodes_;
    isize[1][0] = numElems;
    isize[2][0] = 0;
    cg_zone_write(indexFile_,indexBase_,zoneName_,isize[0],Unstructured,&indexZone_);

    // coordinates
    StdVector< Vector<double> > coords(3);
    Vector<Double> point;

    UInt d = 0;
    for ( ; d < 3; d++ ) {
      coords[d].Resize(numNodes_);
      coords[d].Init();
    }
    UInt dims=ptGrid_->GetDim();
    d = 0;
    for ( ; d < dims; d++ ) {
      for ( UInt i = 0, n=ptGrid_->GetNumNodes(); i < n; i++ ) {
        ptGrid_->GetNodeCoordinate(point,i+1);
        coords[d][i] = point[d];
      }
    }
    
    int indexCoord;
    // write grid coordinates
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateX",&coords[0][0],&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateY",&coords[1][0],&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateZ",&coords[2][0],&indexCoord);

    UInt elemRangeStart = 1;
    StdVector<int> regionIds;
    StdVector<int> origElemNums;
    StdVector<int> elemTypes;
    
    for(UInt i=0, n=ptGrid_->regionData.GetSize(); i<n; i++) 
    {
      RegionIdType regionId = ptGrid_->regionData[i].id;

      StdVector<Elem*> elems;
      
      ptGrid_->GetElems(elems, regionId);
      UInt nElems = elems.GetSize();
      Elem::FEType feType = elems[0]->type;
      UInt j=1;
      regionIds.Resize(regionIds.GetSize()+nElems);
      origElemNums.Resize(regionIds.GetSize()+nElems);
      elemTypes.Resize(regionIds.GetSize()+nElems);
      
      // Check if all elems in regions are of same type
      for( ; j<nElems; j++)
      {
        if(elems[j]->type != feType)
          break;
      }

      if(j < nElems) 
      {
        // Write mixed section.
        WriteMixedSection(elems,
                          ptGrid_->regionData[regionId].name,
                          regionIds,
                          origElemNums,
                          elemTypes,
                          elemRangeStart);
      }
      else
      {
        // Write pure section.
        WritePureSection(elems,
                         ptGrid_->regionData[regionId].name,
                         regionIds,
                         origElemNums,
                         elemTypes,
                         elemRangeStart);
      }
    }

    int indexField;
    std::string solName = "NodeSolution";
    int* indexSol = &indexNodeSol_;
    GridLocation_t location = Vertex;

    if(cg_sol_write(indexFile_, indexBase_, indexZone_,
                    solName.c_str(), location, indexSol))
    {
      cg_close(indexFile_);
      EXCEPTION("Cannot write solution:\n" << cg_get_error());
    }
    
    solName = "ElementSolution";
    indexSol = &indexElemSol_;
    location = CellCenter;

    if(cg_sol_write(indexFile_, indexBase_, indexZone_,
                    solName.c_str(), location, indexSol))
    {
      cg_close(indexFile_);
      EXCEPTION("Cannot write solution:\n" << cg_get_error());
    }

    cg_field_write(indexFile_,indexBase_,indexZone_,indexElemSol_,CGNSLIB_H::Integer,
            "RegionId",&regionIds[0],&indexField);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexElemSol_,CGNSLIB_H::Integer,
            "OrigElemNum",&origElemNums[0],&indexField);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexElemSol_,CGNSLIB_H::Integer,
            "ElemType",&elemTypes[0],&indexField);

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

  void SimOutputCGNS::WriteMixedSection(const StdVector<Elem*>& elems,
                                        const std::string& name,
                                        StdVector<int>& regionIds,
                                        StdVector<int>& origElemNums,
                                        StdVector<int>& elemTypes,
                                        UInt& elemRangeStart) 
  {
    UInt numElems = elems.GetSize();
    Elem::FEType feType;
    
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
      
      TranslateConnectivity(feType, &elemData[idx+1], ptEl->connect);

      idx += numElemNodes+1;

      regionIds[elemRangeStart-1+iElem] = ptEl->regionId;
      origElemNums[elemRangeStart-1+iElem] = ptEl->elemNum;
      elemTypes[elemRangeStart-1+iElem] = feType;
    }

    int indexSection;
    cgsize_t nelemStart = elemRangeStart, nelemEnd = elemRangeStart+numElems-1;
    int nbdyElem=0;
    cg_section_write(indexFile_,indexBase_,indexZone_,name.c_str(),MIXED,nelemStart,
                     nelemEnd,nbdyElem,&elemData[0],&indexSection);

    elemRangeStart = nelemEnd+1;
  }
  
  void SimOutputCGNS::WritePureSection(const StdVector<Elem*>& elems,
                                       const std::string& name,
                                       StdVector<int>& regionIds,
                                       StdVector<int>& origElemNums,
                                       StdVector<int>& elemTypes,
                                       UInt& elemRangeStart)
  {
    UInt numElems = elems.GetSize();
    Elem::FEType feType = elems[0]->type;
    ElementType_t eType = elemTypeMap_[feType];
    UInt numElemNodes = Elem::shapes[feType].numVertices;
    if(writeQuadElems_) 
    {
      numElemNodes = Elem::shapes[feType].numNodes;
    }
    
    StdVector<cgsize_t> elemData(numElems*numElemNodes);

    for (UInt iElem=0; iElem<numElems; iElem++) {
      Elem* ptEl = elems[iElem];
      TranslateConnectivity(feType, &elemData[iElem*numElemNodes], ptEl->connect);

      regionIds[elemRangeStart-1+iElem] = ptEl->regionId;
      origElemNums[elemRangeStart-1+iElem] = ptEl->elemNum;
      elemTypes[elemRangeStart-1+iElem] = feType;
    }

    int indexSection;
    cgsize_t nelemStart = elemRangeStart, nelemEnd = elemRangeStart+numElems-1;
    int nbdyElem=0;
    cg_section_write(indexFile_,indexBase_,indexZone_,name.c_str(),eType,nelemStart,
                     nelemEnd,nbdyElem,&elemData[0],&indexSection);

    elemRangeStart = nelemEnd+1;
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
      elemTypeMap_[Elem::ET_PYRA13]  = CGNSLIB_H::PYRA_13;
      elemTypeMap_[Elem::ET_PYRA14]  = CGNSLIB_H::PYRA_14;
      elemTypeMap_[Elem::ET_WEDGE15] = CGNSLIB_H::PENTA_15;
      elemTypeMap_[Elem::ET_WEDGE18] = CGNSLIB_H::PENTA_18;
      elemTypeMap_[Elem::ET_HEXA20]  = CGNSLIB_H::HEXA_20;
      elemTypeMap_[Elem::ET_HEXA27]  = CGNSLIB_H::HEXA_27;
    }
    else 
    {
      WARN("Quadratic element type will be reduced to linear ones for CGNS!");

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
                                           StdVector<UInt>& connect)
  {
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
      cgnsConn[n] = connect[tr[n]];
    }
  }

  void  SimOutputCGNS::NodeElemDataTransient(const bool definedOnNode,
                                             std::map< std::string, Vector<Double> >& gSol,
                                             const UInt step, 
                                             const Double time)
  {
    std::map< std::string, Vector<Double> >::const_iterator it, end;
    int* indexSol = &indexNodeSol_;
    
    // Check if solution node is already present, and write it if necessary.
    if(!definedOnNode) 
    {
      indexSol = &indexElemSol_;
    }    

    // Write actual solution array.
    it = gSol.begin();
    end = gSol.end();
    for( ; it != end; it++ ) 
    {    
      const Vector<Double>& sol = it->second;
      std::string solName = it->first;
      int indexField;

      if(cg_field_write(indexFile_, indexBase_, indexZone_,
                        *indexSol, RealDouble, solName.c_str(),
                        &sol[0], &indexField)) 
      {
        cg_close(indexFile_);
        EXCEPTION("Cannot write field:\n" <<
                  cg_get_error());
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
    fileName = fileName / (fileName_ + std::string(".cgns"));
    
    // Make sure a previous file with the same name gets removed.
    try {
      if(fs::exists(fileName))
      {
        fs::remove( fileName );
      }
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    // Set type of mid level library to use.
    if(myParam_->Has("mll")) 
    {
      int file_type = CG_FILE_HDF5;

      std::string mll = myParam_->Get("mll")->As<std::string>();

      if(mll == "adf") 
      {
        file_type = CG_FILE_ADF;
      }
#if CGNS_VERSION > 3000
      if(mll == "adf2")
      {
        file_type = CG_FILE_ADF2;
      }
#endif
      
      cg_set_file_type(file_type);
    }    

    // Open CGNS file for write
    if (cg_open(fileName.string().c_str(),CG_MODE_WRITE,&indexFile_)) cg_error_exit();
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
      
      if(  actInfo.definedOn != ResultInfo::NODE &&
           actInfo.definedOn != ResultInfo::ELEMENT &&
           actInfo.definedOn != ResultInfo::SURF_ELEM ) {
        WARN( "CGNS can only write results on element and nodes." );
        continue;
      }
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;      
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
      bool definedOnNode = (actInfo.definedOn == ResultInfo::NODE);

      // Create map of arrays, which are then written to CGNS.
      std::map<std::string, Vector<Double> > gSol;
      FillGlobalVectors(gSol, actResults, entityType );
      NodeElemDataTransient( definedOnNode, gSol, actStep_, actStepVal_);
    } // over all result types
  }
    
  
  void SimOutputCGNS::FinishMultiSequenceStep() {
    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }

  void SimOutputCGNS::FillGlobalVectors(std::map< std::string, Vector<Double> >& gSol, 
                                        const StdVector<shared_ptr<BaseResult> > & solList,
                                        ResultInfo::EntityUnknownType entityType ) {
    static const double H180DEG_OVER_PI = 180.0 / PI;

    BaseMatrix::EntryType entryType = solList[0]->GetEntryType();

    std::map<UInt, UInt> entity2Idx;
    if( entityType == ResultInfo::NODE ) {
      for(UInt i=0, n=numNodes_; i<n; i++) 
      {
        entity2Idx[i+1] = i;
      }
    }
    else {
      UInt idx = 0;
      // Iterate over all regions
      for(UInt i=0, n=ptGrid_->regionData.GetSize(); i<n; i++) 
      {
        StdVector<Elem*> elems;
        ptGrid_->GetElems(elems, ptGrid_->regionData[i].id );
        for( UInt j=0, m=elems.GetSize(); j<m; j++ ) {
          entity2Idx[elems[j]->elemNum] = idx;
          idx++;
        }
      }
    }
    
    UInt numEntities = entity2Idx.size();
    ResultInfo & actResultInfo = *(solList[0]->GetResultInfo());
    std::string resultName = actResultInfo.resultName;
    StdVector<std::string> dofNames = actResultInfo.dofNames;
    UInt numDofs = actResultInfo.dofNames.GetSize();
    LOG_DBG(simOutputCGNS) << "Filling global vector for result '" 
                           << resultName << "' on ";
    for( UInt i = 0; i < solList.GetSize(); i++ ) {
      LOG_DBG(simOutputCGNS) << solList[i]->GetEntityList()->GetName();
    }

    if( entityType == ResultInfo::ELEMENT ||
        entityType == ResultInfo::SURF_ELEM ) {

      // === Element Results ===
      for(UInt iDof=0; iDof<numDofs; iDof++) 
      {
        std::string resultDofName = (numDofs == 1 ? resultName : resultName + "_" + dofNames[iDof]);

        if( entryType == BaseMatrix::DOUBLE ) {
          // real valued results
          gSol[resultDofName].Resize( numEntities );
          gSol[resultDofName].Init();
          for( UInt i = 0; i < solList.GetSize(); i++ ){
            EntityIterator it = solList[i]->GetEntityList()->GetIterator();
            Vector<Double> & actSol = dynamic_cast<Result<Double>&>
                                      (*(solList[i])).GetVector();
            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt elemNum = it.GetElem()->elemNum;
              gSol[resultDofName][entity2Idx[elemNum]] = actSol[it.GetPos()*numDofs+iDof];
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

          gSol[resNameRe].Resize( numEntities ); gSol[resNameRe].Init();
          gSol[resNameIm].Resize( numEntities ); gSol[resNameIm].Init();
          gSol[resNameAmpl].Resize( numEntities ); gSol[resNameAmpl].Init();
          gSol[resNamePhase].Resize( numEntities ); gSol[resNamePhase].Init();

          for( UInt i = 0; i < solList.GetSize(); i++ ){
            EntityIterator it = solList[i]->GetEntityList()->GetIterator();
            Vector<Complex> & actSol = dynamic_cast<Result<Complex>&>
                                      (*(solList[i])).GetVector();
            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt elemNum = it.GetElem()->elemNum;
              if(entity2Idx.find(elemNum) != entity2Idx.end()) 
              {
                UInt idx = entity2Idx[elemNum];
                Complex sol = actSol[it.GetPos()*numDofs+iDof];
                
                gSol[resNameRe][idx] = sol.real();
                gSol[resNameIm][idx] = sol.imag();
                gSol[resNameAmpl][idx] = hypot(sol.real(), sol.imag());
                gSol[resNamePhase][idx] = (std::abs(sol.imag()) > 1e-16) ?
                                          std::atan2( sol.imag(), sol.real() ) * H180DEG_OVER_PI : 
                                          ( sol.real() < 0.0 ) ? 180 : 0 ;
              }
            }
          }
        }
      }
    } else if ( entityType == ResultInfo::NODE ) {

      // === Nodal Results ===
      for(UInt iDof=0; iDof<numDofs; iDof++) 
      {
        std::string resultDofName = (numDofs == 1 ? resultName : resultName + "_" + dofNames[iDof]);
        
        if( entryType == BaseMatrix::DOUBLE ) {
          // real valued results
          gSol[resultDofName].Resize( numEntities );
          gSol[resultDofName].Init();
          for( UInt i = 0; i < solList.GetSize(); i++ ){
            EntityIterator it = solList[i]->GetEntityList()->GetIterator();
            Vector<Double> & actSol = dynamic_cast<Result<Double>&>
                                      (*(solList[i])).GetVector();
            assert( (UInt) (actSol.GetSize()/numDofs) 
                    == solList[i]->GetEntityList()->GetSize());
            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt nodeNum = it.GetNode();
              gSol[resultDofName][entity2Idx[nodeNum]] = actSol[it.GetPos()*numDofs+iDof];
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

          gSol[resNameRe].Resize( numEntities ); gSol[resNameRe].Init();
          gSol[resNameIm].Resize( numEntities ); gSol[resNameIm].Init();
          gSol[resNameAmpl].Resize( numEntities ); gSol[resNameAmpl].Init();
          gSol[resNamePhase].Resize( numEntities ); gSol[resNamePhase].Init();

          for( UInt i = 0; i < solList.GetSize(); i++ ){
            EntityIterator it = solList[i]->GetEntityList()->GetIterator();
            Vector<Complex> & actSol = dynamic_cast<Result<Complex>&>
                                      (*(solList[i])).GetVector();
            for( it.Begin(); !it.IsEnd(); it++ ) {
              UInt idx = entity2Idx[it.GetNode()];
              Complex sol = actSol[it.GetPos()*numDofs+iDof];

              gSol[resNameRe][idx] = sol.real();
              gSol[resNameIm][idx] = sol.imag();
              gSol[resNameAmpl][idx] = hypot(sol.real(), sol.imag());
              gSol[resNamePhase][idx] = (std::abs(sol.imag()) > 1e-16) ?
                                        std::atan2( sol.imag(), sol.real() ) * H180DEG_OVER_PI : 
                                        ( sol.real() < 0.0 ) ? 180 : 0 ;
            }
          }
        }
      }
      
    } else {
      EXCEPTION( "Can only map nodal / element results to grid" );
    }      
  }

}
