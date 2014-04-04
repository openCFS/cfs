// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iomanip>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/Logging/LogConfigurator.hh"

#include "cgnslib.h"
#if CGNS_VERSION < 3100
# define cgsize_t int
#else
# if CG_BUILD_SCOPE
#  error enumeration scoping needs to be off
# endif
#endif

#include "SimOutputCGNS.hh"

namespace CoupledField {


  // ***************
  //   Constructor
  // ***************
  SimOutputCGNS::SimOutputCGNS( std::string& fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode,
                                bool isRestart ) 
    : SimOutput ( fileName, outputNode, infoNode, isRestart ) {
    
    formatName_ = "cgns";
    stepNumOffset_ = 0;
    stepValOffset_ = 0.0;
    
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

    int geodim=ptGrid_->GetDim(), physdim=geodim;
    strcpy(baseName_,"CFS_Simulation");
    cg_base_write(indexFile_,baseName_,geodim,physdim,&indexBase_);

    WriteNodesAndElements();
  }

  void  SimOutputCGNS::WriteNodesAndElements() {

    if (!ptGrid_)
      EXCEPTION("ptGrid_ is not initialized" );

    UInt geodim=ptGrid_->GetDim();
    UInt numNodes=ptGrid_->GetNumNodes();

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
    isize[0][0] = numNodes;
    isize[1][0] = numElems;
    isize[2][0] = 0;
    cg_zone_write(indexFile_,indexBase_,zoneName_,isize[0],Unstructured,&indexZone_);

    // coordinates
    StdVector<double> xCoord(numNodes);
    StdVector<double> yCoord(numNodes);
    StdVector<double> zCoord(numNodes);

    Vector<Double> point;
    for ( UInt i = 0; i < numNodes; i++ ) {
      ptGrid_->GetNodeCoordinate(point,i+1);
      if ( geodim == 3 ) {
        xCoord[i] = point[0];
        yCoord[i] = point[1];
        zCoord[i] = point[2];
      } else if ( geodim == 2 ) {
        xCoord[i] = point[0];
        yCoord[i] = point[1];
      }
    }

    int indexCoord;
    // write grid coordinates
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateX",&xCoord[0],&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateY",&yCoord[0],&indexCoord);
    cg_coord_write(indexFile_,indexBase_,indexZone_,RealDouble,"CoordinateZ",&zCoord[0],&indexCoord);


    for(UInt i=0, n=ptGrid_->regionData.GetSize(); i<n; i++) 
    {
      RegionIdType regionId = ptGrid_->regionData[i].id;

      StdVector<Elem*> elems;
      
      ptGrid_->GetElems(elems, regionId);
      UInt nElems = elems.GetSize();
      Elem::FEType feType = elems[0]->type;
      UInt j=1;
      
      // Check if all elems in regions are of same type
      for( ; j<nElems; j++)
      {
        if(elems[j]->type != feType)
          break;
      }

      if(j < nElems) 
      {
        // Write mixed section.
        WriteMixedSection(elems, ptGrid_->regionData[regionId].name);
      }
      else
      {
        // Write pure section.
        WritePureSection(elems, ptGrid_->regionData[regionId].name);
      }
    }

    if (1 == 0) {
    std::string solname = "DummyPressure1";
    double *press = new double[numNodes];
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 1.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    int indexFlow,indexField;
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 2.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure2";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
    for (UInt node=0; node<numNodes; node++) {
        press[node] = 3.0+xCoord[node]*(1.0 - yCoord[node]*yCoord[node])*exp(1.0-zCoord[node]);
    }
    solname = "DummyPressure3";
    cg_sol_write(indexFile_,indexBase_,indexZone_,solname.c_str(),Vertex,&indexFlow);
    cg_field_write(indexFile_,indexBase_,indexZone_,indexFlow,RealDouble,
            "Pressure",press,&indexField);
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
    }

    cg_close(indexFile_);
  }

  void SimOutputCGNS::WriteMixedSection(const StdVector<Elem*>& elems,
                                        const std::string name) 
  {
    UInt numElems = elems.GetSize();
    Elem::FEType feType;
    
    StdVector<cgsize_t> elemData(numElems*30);

    for (UInt iElem=0, idx=0; iElem<numElems; iElem++) {
      Elem* ptEl = elems[iElem];
      feType = ptEl->type;
      UInt numElemNodes = Elem::shapes[feType].numNodes;
      
      elemData[idx] = elemTypeMap_[feType];
      
      TranslateConnectivity(feType, &elemData[idx+1], ptEl->connect);

      idx += numElemNodes+1;
    }

    int indexSection;
    int nelemStart = 1, nelemEnd = numElems, nbdyElem=0;
    cg_section_write(indexFile_,indexBase_,indexZone_,name.c_str(),MIXED,nelemStart,
                     nelemEnd,nbdyElem,&elemData[0],&indexSection);
  }
  
  void SimOutputCGNS::WritePureSection(const StdVector<Elem*>& elems,
                                       const std::string name)
  {
    UInt numElems = elems.GetSize();
    Elem::FEType feType = elems[0]->type;
    UInt numElemNodes = Elem::shapes[feType].numNodes;
    
    StdVector<cgsize_t> elemData(numElems*numElemNodes);

    for (UInt iElem=0; iElem<numElems; iElem++) {
      Elem* ptEl = elems[iElem];
      TranslateConnectivity(feType, &elemData[iElem*numElemNodes], ptEl->connect);
    }

    int indexSection;
    int nelemStart = 1, nelemEnd = numElems, nbdyElem=0;
    cg_section_write(indexFile_,indexBase_,indexZone_,name.c_str(),elemTypeMap_[feType],nelemStart,
                     nelemEnd,nbdyElem,&elemData[0],&indexSection);
  }

  void SimOutputCGNS::InitElemTypeMap(){
    elemTypeMap_.clear();
    elemTypeMap_[Elem::ET_UNDEF]   = CGNSLIB_H::ElementTypeNull;
    elemTypeMap_[Elem::ET_POINT]   = CGNSLIB_H::NODE;
    elemTypeMap_[Elem::ET_LINE2]   = CGNSLIB_H::BAR_2;
    elemTypeMap_[Elem::ET_LINE3]   = CGNSLIB_H::BAR_3;
    elemTypeMap_[Elem::ET_TRIA3]   = CGNSLIB_H::TRI_3;
    elemTypeMap_[Elem::ET_TRIA6]   = CGNSLIB_H::TRI_6;
    elemTypeMap_[Elem::ET_QUAD4]   = CGNSLIB_H::QUAD_4;
    elemTypeMap_[Elem::ET_QUAD8]   = CGNSLIB_H::QUAD_8;
    elemTypeMap_[Elem::ET_QUAD9]   = CGNSLIB_H::QUAD_9;
    elemTypeMap_[Elem::ET_TET4]    = CGNSLIB_H::TETRA_4;
    elemTypeMap_[Elem::ET_TET10]   = CGNSLIB_H::TETRA_10;
    elemTypeMap_[Elem::ET_PYRA5]   = CGNSLIB_H::PYRA_5;
    elemTypeMap_[Elem::ET_PYRA13]  = CGNSLIB_H::PYRA_13;
    elemTypeMap_[Elem::ET_PYRA14]  = CGNSLIB_H::PYRA_14;
    elemTypeMap_[Elem::ET_WEDGE6]  = CGNSLIB_H::PENTA_6;
    elemTypeMap_[Elem::ET_WEDGE15] = CGNSLIB_H::PENTA_15;
    elemTypeMap_[Elem::ET_WEDGE18] = CGNSLIB_H::PENTA_18;
    elemTypeMap_[Elem::ET_HEXA8]   = CGNSLIB_H::HEXA_8;
    elemTypeMap_[Elem::ET_HEXA20]  = CGNSLIB_H::HEXA_20;
    elemTypeMap_[Elem::ET_HEXA27]  = CGNSLIB_H::HEXA_27;
  }

  void SimOutputCGNS::TranslateConnectivity(Elem::FEType feType,
                                           cgsize_t* cgnsConn,
                                           StdVector<UInt>& connect)
  {
    UInt numElemNodes = Elem::shapes[feType].numNodes;

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
      cgnsConn[tr[n]] = connect[n];
    }
  }

  void  SimOutputCGNS::NodeElemDataTransient(const UInt dataSetNr,
                                                  const std::string & title, 
                                                  const Vector<Double> & x, 
                                                  const UInt step, 
                                                  const Double time, 
                                                  const UInt nrNodes,
                                                  const UInt nrDofs)
  {
    //    EXCEPTION("NodeElemDataTransient not implemented!");
  }  

  void SimOutputCGNS::NodeElemDataHarmonic(const UInt dataSetNr,
                                                const std::string & title, 
                                                const Vector<Complex> & x, 
                                                const UInt step,
                                                const Double frequency, 
                                                const ComplexFormat format,
                                                const UInt nrNodes,
                                                const UInt nrDofs)
  {
    //    EXCEPTION("NodeElemDataHarmonic not implemented!");
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
      if(mll == "adf")
      {
        file_type = CG_FILE_ADF2;
      }
      
      cg_set_file_type(file_type);
    }    

    //     open CGNS file for write
    outputFileOK_ = false;
    if (cg_open(fileName.c_str(),CG_MODE_WRITE,&indexFile_)) cg_error_exit();
    outputFileOK_ = true;

    WriteGrid();
  }

  void SimOutputCGNS::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory ) {
    
  }

  void SimOutputCGNS::BeginStep( UInt stepNum, Double stepVal ) {
    
    actStep_ = stepNum + stepNumOffset_;
    actStepVal_ = stepVal + stepValOffset_;
    resultMap_.clear();
  }

  void SimOutputCGNS::AddResult( shared_ptr<BaseResult>  sol ) {
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputCGNS::FinishStep( ) {

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
        WARN( "Unv can only write results on element and nodes" );
        continue;
      }
      
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;
      std::string title =  SolutionTypeToString( actInfo.resultType );
      
      // if result could not be mapped, omit it
      if( title == "") {
        std::stringstream warning;
        warning  <<  "Result '" << actInfo.resultName
                 << "' coul not be mappted to unv result type and is "
                 << "omitted!";
        WARN( warning.str().c_str() );
        continue;
      }
      
      StdVector<std::string> & dofNames = actInfo.dofNames;
      UInt numDofs = dofNames.GetSize();
      
      // Determine type of entity the result is defined on
      UInt mapTo = 0;
      UInt numEntities = 0;
      if( actInfo.definedOn == ResultInfo::NODE ) {
        mapTo = 55;
        numEntities = ptGrid_->GetNumNodes();;
      } else {
        mapTo = 56;
        numEntities = ptGrid_->GetNumVolElems();
      }
      
      if( actResults[0]->GetEntryType() == BaseMatrix::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        
        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataTransient( mapTo, title, gSol, actStep_, 
                               actStepVal_ , numEntities ,
                               numDofs );
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );

        // Special case: If result is mechStress, restort entries 
        // according to capa notation
        if( actInfo.resultType == MECH_STRESS ) { 
          SortStresses( gSol, dofNames);
          numDofs = 6;
        }
        
        NodeElemDataHarmonic( mapTo, title, gSol, actStep_, 
                              actStepVal_, actInfo.complexFormat,
                              numEntities, numDofs );
      }
    } // over all result types
  }
  
  
  void SimOutputCGNS::FinishMultiSequenceStep() {

    // set offset for step value and number to last values
    stepNumOffset_ = actStep_;
    stepValOffset_ = actStepVal_;
  }
  
  
  template<class TYPE>
  void SimOutputCGNS::SortStresses( Vector<TYPE> & vec,

                                   const StdVector<std::string>& dofNames ){


    //soring according to capa (unv) notation
    //our notation is Voigt: Txx Tyy Tzz Tyz Txz Txy


    // check size of dofNames:
    // 6: 3D
    // 4: axi
    // 3: plane 

    // ensure, that number of entries is multiples of number of dofs
    assert( vec.GetSize() % dofNames.GetSize() == 0 );
    UInt numDofs = dofNames.GetSize();
    UInt numEntities = (UInt) vec.GetSize() / dofNames.GetSize();

    // CAPA always expects a stress vector with 6 entries
    Vector<TYPE> sorted;
    sorted.Resize( numEntities * 6 );
    
    // a) 3D-case (6 entries for stress vector)
    if( numDofs == 6 ) {
      
      //Txx Txy Tyy Txz Tyz Tzz
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=6 ) {
        sorted[j++] = vec[i];
        sorted[j++] = vec[i+5];
        sorted[j++] = vec[i+1];
        sorted[j++] = vec[i+4];
        sorted[j++] = vec[i+3];
        sorted[j++] = vec[i+2];
      }
    }
  
    // b) axisymmetric-case (4 entries for stress vector)  
    else if( numDofs == 4 ) {
      //Tphiphi 0 Trr 0 Trz Tzz
      
      
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=4 ) {
        sorted[j++] = vec[i+3];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      }
    }

    // c) planestress/strain-case (3 entries for stress vector)  
    else if( numDofs == 3 ) {
      //0 0 Tyy 0 Tyz Tzz
    
      UInt j = 0;
      for( UInt i = 0; i < vec.GetSize(); i+=3 ) {
        sorted[j++] = 0.0;
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+0];
        sorted[j++] = 0.0;
        sorted[j++] = vec[i+2];
        sorted[j++] = vec[i+1];
      } 

    } else {
      EXCEPTION( "Can not resort a stress vector with " << numDofs << " entries" );
    }
   
    // store sorted vector back to original one
    vec = sorted;
  }
  
  std::string SimOutputCGNS::SolutionTypeToString(const SolutionType type) const
  {

    //   std::string warnMsg;

    switch (type)
      {
      case MECH_DISPLACEMENT:
        return "displacement";
        break;
      case MECH_ACCELERATION:
        return "acceleration";
        break;
      case MECH_VELOCITY:
        return "velocity";
        break;
      case MECH_FORCE:
        break;
      case MECH_STRESS:
        return "stress";
        break;
      case MECH_STRAIN:
        EXCEPTION("Not implemented" );
        break;
      case MECH_RHS_LOAD:
        return "mechRhsLoad";
         break;
      case ELEC_POTENTIAL:
        return "electric potential";
        break;
      case ELEC_FIELD_INTENSITY:
        return "electric field";
        break;
      case ELEC_FORCE_VWP: 
        EXCEPTION("Not implemented" );
        break;
      case ELEC_CHARGE:
        return "electric charge";
        break;
      case ELEC_FLUX_DENSITY:
        EXCEPTION("Not implemented");
        break; 
      case ELEC_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case ELEC_RHS_LOAD:
        return "elecRhsLoad";
        break;
      case SMOOTH_DISPLACEMENT:
        return "displacement";
        break;
      case ACOU_POTENTIAL:
        return "fluid potential";
        break;

      case ACOU_RHS_LOAD:
        return "fluid potential";
        break;
      case ACOU_PRESSURE:
        //       warnMsg = "Due to the restrictions in the .unv file format, the ";
        //       warnMsg += "acoustic pressure is written as acoustic (fluid) potential!";
        //       WARN(warnMsg.c_str(), __FILE__, __LINE__);
        return "fluid potential";
        break;
      case ACOU_FORCE:
        EXCEPTION("Not implemented" );
        break;
      case ACOU_POTENTIAL_DERIV_1:
        return "fluid potential, 1st deriv.";
        break;
      case ACOU_POTENTIAL_DERIV_2:
        return "fluid potential, 2nd deriv.";
        break;
      case MAG_POTENTIAL:
        return "mag. vector potential";
        break;
      case MAG_FLUX_DENSITY:
        return "mag. flux density";
        break;
      case MAG_FORCE_VWP:
        EXCEPTION("Not implemented");
        break;
      case MAG_FORCE_LORENTZ:
        EXCEPTION("Not implemented");
        break;
      case MAG_ENERGY:
        EXCEPTION("Not implemented");
        break;
      case MAG_RHS_LOAD:
        return "magRhsLoad";
        break;

      case HEAT_TEMPERATURE:
        return "temperature";
        break;
      case HEAT_RHS_LOAD:
        return "temperatureRhsLoad";
        break;

      default:
        WARN( "Wrong type of solution or 'SolutionType2String'" );
        return "UNKNOWN";
        break;
      }
    return std::string();
  }

}
