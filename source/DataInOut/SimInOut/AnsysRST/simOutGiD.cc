// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <complex>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/elem.hh"
#include "DataInOut/Logging/cfslog.hh"

#include "simOutGiD.hh"


namespace CoupledField {
  
  // declare logging stream
  DECLARE_LOG(simOutputGiD)
  DEFINE_LOG(simOutputGiD, "SimOutputGiD")
 
    SimOutputGiD::SimOutputGiD( const std::string& fileName,
                                ParamNode * outputNode )
      : SimOutput( fileName, outputNode ) {

    ENTER_FCN( "SimOutputGiD::SimOutputGiD" );

    // Initialize variables
    formatName_ = "gid";
    fileName_ = fileName;
    dirName_ = ".";
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    dim_ = GiD_2D;
    isAscii_ = true;
    degen3DElems_ = true;

    // Determine, if binary result file should be written
    isAscii_ = !(myParam_->Get("binaryFormat")->AsBool() );

    // Open result file
    std::string postFileName = fileName_;
    if ( isAscii_ == true) {
      postFileName += ".post.res";
      GiD_OpenPostResultFile( postFileName.c_str(), GiD_PostAscii );
    } else {
      postFileName += ".post.bin";
      GiD_OpenPostResultFile( postFileName.c_str(), GiD_PostBinary );
    }
  }


  SimOutputGiD::~SimOutputGiD() {

    ENTER_FCN( "SimOutputGiD::~SimOutputGiD" );

    // Close result file
    GiD_ClosePostResultFile();
  }


  void SimOutputGiD::WriteGrid() {
    ENTER_FCN( "SimOutputGiD::WriteGrid" );
    
    LOG_TRACE(simOutputGiD) << "Writing mesh";
    
    // open mesh file (only needed in ASCII case
    if ( isAscii_ == true) {
      std::string meshFileName = fileName_ + ".post.msh";
      GiD_OpenPostMeshFile( meshFileName.c_str(), GiD_PostAscii );
    }

    // write nodes
    WriteNodes();
    
    // write elements
    WriteElements();

    // write definition for Gauss points
    if ( dim_ == GiD_2D ) {
      GiD_BeginGaussPoint( "mid-2D", GiD_Quadrilateral, NULL, 1, 1, 1);
      GiD_EndGaussPoint();
      if ( ptGrid_->GetNumSurfElems() > 0 ) {
        GiD_BeginGaussPoint( "mid-2D", GiD_Linear, NULL, 1, 0, 1);
        GiD_EndGaussPoint();
      }
    } else {
      if ( degen3DElems_ == true ) {
        GiD_BeginGaussPoint( "mid-3D", GiD_Hexahedra, NULL, 1, 1, 1);
      } else {
        GiD_BeginGaussPoint( "mid-3D", GiD_Tetrahedra, NULL, 1, 1, 1);
      }
      GiD_EndGaussPoint();
      if ( ptGrid_->GetNumSurfElems() > 0 ) {
        GiD_BeginGaussPoint( "mid-3D", GiD_Quadrilateral, NULL, 1, 1, 1);
        GiD_EndGaussPoint();
      }
    }
      
    // close mesh file (only needed in ASCII case)
    if ( isAscii_ == true ) {
      GiD_ClosePostMeshFile();
    }
    LOG_TRACE(simOutputGiD)<< "Finished writing of mesh" << std::endl;

  }

  void SimOutputGiD::WriteNodes() {
    ENTER_FCN( "SimOutputGiD::WriteNodes" );

    LOG_TRACE(simOutputGiD) << "Writing nodes";
    
    // get number of nodes
    UInt numNodes = ptGrid_->GetNumNodes();

    // write nodal coordinates
    GiD_BeginMesh( "Nodes", dim_, GiD_Point, 1 );

    GiD_BeginCoordinates();    
    if ( dim_ == GiD_2D ) {

      Point point;
      for ( UInt i = 1; i <= numNodes; i++ ) {
        ptGrid_->GetNodeCoordinate(point,i);
        GiD_WriteCoordinates2D(i, point[0], point[1]);
      }
    } else {
      Point point;
      for ( UInt i = 1; i <= numNodes; i++ ) {
        ptGrid_->GetNodeCoordinate(point,i);
        GiD_WriteCoordinates(i, point[0], point[1], point[2]);
      } 
    }
    GiD_EndCoordinates();
    
    // dummy call for writing elements
    GiD_BeginElements();    
    GiD_EndElements();
    
    GiD_EndMesh();

    LOG_TRACE(simOutputGiD) << "Finished writing nodes" << std::endl;

  }

  void SimOutputGiD::WriteElements() {
    ENTER_FCN( "SimOutputGiD::WriteElements" );
   
    LOG_TRACE(simOutputGiD) << "Writing elements";
    
    // iterate over all regions
    StdVector<RegionIdType> regionIds;
    StdVector<Elem*> elemVec;
    ptGrid_->GetRegionIds( regionIds );

    for ( UInt iReg = 0; iReg < regionIds.GetSize(); iReg++ ) {

      // get region name and elements
      std::string regionName = ptGrid_->RegionIdToName(regionIds[iReg]);
      ptGrid_->GetElems(elemVec,regionIds[iReg]);

    
      GiD_ElementType eType = GiD_NoElement;
      UInt numElemNodes = 0;

      // determine element type and number of nodes
      if ( elemVec[0]->ptElem->GetDim() == 1 ) {
        eType = GiD_Linear;
        if ( ptGrid_->IsQuadratic() == true ) {
          numElemNodes = 3;
        } else {
          numElemNodes = 2;
        }
      } else if ( elemVec[0]->ptElem->GetDim() == 2 ) {
        eType = GiD_Quadrilateral;
        if ( ptGrid_->IsQuadratic() == true ) {
          numElemNodes = 8;
        } else {
          numElemNodes = 4;
        }
      } else if ( elemVec[0]->ptElem->GetDim() == 3 ) {
        if ( degen3DElems_ == true ) {
          eType = GiD_Hexahedra;
          if ( ptGrid_->IsQuadratic() == true ) {
            numElemNodes = 20;
          } else {
            numElemNodes = 8;
          }
        } else {
          eType = GiD_Tetrahedra;
          if ( ptGrid_->IsQuadratic() == true ) {
            numElemNodes = 10;
          } else {
            numElemNodes = 4;
          }
        }
      }
      
      GiD_BeginMesh( regionName.c_str(), dim_, eType, numElemNodes );
    
      // write dummy coordinate section
      GiD_BeginCoordinates();
      GiD_EndCoordinates();
      
      // write alle element declarations
      GiD_BeginElements();
      for ( UInt iElem = 0; iElem < elemVec.GetSize(); iElem++ ) {
        Elem * ptEl = elemVec[iElem];
    
        WriteElement( ptEl );
      }
      GiD_EndElements();

      GiD_EndMesh();
    
    }

    LOG_TRACE(simOutputGiD) << "Writing named elements/nodes";
    
    // write named nodes to mesh file
    // Note:  Named nodes are written as 'point' elements, so we need
    // for each set of named nodes a separate MESH section
    StdVector<std::string> nodeNames, elemNames;
    StdVector<UInt> nodeNumbers, elemNumbers;

    ptGrid_->GetListNodeNames(nodeNames);
    
    // check if there are any named nodes / elems in the grid
    if( nodeNames.GetSize() == 0 ) {
      return;
    }

    // get number of 'normal' regions within the grid
    UInt numRegions = regionIds.GetSize();
    
    // remember number of 'normal' elements
    UInt numelem = ptGrid_->GetNumElems() + 1;
    
    for( UInt i = 0; i < nodeNames.GetSize(); i++ ) {

      // get named nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);
      GiD_BeginMesh( nodeNames[i].c_str(), dim_, GiD_Point, 1 );

      // write empty coordinate section
      GiD_BeginCoordinates();
      GiD_EndCoordinates();

      // write 'nodal' elements
      GiD_BeginElements();
      
      for ( UInt iNode = 0; iNode < nodeNumbers.GetSize(); iNode ++ ) {
        UInt nodeNumber[2] = {nodeNumbers[iNode],i+numRegions+1};
        GiD_WriteElementMat( numelem++, (int*)(&nodeNumber) ); 
      }
      GiD_EndElements();
      GiD_EndMesh();
    }

    LOG_TRACE(simOutputGiD) << "Finished writing elements" << std::endl;
  }
  
  void SimOutputGiD::
  WriteElement( Elem * ptEl ) {
    
    // determine element dimension
    BaseFE * ptFE = ptEl->ptElem;
    StdVector<UInt> const & connect = ptEl->connect;
    StdVector<UInt> connectM;
    
    // --- LINEAR ELEMENTS ---
    if ( ptGrid_->IsQuadratic() == false ) {
      if ( ptFE == ptTr1 ) {
        connectM = connect;
        connectM.Push_back(connectM[2]);
      } else  if ( ptFE == ptTet1 ) {
        if ( degen3DElems_ == true ) {
          connectM.Resize(8);
          connectM[0]= connect[0];
          connectM[1]= connect[1];
          connectM[2]= connect[2];
          connectM[3]= connect[2];
          connectM[4]= connect[3];
          connectM[5]= connect[3];
          connectM[6]= connect[3];
          connectM[7]= connect[3];
        } else {
          connectM = connect;
        }
      } else  if ( ptFE == ptPyra1 ) {
        connectM.Resize(8);
        connectM[0]= connect[0];
        connectM[1]= connect[1];
        connectM[2]= connect[2];
        connectM[3]= connect[3];
        connectM[4]= connect[4];
        connectM[5]= connect[4];
        connectM[6]= connect[4];
        connectM[7]= connect[4];
      } else  if ( ptFE == ptWedge1 ) {
        connectM.Resize(8);
        connectM[0]= connect[0];
        connectM[1]= connect[1];
        connectM[2]= connect[2];
        connectM[3]= connect[2];
        connectM[4]= connect[3];
        connectM[5]= connect[4];
        connectM[6]= connect[5];
        connectM[7]= connect[5];
        
      } else {
        connectM = connect;
      }
    } else {
      // --- QUADRATIC ELEMENTS ---
      if ( ptFE == ptTr2 ) {
        connectM.Resize(8);
        connectM[0]= connect[0];
        connectM[1]= connect[1];
        connectM[2]= connect[2];
        connectM[3]= connect[2];
        connectM[4]= connect[3];
        connectM[5]= connect[4];
        connectM[6]= connect[2];
        connectM[7]= connect[5];
      } else  if ( ptFE == ptTet2 ) {
        if ( degen3DElems_ == true ) {
          connectM.Resize(20);
          connectM[0]= connect[0];
          connectM[1]= connect[1];
          connectM[2]= connect[2];
          connectM[3]= connect[2];
          connectM[4]= connect[3];
          connectM[5]= connect[3];
          connectM[6]= connect[3];
          connectM[7]= connect[3];
          connectM[8]= connect[4];
          connectM[9]= connect[5];
          connectM[10]= connect[6];
          connectM[11]= connect[2];
          connectM[12]= connect[7];
          connectM[13]= connect[8];
          connectM[14]= connect[9];
          connectM[15]= connect[2];
          connectM[16]= connect[3];
          connectM[17]= connect[3];
          connectM[18]= connect[3];
          connectM[19]= connect[3];
        } else {
          connectM = connect;
        }
        // NOTE: The numberings of hexehadras in gid differs for
        // the quadratic case!
      } else  if ( ptFE == ptHexa2 ) {
        connectM=connect;
        connectM[12] = connect[16];
        connectM[13] = connect[17];
        connectM[14] = connect[18];
        connectM[15] = connect[19];
        connectM[16] = connect[12];
        connectM[17] = connect[13];
        connectM[18] = connect[14];
        connectM[19] = connect[15];
      } else  if ( ptFE == ptPyra2 ) {
        connectM.Resize(20);
        connectM[0]= connect[0];
        connectM[1]= connect[1];
        connectM[2]= connect[2];
        connectM[3]= connect[3];
        connectM[4]= connect[4];
        connectM[5]= connect[4];
        connectM[6]= connect[4];
        connectM[7]= connect[4];
        connectM[8]= connect[5];
        connectM[9]= connect[6];
        connectM[10]= connect[7];
        connectM[11]= connect[8];
        connectM[12]= connect[9];
        connectM[13]= connect[10];
        connectM[14]= connect[11];
        connectM[15]= connect[12];
        connectM[16]= connect[4];
        connectM[17]= connect[4];
        connectM[18]= connect[4];
        connectM[19]= connect[4];
      } else  if ( ptFE == ptWedge2 ) {
        connectM.Resize(20);
        connectM[0]= connect[0];
        connectM[1]= connect[1];
        connectM[2]= connect[2];
        connectM[3]= connect[2];
        connectM[4]= connect[3];
        connectM[5]= connect[4];
        connectM[6]= connect[5];
        connectM[7]= connect[5];
        connectM[8]= connect[6];
        connectM[9]= connect[7];
        connectM[10]= connect[8];
        connectM[11]= connect[2];
        connectM[12]= connect[12];
        connectM[13]= connect[13];
        connectM[14]= connect[14];
        connectM[15]= connect[14];
        connectM[16]= connect[9];
        connectM[17]= connect[10];
        connectM[18]= connect[11];
        connectM[19]= connect[5];
      } else {
        connectM = connect;
      }
    }
    connectM.Push_back(ptEl->regionId+1);
    GiD_WriteElementMat( ptEl->elemNum, (int*)(connectM.GetPointer()) ); 
  }

 
  void SimOutputGiD::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd ) {
    ENTER_FCN( "SimOutputGiD::RegisterResult" );
    
    ResultInfo & actDof = *(sol->GetResultInfo());
    LOG_DBG(simOutputGiD) << "Registering output '" << actDof.resultName
                          << "' with saveBegin " << saveBegin
                          << ", saveInc " << saveInc 
                          << ", saveEnd " << saveEnd;

  }

  void SimOutputGiD::BeginStep( UInt stepNum, Double stepVal ) {
    ENTER_FCN( "SimOutputGiD::BeginStep" );

    actStep_ = stepNum;
    actStepVal_ = stepVal;
    resultMap_.clear();
  }
 
  void SimOutputGiD::AddResult( shared_ptr<BaseResult> sol ) {
    ENTER_FCN( " SimOutputGiD::AddResult" );
    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGiD) << "Adding result '" 
                          << actDof.resultName  << "'";
      
    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }
  
  void SimOutputGiD::FinishStep( ) {
    ENTER_FCN( "SimOutputGiD::FinishStep" );

    LOG_TRACE(simOutputGiD) << "Starting to finish Step";
    
    std::string title;
    
    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // get result info object and results for current result type
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      title = actInfo.resultName;
      if( actInfo.unit.size() != 0 ) {
        title += " (";
        title += actInfo.unit;
        title += ")";
      }
      ResultInfo::EntryType entryType =  actInfo.entryType;
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;

      // check if result is defined on nodes or elements
      StdVector<std::string> & dofNames = actInfo.dofNames;  
      if(  actInfo.definedOn != ResultInfo::NODE &&
           actInfo.definedOn != ResultInfo::PFEM &&
           actInfo.definedOn != ResultInfo::ELEMENT &&
           actInfo.definedOn != ResultInfo::SURF_ELEM ) {
        Warning( "GiD can only write results on element and nodes",
                 __FILE__, __LINE__ );
        continue;
      }

      LOG_DBG(simOutputGiD) << "Writing result '" << title << "'";
      
      if( actResults[0]->GetEntryType() == EntryType::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        WriteNodeElemDataTrans( gSol, dofNames, title, entryType,
                                entityType, actStepVal_ );        
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );
        WriteNodeElemDataHarm( gSol, dofNames, title, entryType, entityType,
                               actStepVal_, actInfo.complexFormat );        
      }
      
    }
  }
  

  void SimOutputGiD::
  WriteNodeElemDataTrans( const Vector<Double> & var, 
                          const StdVector<std::string> & dofNames,
                          const std::string& name, 
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {
    ENTER_FCN ( "SimOutputGiD::WriteNodeElemDataTrans" );
    
    // get number of entities
    UInt numEnt = 0;
    char * dummy;
    GiD_ResultLocation loc;
    if ( entityType == ResultInfo::NODE ||
         entityType == ResultInfo::PFEM ) {
      loc = GiD_OnNodes;
      numEnt = ptGrid_->GetNumNodes();
      dummy = NULL;
    } else {
      loc = GiD_OnGaussPoints;
      numEnt = ptGrid_->GetNumElems();
      if ( dim_ == GiD_2D) {
        dummy = "mid-2D"; 
      } else {
        dummy = "mid-3D";
      }
    }
    
    
    UInt numDofs = dofNames.GetSize();
    
    // write Result header
    if ( entryType == ResultInfo::SCALAR )  {
      GiD_BeginResultHeader( name.c_str(), "transient", time,
                             GiD_Scalar, loc, dummy );
      GiD_ResultValues();
      for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
        GiD_WriteScalar( iEnt, var[iEnt-1]); 
      }
    } else if( entryType == ResultInfo::VECTOR ) {
      GiD_BeginResultHeader( name.c_str(), "transient", time,
                             GiD_Vector, loc, dummy );
      const char *names[3] = {"___", "__", "_"};
      for( UInt i = 0; i < numDofs; i++ ) {
        names[i] = dofNames[i].c_str();
      }
      GiD_ResultComponents( 3, names );
      
      // write results
      GiD_ResultValues();
      UInt offset = 0;
      
      for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
        offset = (iEnt-1) * numDofs;
        
        switch( numDofs ) {
        case 1:
          GiD_WriteVector( iEnt, var[offset], 0.0, 0.0);
          break;
        case 2:
          GiD_WriteVector( iEnt, var[offset], var[offset+1], 0.0);
          break;
        case 3:
          GiD_WriteVector( iEnt, var[offset], var[offset+1], var[offset+2] );
          break;
        default:
          EXCEPTION( "GiD can only write vectors with 3 entries. Your result has "
                     << numDofs << " components!");
        }
      }
      
  } else if( entryType == ResultInfo::TENSOR ) {
    GiD_BeginResultHeader( name.c_str(), "transient", time,
                           GiD_Matrix, loc, dummy );
      const char *names[6] = {"______", "_____", "____", 
                              "___", "__", "_" };
      for( UInt i = 0; i < numDofs; i++ ) {
        names[i] = dofNames[i].c_str();
      }
      GiD_ResultComponents( 6, names );
      GiD_ResultValues();
      
      // write results
      UInt offset = 0;
#define LOOP(e1, e2, e3, e4, e5, e6 )                   \
for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {         \
  offset = (iEnt-1) * numDofs;                          \
  GiD_Write3DMatrix( iEnt, e1, e2, e3, e4 ,e5, e6 );    \
}
      if ( numDofs == 3 ) {
        LOOP( var[offset], var[offset+1], var[offset+2], 
              0.0, 0.0, 0.0 );
      } else if( numDofs == 4 ) {
        LOOP( var[offset], var[offset+1], var[offset+2], 
              var[offset+3], 0.0, 0.0 );
      } else if( numDofs == 5 ) {
        LOOP( var[offset], var[offset+1], var[offset+2], 
              var[offset+3], var[offset+4], 0.0 );
      } else if( numDofs == 6 ) {
        LOOP( var[offset], var[offset+1], var[offset+2], 
              var[offset+3], var[offset+4],  var[offset+5] );
      }
#undef LOOP
    }
    
    GiD_EndResult();
  }     
  
  void SimOutputGiD::
  WriteNodeElemDataHarm( const Vector<Complex> & var, 
                         const StdVector<std::string> & dofNames,
                         const std::string name, 
                         const ResultInfo::EntryType entryType,
                         ResultInfo::EntityUnknownType entityType,
                         const Double freq, 
                         const ComplexFormat outputFormat ) {
    ENTER_FCN ( "SimOutputGiD::WriteNodeElemDataHarm" );
    
   // get number of entities
    UInt numEnt = 0;
    char * dummy;
    GiD_ResultLocation loc;
    if ( entityType == ResultInfo::NODE ||
         entityType == ResultInfo::PFEM ) {
      loc = GiD_OnNodes;
      numEnt = ptGrid_->GetNumNodes();
      dummy = NULL;
    } else {
      loc = GiD_OnGaussPoints;
      numEnt = ptGrid_->GetNumElems();
      if ( dim_ == GiD_2D) {
        dummy = "mid-2D"; 
      } else {
        dummy = "mid-3D";
      }
    }
    
    // determine complete name of output quantities
    std::string outName1, outName2;
    if (outputFormat == REAL_IMAG) {
      outName1 = name + "-real";
      outName2 = name + "-imag";
    } else {      
      outName1 = name + "-amp";
      outName2 = name + "-phase";
    }
    
    UInt numDofs = dofNames.GetSize();
    
    // === Scalar entries ===
    if ( entryType == ResultInfo::SCALAR ) {
      
      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), "harmonic", freq,
                             GiD_Scalar, loc, dummy );
      GiD_ResultValues();
      
      if (outputFormat == REAL_IMAG) {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          GiD_WriteScalar( iEnt, var[iEnt-1].real()); 
        }
      } else {
        
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          GiD_WriteScalar( iEnt, std::abs(var[iEnt-1])); 
        }
      }
      GiD_EndResult();

      // --- Second component ---
      GiD_BeginResultHeader( outName2.c_str(), "harmonic", freq,
                             GiD_Scalar, loc, dummy );
      GiD_ResultValues();
      
      if (outputFormat == REAL_IMAG) {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          GiD_WriteScalar( iEnt, var[iEnt-1].imag()); 
        }
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          if ( std::abs(var[iEnt-1].imag()) > 1e-16 ) {
            GiD_WriteScalar( iEnt, CPhase( var[iEnt-1] ) );
          } else {
            GiD_WriteScalar( iEnt, 0.0); 
          }
        }
      }
      GiD_EndResult();
      
    } else if ( entryType == ResultInfo::VECTOR ) {
      // === Vectorial entries ===
      
      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), "harmonic", freq,
                             GiD_Vector, loc, dummy );
      const char *names[] = {"___", "__", "_"};
      for( UInt i = 0; i < numDofs; i++ ) {
        names[i] = dofNames[i].c_str();
      }
      GiD_ResultComponents( 3, names );
      GiD_ResultValues();
      
      UInt offset = 0;
      if (outputFormat == REAL_IMAG) {
        
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          switch( numDofs ) {
          case 1:
            GiD_WriteVector( iEnt, var[offset].real(), 
                             0.0, 0.0);
            break;
          case 2:
            GiD_WriteVector( iEnt, var[offset].real(), 
                             var[offset+1].real(), 0.0);
            break;
          case 3:
            GiD_WriteVector( iEnt, var[offset].real(), 
                             var[offset+1].real(), 
                             var[offset+2].real() );
            break;
          default:
            EXCEPTION( "GiD can only write vectors with 3 entries. Your result has "
                       << numDofs << " components!");
          }
        }
        GiD_EndResult(); 
        
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          switch( numDofs ) {
          case 1:
            GiD_WriteVector( iEnt, 
                             std::abs(var[offset]), 
                             0.0, 0.0 );
            break;
          case 2:
            GiD_WriteVector( iEnt, 
                             std::abs(var[offset]), 
                             std::abs(var[offset+1]), 
                             0.0 );
            break;
          case 3:
            GiD_WriteVector( iEnt, std::abs(var[offset]),
                             std::abs(var[offset+1]), 
                             std::abs(var[offset+2]) );
            break;
          default:
            EXCEPTION( "GiD can only write vectors with 3 entries. Your result has "
                       << numDofs << " components!");

          }
        }
        GiD_EndResult(); 
      }
      
      // --- Second component ---
      GiD_BeginResultHeader( outName2.c_str(), "harmonic", freq,
                             GiD_Vector, loc, dummy );
      GiD_ResultComponents( 3, names );
      GiD_ResultValues();
      
      offset = 0;
      if (outputFormat == REAL_IMAG) {
        
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          if ( numDofs == 3 ) {
            GiD_WriteVector( iEnt, var[offset].imag(), 
                             var[offset+1].imag(), 
                             var[offset+2].imag() );
          } else { 
            GiD_WriteVector( iEnt, var[offset].imag(), 
                             var[offset+1].imag(), 
                             0.0 );
          }
        }
        GiD_EndResult(); 
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          if ( numDofs == 3 ) {
            GiD_WriteVector( iEnt, CPhase( var[offset] ), 
                             CPhase( var[offset+1] ), 
                             CPhase( var[offset+2] ) );
            
          } else { 
            GiD_WriteVector( iEnt, CPhase( var[offset] ),
                             CPhase( var[offset+1]),
                             0.0 );
          }
        }
        GiD_EndResult(); 
      }
    } else if ( entryType == ResultInfo::TENSOR ) {
      // === Vectorial entries ===
      
      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), "harmonic", freq,
                             GiD_Matrix, loc, dummy );
      const char *names[6] = {"______", "_____", "____", 
                              "___", "__", "_" };
      for( UInt i = 0; i < numDofs; i++ ) {
        names[i] = dofNames[i].c_str();
      }
      GiD_ResultComponents( 6, names );
      GiD_ResultValues();
      
     
#define LOOP(e1, e2, e3, e4, e5, e6 )                   \
for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {         \
  offset = (iEnt-1) * numDofs;                          \
  GiD_Write3DMatrix( iEnt, e1, e2, e3, e4 ,e5, e6 );    \
 }
      
      UInt offset = 0;
      if (outputFormat == REAL_IMAG) {
        if( numDofs == 3 ) {
          LOOP( var[offset].real(), var[offset+1].real(), 
                var[offset+2].real(), 0.0 , 0.0, 0.0 );
        } else if( numDofs == 4 ) {
          LOOP( var[offset].real(), var[offset+1].real(), 
                var[offset+2].real(), var[offset+3].real(),
                0.0 , 0.0 );
        } else if( numDofs == 5 ) {
          LOOP( var[offset].real(), var[offset+1].real(), 
                var[offset+2].real(), var[offset+3].real(),
                var[offset+4].real(), 0.0 );
        } else if( numDofs == 6 ) {
          LOOP( var[offset].real(), var[offset+1].real(), 
                var[offset+2].real(), var[offset+3].real(),
                var[offset+4].real(), var[offset+5].real() );
        }
        GiD_EndResult(); 
        
      } else {
        if( numDofs == 3 ) {
          LOOP( std::abs(var[offset]), std::abs(var[offset+1]),
                std::abs(var[offset+2]), 0.0, 0.0, 0.0 );
        } else if( numDofs == 4 ) {
          LOOP( std::abs(var[offset]), std::abs(var[offset+1]),
                std::abs(var[offset+2]), std::abs(var[offset+3]),
                0.0, 0.0 );
        } else if( numDofs == 5 ) {
          LOOP( std::abs(var[offset]), std::abs(var[offset+1]),
                std::abs(var[offset+2]), std::abs(var[offset+3]),
                std::abs(var[offset+4]), 0.0);
        } else if( numDofs == 6 ) {
          LOOP( std::abs(var[offset]), std::abs(var[offset+1]),
                std::abs(var[offset+2]), std::abs(var[offset+3]),
                std::abs(var[offset+4]), std::abs(var[offset+5]) );
        }
        GiD_EndResult(); 
        
      }
      
      // --- Second component ---
      GiD_BeginResultHeader( outName2.c_str(), "harmonic", freq,
                             GiD_Vector, loc, dummy );
      GiD_ResultComponents( 6, names );
      GiD_ResultValues();
      
      offset = 0;
      if (outputFormat == REAL_IMAG) {
        if( numDofs == 3 ) {
          LOOP( var[offset].imag(), var[offset+1].imag(), 
                var[offset+2].imag(), 0.0 , 0.0, 0.0 );
        } else if( numDofs == 4 ) {
          LOOP( var[offset].imag(), var[offset+1].imag(), 
                var[offset+2].imag(), var[offset+3].imag(),
                0.0 , 0.0 );
        } else if( numDofs == 5 ) {
          LOOP( var[offset].imag(), var[offset+1].imag(), 
                var[offset+2].imag(), var[offset+3].imag(),
                var[offset+4].imag(), 0.0 );
        } else if( numDofs == 6 ) {
          LOOP( var[offset].imag(), var[offset+1].imag(), 
                var[offset+2].imag(), var[offset+3].imag(),
                var[offset+4].imag(), var[offset+5].imag() );
        }
      } else {
        if( numDofs == 3 ) {
          LOOP( CPhase(var[offset]), CPhase(var[offset+1]),
                CPhase(var[offset+2]), 0.0, 0.0, 0.0 );
        } else if( numDofs == 4 ) {
          LOOP( CPhase(var[offset]), CPhase(var[offset+1]),
                CPhase(var[offset+2]), CPhase(var[offset+3]),
                0.0, 0.0 );
        } else if( numDofs == 5 ) {
          LOOP( CPhase(var[offset]), CPhase(var[offset+1]),
                CPhase(var[offset+2]), CPhase(var[offset+3]),
                CPhase(var[offset+4]), 0.0);
        } else if( numDofs == 6 ) {
          LOOP( CPhase(var[offset]), CPhase(var[offset+1]),
                CPhase(var[offset+2]), CPhase(var[offset+3]),
                CPhase(var[offset+4]), CPhase(var[offset+5]) );
        }
      }
      GiD_EndResult(); 
    }
    
  }
  
  // ********
  //   Init
  // ********
  void SimOutputGiD::Init( Grid* ptGrid) {
    
    ENTER_FCN( "SimOutputGiD::OpenFile" );

    ptGrid_ = ptGrid;
    // initialize history files
    //InitHistoryFiles();

    // determine dimension of grid
    UInt dim = ptGrid_->GetDim();
    if ( dim == 2 ) {
      dim_ = GiD_2D;
    } else if ( dim == 3 ) {
      dim_ = GiD_3D;
    } else {
      EXCEPTION( "SimOutputGiD: Grid dimension " << dim
               << " is not supported by GiD mesh format." );
    }

    // check, if only tetrahedra are present
    // -> only in this case we can treat tetrahedras as elements with 4 nodes.
    // In all other cases (hexa and tets mixed) tetrahedras are treated as
    // degenerated hexahedras.
    if ( dim_ == GiD_3D ) {
      if ( ptGrid_->GetNumElemOfType( ET_TET4 ) > 0  &&
           ptGrid_->GetNumElemOfType( ET_HEXA8 ) ==  0  &&
           ptGrid_->GetNumElemOfType( ET_PYRA5 ) ==  0  &&
           ptGrid_->GetNumElemOfType( ET_WEDGE6 ) ==  0  ) {
        degen3DElems_ = false;
      }
    }
    
  }
} // end of namespace
