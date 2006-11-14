#include <complex>

#include "outGiD.hh"

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/elem.hh"


namespace CoupledField {


 
  WriteResultsGiD::WriteResultsGiD( const Char *const filename)
    : WriteResults(filename) {

    ENTER_FCN( "WriteResultsGiD::WriteResultsGiD" );

    // Initialize variables
    dim_ = GiD_2D;
    isAscii_ = true;
    degen3DElems_ = true;

    // Determine, if binary result file should be written
    isAscii_ = !params->IsSet( "binaryFormat", "gid" );

    // Open result file
    std::string postFileName = namefile_;
    if ( isAscii_ == true) {
      postFileName += ".post.res";
      GiD_OpenPostResultFile( postFileName.c_str(), GiD_PostAscii );
    } else {
      postFileName += ".post.bin";
      GiD_OpenPostResultFile( postFileName.c_str(), GiD_PostBinary );
    }
  }


  WriteResultsGiD::~WriteResultsGiD() {

    ENTER_FCN( "WriteResultsGiD::~WriteResultsGiD" );

    // Close result file
    GiD_ClosePostResultFile();
  }


  void WriteResultsGiD::WriteGrid() {

    ENTER_FCN( "WriteResultsGiD::WriteGrid" );
    
    // open mesh file (only needed in ASCII case
    if ( isAscii_ == true) {
      std::string meshFileName = namefile_ + ".post.msh";
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
    } else {
      if ( degen3DElems_ == true ) {
        GiD_BeginGaussPoint( "mid-3D", GiD_Hexahedra, NULL, 1, 1, 1);
      } else {
        GiD_BeginGaussPoint( "mid-3D", GiD_Tetrahedra, NULL, 1, 1, 1);
      }
      GiD_EndGaussPoint();
    }
      
    // close mesh file (only needed in ASCII case)
    if ( isAscii_ == true ) {
      GiD_ClosePostMeshFile();
    }

  }

  void WriteResultsGiD::WriteNodes() {
    ENTER_FCN( "WriteResultsGiD::WriteNodes" );

  

    // get number of nodes
    UInt numNodes = ptGrid_->GetNumNodes();

    // write nodal coordinates
    GiD_BeginMesh( "Nodes", dim_, GiD_Point, 1 );

    GiD_BeginCoordinates();    
    if ( dim_ == GiD_2D ) {

      Point<2> point;
      for ( UInt i = 1; i <= numNodes; i++ ) {
        ptGrid_->GetNodeCoordinate(point,i);
        GiD_WriteCoordinates2D(i, point[0], point[1]);
      }
    } else {

      Point<3> point;
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

  }

  void WriteResultsGiD::WriteElements() {
    ENTER_FCN( "WriteResultsGiD::WriteElements" );
   
    // iterate over all regions
    StdVector<RegionIdType> regionIds;
    StdVector<Elem*> elemVec;
    ptGrid_->GetRegionIds( regionIds );

    for ( UInt iReg = 0; iReg < regionIds.GetSize(); iReg++ ) {

      // get region name and elements
      ptGrid_->GetElems(elemVec,regionIds[iReg]);
      std::string regionName = ptGrid_->RegionIdToName(regionIds[iReg]);
    
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
    
    // remember number of 'normal' elements
    UInt numelem =ptGrid_->GetNumElems();
    
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
        UInt nodeNumber = nodeNumbers[iNode];
        GiD_WriteElement( numelem++, (int*)(&nodeNumber) ); 
      }
      GiD_EndElements();
      GiD_EndMesh();
    }
    
  }
  
  void WriteResultsGiD::
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
    GiD_WriteElement( ptEl->elemNum, (int*)(connectM.GetPointer()) ); 
  }
  
  
  void WriteResultsGiD::
  WriteNodeSolutionTransient( const NodeStoreSol<Double> &sol, 
                              const UInt step, const Double time ) {
    ENTER_FCN( "WriteResultsGiD::WriteNodeSolutionTransient" );
    
    Vector<Double> gSol;
    StdVector<SolutionType> solTypes;
    UInt dof = 0;
    std::string title;
    sol.GetSolutionTypes(solTypes);
    
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
      sol.GetGlobalSolVector(solTypes[iSol], gSol );
      Enum2String( solTypes[iSol], title );
      dof = sol.GetDof(solTypes[iSol]);
      
      WriteNodeElemDataTrans( gSol, dof, title, GiD_OnNodes, time );
    }
  }
  
  
  void WriteResultsGiD::
  WriteElemSolutionTransient( const ElemStoreSol<Double> &sol, 
                              const UInt step, const Double time ) {

    ENTER_FCN ( "WriteResultsGiD::WriteElemSolutionTransient" );
    
    
    Vector<Double> gSol;
    StdVector<SolutionType> solTypes;
    UInt dof = 0;
    std::string title;
    sol.GetSolutionTypes(solTypes);
    
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
      sol.GetGlobalSolVector(solTypes[iSol], gSol );
      Enum2String( solTypes[iSol], title );
      dof = sol.GetDof(solTypes[iSol]);
      
      WriteNodeElemDataTrans( gSol, dof, title, GiD_OnGaussPoints, time );
    }
  }
  
  void WriteResultsGiD::
  WriteNodeSolutionHarmonic(const NodeStoreSol<Complex> & sol, 
                            const UInt step,
                            const Double frequency, 
                            const ComplexFormat format) {
    
    ENTER_FCN( "WriteResultsGiD::WriteNodeSolutionHarmonic" );
    
    Vector<Complex> gSol;
    StdVector<SolutionType> solTypes;
    UInt dof = 0;
    std::string title;
    sol.GetSolutionTypes(solTypes);
    
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
      sol.GetGlobalSolVector(solTypes[iSol], gSol );
      Enum2String( solTypes[iSol], title );
      dof = sol.GetDof(solTypes[iSol]);
      
      WriteNodeElemDataHarm( gSol, dof, title, GiD_OnNodes, 
                             frequency, format );
    }
    
  }
  
  void WriteResultsGiD::
  WriteElemSolutionHarmonic(const ElemStoreSol<Complex> & sol, 
                            const UInt step,
                            const Double frequency,
                            const ComplexFormat format) {
    
    ENTER_FCN( "WriteResultsGiD::WriteElemSolutionHarmonic" );
 
    Vector<Complex> gSol;
    StdVector<SolutionType> solTypes;
    UInt dof = 0;
    std::string title;
    sol.GetSolutionTypes(solTypes);
    
    for (UInt iSol=0; iSol<solTypes.GetSize(); iSol++) {
      sol.GetGlobalSolVector(solTypes[iSol], gSol );
      Enum2String( solTypes[iSol], title );
      dof = sol.GetDof(solTypes[iSol]);
      WriteNodeElemDataHarm( gSol, dof, title, GiD_OnGaussPoints, 
                             frequency, format );
    }
    
  }
  
  void WriteResultsGiD::
  WriteNodeElemDataTrans( const Vector<Double> & var, const UInt dof,
                          const std::string name, 
                          const GiD_ResultLocation entType, 
                          const Double time ) {
    ENTER_FCN ( "WriteResultsGiD::WriteNodeElemDataTrans" );
    
    // get number of entities
    UInt numEnt = 0;
    char * dummy;
    if ( entType == GiD_OnNodes ) {
      numEnt = ptGrid_->GetNumNodes();
      dummy = NULL;
    } else {
      numEnt = ptGrid_->GetNumElems();
      if ( dim_ == GiD_2D) {
        dummy = "mid-2D"; 
      } else {
        dummy = "mid-3D";
      }
    }

    // write Result header
    if ( dof == 1 ) {
      GiD_BeginResultHeader( name.c_str(), "transient", time,
                             GiD_Scalar, entType, dummy );
      GiD_ResultValues();
      for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
        GiD_WriteScalar( iEnt, var[iEnt-1]); 
      }
    } else {
      // Special treatment for mechanic stresses:
      // As up to now we can not determine, if the result is of stress type
      // we have to determine this by hand from the xml-file
      std::string stressName, subType;
      GiD_ResultType type = GiD_Vector;
      Enum2String( MECH_STRESS, stressName );
      if ( name == stressName ) {
        GiD_BeginResultHeader( name.c_str(), "transient", time,
                               GiD_Matrix, entType, dummy );
        // Check which mechanic subtype we have to treat
        params->Get( "subType", subType, "mechanic" );
        
        if( subType == "3d" ) {
          const char *names[] = {"xx", "xy", "yy", "xz", "yz", "zz"};
          GiD_ResultComponents( 6, names );
        } else if( subType == "planeStrain" ) {
          const char *names[] = {"_", "__", "xx", "___", "xy", "yy"};
          GiD_ResultComponents( 6, names );
        } else if( subType == "axi" ) {
          const char *names[] = {"phiphi", "_", "rr", "__", "rz", "zz"};
          GiD_ResultComponents( 6, names );
        } else {
          Error (" subType not known", __FILE__, __LINE__ );
        }
      } else {
        GiD_BeginResultHeader( name.c_str(), "transient", time,
                               type, entType, dummy );
        const char *names[] = {"x", "y", "z"};
        GiD_ResultComponents( 3, names );
      }
        
      GiD_ResultValues();
     
      // write results
      UInt offset = 0;
      for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
        offset = (iEnt-1) * dof;
        if ( dof == 6 ) {
          GiD_Write3DMatrix( iEnt, var[offset], var[offset+1], var[offset+2],
                             var[offset+3], var[offset+4],var[offset+5]  );
        } else if ( dof == 3 ) {
          GiD_WriteVector( iEnt, var[offset], var[offset+1], var[offset+2] );
        } else { 
          GiD_WriteVector( iEnt, var[offset], var[offset+1], 0.0);
        }
      }
    }
        
    GiD_EndResult();
    
  }     
  
  void WriteResultsGiD::
  WriteNodeElemDataHarm( const Vector<Complex> & var, const UInt dof,
                         const std::string name, 
                         const GiD_ResultLocation entType, 
                         const Double freq, 
                         const ComplexFormat outputFormat ) {
    ENTER_FCN ( "WriteResultsGiD::WriteNodeElemDataHarm" );
    
    // get number of entities
    UInt numEnt = 0;
    char * dummy;
    if ( entType == GiD_OnNodes ) {
      numEnt = ptGrid_->GetNumNodes();
      dummy = NULL;
    } else {
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
    
    // === Scalar entries ===
    if ( dof == 1 ) {
      
      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), "harmonic", freq,
                             GiD_Scalar, entType, dummy );
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
                             GiD_Scalar, entType, dummy );
      GiD_ResultValues();
      
      if (outputFormat == REAL_IMAG) {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          GiD_WriteScalar( iEnt, var[iEnt-1].imag()); 
        }
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          if ( std::abs(var[iEnt-1].imag()) > 1e-16 ) {
            GiD_WriteScalar( iEnt, std::atan2(var[iEnt-1].imag(),var[iEnt-1].real())*180/PI);             
          } else {
            GiD_WriteScalar( iEnt, 0.0); 
          }
        }
      }
      GiD_EndResult();
      
    } else {
      // === Vectorial entries ===

      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), "harmonic", freq,
                             GiD_Vector, entType, dummy );
      const char *names[] = {"x", "y", "z"};
      GiD_ResultComponents( 3, names );
      GiD_ResultValues();

      UInt offset = 0;
      if (outputFormat == REAL_IMAG) {

        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * dof;
          if ( dof == 3 ) {
            GiD_WriteVector( iEnt, var[offset].real(), 
                             var[offset+1].real(), 
                             var[offset+2].real() );
          } else { 
            GiD_WriteVector( iEnt, var[offset].real(), 
                             var[offset+1].real(), 0.0);
          }
        }
        GiD_EndResult(); 
        
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * dof;
          if ( dof == 3 ) {
            GiD_WriteVector( iEnt, std::abs(var[offset]),
                             std::abs(var[offset+1]), 
                             std::abs(var[offset+2]) );
          } else { 
            GiD_WriteVector( iEnt, 
                             std::abs(var[offset]), 
                             std::abs(var[offset+1]), 
                             0.0 );
          }
        }
        GiD_EndResult(); 
      }

      // --- Second component ---
      GiD_BeginResultHeader( outName2.c_str(), "harmonic", freq,
                             GiD_Vector, entType, dummy );
      GiD_ResultValues();

      offset = 0;
      if (outputFormat == REAL_IMAG) {

        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * dof;
          if ( dof == 3 ) {
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
          offset = (iEnt-1) * dof;
          if ( dof == 3 ) {
            GiD_WriteVector( iEnt, 
                             (std::abs(var[offset].imag()) > 1e-16) ? 
                             std::atan2(var[offset].imag(),var[offset].real() )*180/PI : 0.0,
                             (std::abs(var[offset+1].imag()) > 1e-16) ? 
                             std::atan2(var[offset+1].imag(),var[offset+1].real() )*180/PI : 0.0,
                             (std::abs(var[offset+2].imag()) > 1e-16) ? 
                             std::atan2(var[offset+2].imag(),var[offset+2].real() )*180/PI : 0.0 );
          } else { 
            GiD_WriteVector( iEnt, 
                             (std::abs(var[offset].imag()) > 1e-16) ? 
                             std::atan2(var[offset].imag(),var[offset].real() )*180/PI : 0.0,
                             (std::abs(var[offset+1].imag()) > 1e-16) ? 
                             std::atan2(var[offset+1].imag(),var[offset+1].real() )*180/PI : 0.0,
                             0.0 );
          }
        }
        GiD_EndResult(); 
      }
    }
  }     
  
  // ********
  //   Init
  // ********
  void WriteResultsGiD::Init( Grid *aptgrid ) {
    
    ENTER_FCN( "WriteResultsGiD::OpenFile" );

    ptGrid_ = aptgrid;

    // initialize history files
    InitHistoryFiles();

    // determine dimension of grid
    UInt dim = ptGrid_->GetDim();
    if ( dim == 2 ) {
      dim_ = GiD_2D;
    } else if ( dim == 3 ) {
      dim_ = GiD_3D;
    } else {
      (*error) << "WriteResultsGiD: Grid dimension " << dim
               << " is not supported by GiD mesh format.";
      Error( __FILE__, __LINE__ );
    }

    // check, if only tetrahedra are present
    // -> only in this case we can treat tetrahedras as elements with 4 nodes.
    // In all other cases (hexa and tets mixed) tetrahedras are treated as
    // degenerated hexahedras.
    if ( dim_ == GiD_3D ) {
      if ( ptGrid_->GetNumElemOfType( TET ) > 0  &&
           ptGrid_->GetNumElemOfType( HEX ) ==  0  &&
           ptGrid_->GetNumElemOfType( PYR ) ==  0  &&
           ptGrid_->GetNumElemOfType( WED ) ==  0  ) {
        degen3DElems_ = false;
      }
    }
    
  }
} // end of namespace
