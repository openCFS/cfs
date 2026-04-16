// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <complex>

#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
using namespace boost::assign;
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "FeBasis/BaseFE.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "SimOutGiD.hh"

#define USE_CONST
// Header of postprocessing library
#pragma GCC diagnostic push // pop at end of file only
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <gidpost.h>

#undef USE_CONST

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simOutputGiD, "SimOutputGiD")

    SimOutputGiD::SimOutputGiD( const std::string& fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode, 
                                bool isRestart )
  : SimOutput( fileName, outputNode, infoNode, isRestart ),
    isInitialized_(false)
  {
    // Initialize variables
    formatName_ = "gid";
    fileName_ = fileName;
    std::string dirString = "results_" + formatName_; 
    outputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString; 
        
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    dim_ = GiD_2D;
    actMeshId_ = 1;
    lastStepVal_ = -1.0;
    lastStepRepeated_ = 0;
    isAscii_ = true;
    groupEigenFreqs_ = true;
    printGridOnly_ = false;
    gridWritten_ = false;

    // Determine, if binary result file should be written
    isAscii_ = !(myParam_->Get("binaryFormat")->As<bool>() );
    
    // Determine, if eigenfrequencies should be grouped
    if( myParam_->Has("groupEigenFreqs") ) {
      groupEigenFreqs_= myParam_->Get("groupEigenFreqs")->As<bool>();
    }
    
    //! Determine if sequenceSteps sould be merged
    if (myParam_->Has("mergeSequenceSteps")) {
      mergeSequenceSteps_= myParam_->Get("mergeSequenceSteps")->As<bool>();
    }

    // concatenate output file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    // store complete name including directory into fileName_
    fileName_ = fs::path(dirName_ / fileName_ ).string();
    if(isRestart) {
      WARN("Restart support for gidpost was removed with upgrade to gidpost 2.11 (around Feb. 2025). This is untested ...")
    }
  }


  SimOutputGiD::~SimOutputGiD() {
    if(isInitialized_) 
    {
      // Close result file
      GiD_ClosePostResultFile();
    }
  }


  void SimOutputGiD::WriteGrid() {

    LOG_TRACE(simOutputGiD) << "Writing mesh";

    // Leave, if grid was already written or restarted simulation
    if(gridWritten_ || isRestart_) 
      return;
    
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
      if ( ptGrid_->GetNumElemOfDim(1) > 0 ) {
        GiD_BeginGaussPoint( "mid-2D", GiD_Linear, NULL, 1, 0, 1);
        GiD_EndGaussPoint();
      }
    } else {
      GiD_BeginGaussPoint( "mid-3D", GiD_Hexahedra, NULL, 1, 1, 1);
      GiD_EndGaussPoint();
      if ( ptGrid_->GetNumElemOfDim(2) > 0 ) {
        GiD_BeginGaussPoint( "mid-3D", GiD_Quadrilateral, NULL, 1, 1, 1);
        GiD_EndGaussPoint();
      }
      if ( ptGrid_->GetNumElemOfDim(1) > 0 ) {
        GiD_BeginGaussPoint( "mid-3D", GiD_Linear, NULL, 1, 0, 1);
        GiD_EndGaussPoint();
      }
    }

    // close mesh file (only needed in ASCII case)
    if ( isAscii_ == true ) {
      GiD_ClosePostMeshFile();
    }

    // close result file, if binary format and only print grid
//    if( isAscii_ == false && printGridOnly_ == true ) {
//      GiD_ClosePostResultFile();
//    }

    gridWritten_ = true;
    LOG_TRACE(simOutputGiD)<< "Finished writing of mesh" << std::endl;

  }

  void SimOutputGiD::WriteNodes() {

    LOG_TRACE(simOutputGiD) << "Writing nodes";

    // get number of nodes
    UInt numNodes = ptGrid_->GetNumNodes();

    // write nodal coordinates
    GiD_BeginMesh( "Nodes", static_cast<GiD_Dimension>(dim_), GiD_Point, 1 );

    GiD_BeginCoordinates();
    if ( dim_ == GiD_2D ) {

      Vector<Double> point(2);
      for ( UInt i = 1; i <= numNodes; i++ ) {
        ptGrid_->GetNodeCoordinate(point,i);
        GiD_WriteCoordinates2D(i, point[0], point[1]);
      }
    } else {
      Vector<Double> point(3);
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

  void SimOutputGiD::WriteElementMesh( const std::string& name ,
                                       StdVector<Elem*> & elemVec ) {
    UInt numElemNodes;
    GiD_ElementType eType;
    numElemNodes = 0;
    eType = GiD_NoElement;

    // determine element type and number of nodes
    UInt dim = Elem::shapes[elemVec[0]->type].dim;
    if ( dim == 1 ) {

      eType = GiD_Linear;
      numElemNodes = elemVec[0]->connect.GetSize();

    } else if ( Elem::shapes[elemVec[0]->type].dim == 2 ) {

      eType = GiD_Quadrilateral;
      numElemNodes = elemVec[0]->connect.GetSize();

      switch(numElemNodes) {
      case 3:
      case 4:
        numElemNodes = 4;
        break;
      case 6:
      case 8:
        numElemNodes = 8;
        break;
      }

    } else if ( Elem::shapes[elemVec[0]->type].dim == 3 ) {

      eType = GiD_Hexahedra;
      numElemNodes = elemVec[0]->connect.GetSize();

      switch(numElemNodes) {
      case 4:
      case 5:
      case 6:
      case 8:
        numElemNodes = 8;
        break;
      case 9:
        numElemNodes = 9;
        break;
      case 10:
      case 13:
      case 14:
      case 15:
      case 18:
      case 20:
        numElemNodes = 20;
        break;
      case 27:
        numElemNodes = 27;
        break;
      }

    }

    GiD_BeginMesh( name.c_str(), static_cast<GiD_Dimension>(dim_), eType, numElemNodes );

    // write dummy coordinate section
    GiD_BeginCoordinates();
    GiD_EndCoordinates();

    // write all element declarations
    GiD_BeginElements();
    for ( UInt iElem = 0; iElem < elemVec.GetSize(); iElem++ ) {
      Elem * ptEl = elemVec[iElem];

      WriteElement( ptEl, numElemNodes);
    }
    GiD_EndElements();

    GiD_EndMesh();
  }

  void SimOutputGiD::WriteElements() {
 
     LOG_TRACE(simOutputGiD) << "Writing elements";

    // iterate over all regions
    StdVector<RegionIdType> regionIds;
    StdVector<Elem*> elemVec;

    // =================
    //  REGION ELEMENTS
    // =================
    for ( UInt iReg = 0; iReg < ptGrid_->GetNumRegions(); iReg++ ) {
      // get region name and elements
      std::string regionName = ptGrid_->GetRegion().ToString(iReg);
      ptGrid_->GetElems(elemVec,iReg);

      // write list of elements
      WriteElementMesh( regionName, elemVec );

      // increment current mesh id
      actMeshId_++;
    }

    LOG_TRACE(simOutputGiD) << "Writing named elements/nodes";

    // ================
    //  ELEMENT GROUPS
    // ================
    StdVector<std::string> elemNames;
    StdVector<UInt> elemNumbers;
    ptGrid_->GetListElemNames(elemNames);

    for( UInt i = 0; i < elemNames.GetSize(); i++ ) {

      ptGrid_->GetElemsByName( elemVec, elemNames[i] );

      // write list of elements
      WriteElementMesh( elemNames[i], elemVec );

      // increment current mesh id
      actMeshId_++;
    }

    // =============
    //  NODE GROUPS
    // =============
    StdVector<std::string> nodeNames;
    StdVector<UInt> nodeNumbers;

    // write named nodes to mesh file
    // Note:  Named nodes are written as 'point' elements, so we need
    // for each set of named nodes a separate MESH section


    // remember number of 'normal' elements
    UInt numelem = ptGrid_->GetNumElems() + 1;

    ptGrid_->GetListNodeNames(nodeNames);

    // check if there are any named nodes / elems in the grid
    if( nodeNames.GetSize() == 0 ) {
      return;
    }

    for( UInt i = 0; i < nodeNames.GetSize(); i++ ) {

      // get named nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);
      GiD_BeginMesh( nodeNames[i].c_str(), static_cast<GiD_Dimension>(dim_), GiD_Point, 1 );

      // write empty coordinate section
      GiD_BeginCoordinates();
      GiD_EndCoordinates();

      // write 'nodal' elements
      GiD_BeginElements();

      for ( UInt iNode = 0; iNode < nodeNumbers.GetSize(); iNode ++ ) {
        UInt nodeNumber[2] = { nodeNumbers[iNode], UInt(actMeshId_) };
        GiD_WriteElementMat( numelem++, (int*)(&nodeNumber) );
      }
      GiD_EndElements();
      GiD_EndMesh();

      // increment current mesh id
      actMeshId_++;
    }

    LOG_TRACE(simOutputGiD) << "Finished writing elements" << std::endl;
  }

  void SimOutputGiD::
  WriteElement( Elem * ptEl, UInt numNodes) {

    UInt elemNum;
    Elem::FEType eType;

    eType = ptEl->type;
    elemNum = ptEl->elemNum;
    StdVector<UInt> const & elConn = ptEl->connect;
    std::vector<UInt> reorderIdx;

    // We 
    switch(eType)
    {
    default:
      reorderIdx.resize(numNodes);
      for(UInt i=0; i<numNodes; i++) reorderIdx[i] = i;
      break;

    case Elem::ET_TRIA3:
    case Elem::ET_TRIA6:
      reorderIdx += 0,1,2,2;
      if(eType != Elem::ET_TRIA3) {
        reorderIdx += 3,4,5,5;
      }
      break;

    case Elem::ET_TET4:
    case Elem::ET_TET10:
      reorderIdx += 0,1,2,2,3,3,3,3;
      if(eType != Elem::ET_TET4) {
        reorderIdx += 4,5,6,6,7,8,9,9,3,3,3,3;
      }
      break;

    case Elem::ET_HEXA8:
      reorderIdx += 4,5,6,7,0,1,2,3;
      break;

      // NOTE: The numberings of hexehadras in gid differs for
      // the quadratic case!
    case Elem::ET_HEXA20:
    case Elem::ET_HEXA27:
      reorderIdx += 0,1,2,3,4,5,6,7,8,9,10,11,16,17,18,19,12,13,14,15;
      if(eType == Elem::ET_HEXA27) {
        reorderIdx += 24,20,21,22,23,25,26;
      }
      break;
      
    case Elem::ET_PYRA5:
    case Elem::ET_PYRA13:
    case Elem::ET_PYRA14:
      reorderIdx += 0,1,2,3,4,4,4,4;
      if(eType != Elem::ET_PYRA5) {
        reorderIdx += 5,6,7,8,9,10,11,12,4,4,4,4;
      }
      break;

    case Elem::ET_WEDGE6:
    case Elem::ET_WEDGE15:
    case Elem::ET_WEDGE18:
      reorderIdx += 0,1,2,2,3,4,5,5;
      if(eType != Elem::ET_WEDGE6) {
        reorderIdx += 6,7,8,8,12,13,14,14,9,10,11,11;
      }
      break;
    }

    std::vector<int> connect(numNodes+1);
    for(UInt i=0; i<numNodes; i++) {
      connect[i] = elConn[reorderIdx[i]];
    }

    connect[numNodes] = actMeshId_;
    GiD_WriteElementMat( elemNum, &connect[0] );
  }


  void SimOutputGiD::BeginMultiSequenceStep( UInt step,
                                             BasePDE::AnalysisType type,
                                             UInt numSteps ) {
    // Write grid at the beginning of a multi sequence step.
    WriteGrid();

    actAnalysis_ = type;
    actMSStep_ = step;
  }

  void SimOutputGiD::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory ) {

    ResultInfo & actDof = *(sol->GetResultInfo());
    LOG_DBG(simOutputGiD) << "Registering output '" << actDof.resultName
                          << "' with saveBegin " << saveBegin
                          << ", saveInc " << saveInc
                          << ", saveEnd " << saveEnd;

  }

  void SimOutputGiD::BeginStep( UInt stepNum, Double stepVal ) {

    lastStepVal_ = actStepVal_;
    
    actStep_ = stepNum;
    actStepVal_ = stepVal;

    // Check, if current step value is the same as the previous step value.
    // This can happen e.g. in an eigenfrequency analysis, where there maybe
    // more than one mode shape associated with an frequency. We count the 
    // number of occurences, to include it in the the resultname.
    
    if( actAnalysis_ == BasePDE::EIGENFREQUENCY ) {
      if( std::abs(actStepVal_) > EPS ) {
        if( std::abs(actStepVal_ - lastStepVal_) / actStepVal_ < 1e-4) {
          lastStepRepeated_++;
        } else {
          lastStepRepeated_ = 0;
        }
      }
    } else {
      lastStepRepeated_ = 0;
    }
    resultMap_.clear();
  }

  void SimOutputGiD::AddResult( shared_ptr<BaseResult> sol ) {
    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGiD) << "Adding result '"
                          << actDof.resultName  << "'";

    resultMap_[sol->GetResultInfo()->resultName].Push_back( sol );
  }

  void SimOutputGiD::FinishStep( ) {

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
      if(!ValidateNodesAndElements(actInfo)) continue;
      LOG_DBG(simOutputGiD) << "Writing result '" << title << "'";

      if( actResults[0]->GetEntryType() == BaseMatrix::DOUBLE ) {
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
  FinishMultiSequenceStep( ) {
  }


  void SimOutputGiD::
  WriteNodeElemDataTrans( const Vector<Double> & var,
                          const StdVector<std::string> & dofNames,
                          const std::string& name,
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {

    // assemble name for analysis step
    std::string analysisName = "transient";
    if( !mergeSequenceSteps_)
      analysisName += "_" + std::to_string( actMSStep_ );


    // get number of entities
    UInt numEnt = 0;
    const char * dummy;
    GiD_ResultLocation loc;
    if ( entityType == ResultInfo::NODE ) {
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

    // In case of an eigenfrequency analysis we can optionally group the eigenfrequencies
    std::string outName;
    if (actAnalysis_ == BasePDE::EIGENFREQUENCY  &&
        groupEigenFreqs_ == true )  {
      std::string temp;
      temp = "ef" + std::to_string(lastStepRepeated_+1);
      temp += "//";
      temp += name;
      outName = temp;
    } else {
      outName = name;
    }
    
    // write Result header
    if ( entryType == ResultInfo::SCALAR )  {
      GiD_BeginResultHeader( outName.c_str(), analysisName.c_str(), time,
                             GiD_Scalar, loc, dummy );
      GiD_ResultValues();
      for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
        GiD_WriteScalar( iEnt, var[iEnt-1]);
      }
    } else if( entryType == ResultInfo::VECTOR ) {
      GiD_BeginResultHeader( outName.c_str(), analysisName.c_str(), time,
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
    GiD_BeginResultHeader( outName.c_str(), analysisName.c_str(), time,
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

    // assemble name for analysis step
    std::string analysisName = "harmonic";
    if( !mergeSequenceSteps_)
      analysisName += "_" + std::to_string( actMSStep_ );

   // get number of entities
    UInt numEnt = 0;
    const char * dummy;
    GiD_ResultLocation loc;
    if ( entityType == ResultInfo::NODE ) {
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
    
    // In case of an eigenfrequency analysis we can optionally group the eigenfrequencies
    std::string outBase;
    if (actAnalysis_ == BasePDE::EIGENFREQUENCY  &&
        groupEigenFreqs_ == true )  {
      std::string temp;
      temp = "ef" + std::to_string(lastStepRepeated_+1);
      temp += "//";
      temp += name;
      outBase = temp;
    } else {
      outBase = name;
    }

    // determine complete name of output quantities
    std::string outName1, outName2;
    if (outputFormat == REAL_IMAG) {
      outName1 = outBase + "-real";
      outName2 = outBase + "-imag";
    } else {
      outName1 = outBase + "-amp";
      outName2 = outBase + "-phase";
    }

    UInt numDofs = dofNames.GetSize();
    


    // === Scalar entries ===
    if ( entryType == ResultInfo::SCALAR ) {

      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), analysisName.c_str(), freq,
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
      GiD_BeginResultHeader( outName2.c_str(), analysisName.c_str(), freq,
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
      GiD_BeginResultHeader( outName1.c_str(), analysisName.c_str(), freq,
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
      GiD_BeginResultHeader( outName2.c_str(), analysisName.c_str(), freq,
                             GiD_Vector, loc, dummy );
      GiD_ResultComponents( 3, names );
      GiD_ResultValues();

      offset = 0;
      if (outputFormat == REAL_IMAG) {

        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          switch ( numDofs ) {
            case 3:
              GiD_WriteVector( iEnt, var[offset].imag(), 
                               var[offset+1].imag(), 
                               var[offset+2].imag() );
              break;
            case 2:
              GiD_WriteVector( iEnt, var[offset].imag(), 
                               var[offset+1].imag(), 
                               0.0 );
              break;
            case 1:
              GiD_WriteVector( iEnt, var[offset].imag(), 0.0, 0.0 );
              break;
          }
        }
        GiD_EndResult();
      } else {
        for ( UInt iEnt = 1; iEnt <= numEnt; iEnt++ ) {
          offset = (iEnt-1) * numDofs;
          switch ( numDofs ) {
            case 3:
              GiD_WriteVector( iEnt, CPhase( var[offset] ), 
                               CPhase( var[offset+1] ), 
                               CPhase( var[offset+2] ) );
              break;
            case 2:
              GiD_WriteVector( iEnt, CPhase( var[offset] ),
                               CPhase( var[offset+1]),
                               0.0 );
              break;
            case 1:
              GiD_WriteVector( iEnt, CPhase( var[offset] ), 0.0, 0.0 );
              break;
          }
        }
        GiD_EndResult();
      }
    } else if ( entryType == ResultInfo::TENSOR ) {
      // === Tensor entries ===

      // --- First component ---
      GiD_BeginResultHeader( outName1.c_str(), analysisName.c_str(), freq,
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
      GiD_BeginResultHeader( outName2.c_str(), analysisName.c_str(), freq,
                             GiD_Matrix, loc, dummy );
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
  void SimOutputGiD::Init( Grid* ptGrid, bool printGridOnly) {


    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;

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

    // All meshes will be treated as degenerated quad and hexa meshes!

    std::string postFileName;
    
    
    // Open result file
    if ( isAscii_ == true) {
      postFileName = fileName_ + ".post.res";
      isInitialized_ = !GiD_OpenPostResultFile( postFileName.c_str(),
                                                GiD_PostAscii );
    } else {
      postFileName = fileName_ + ".post.bin";
      isInitialized_ = !GiD_OpenPostResultFile( postFileName.c_str(),
                                                GiD_PostBinary );
    }

    // print grid
    if(printGridOnly) {
      WriteGrid();

      if(isAscii_) 
      {      
        try 
        {
          fs::remove( postFileName );
        } catch (std::exception &ex) {
          EXCEPTION(ex.what());
        }
      }
    }
  }


} // end of namespace

#pragma GCC diagnostic pop

