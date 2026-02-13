// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <ctime>
#include <fstream>

#include "DataInOut/SimInput.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#include "GmshHelper.hh"
#include "SimOutputGmsh.hh"

namespace CoupledField {

  // declare logging stream
  DEFINE_LOG(simOutputGmsh, "SimOutputGmsh")

   //*****************
  //   Constructor
  //*****************
  SimOutputGmsh::SimOutputGmsh( const std::string fileName,
                                PtrParamNode outputNode,
                                PtrParamNode infoNode, bool isRestart ) :
  SimOutput( fileName, outputNode, infoNode, isRestart ),
    ascii_(true),
    bigEndian_(false),
    stepNumOffset_(0),
    numRegions_(0)
  {    
    
    // The restart case is currently not implemented, i.e. resuls from a 
    // partial simulation get overwritten.
    if( isRestart_ ) {
      WARN( "The GMSH-Writer is currently not adapted to write restarted "
          "results correctly, thus the results of the previous run get"
          " overwritten." );
    }
    
    // Initialize variables
    formatName_ = "gmsh";
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    std::ostringstream strBuffer;

    std::string dirString = "results_" + formatName_; 
    outputNode->GetValue("directory", dirString, ParamNode::PASS );
    dirName_ = dirString;  
    
    fileName_ = fileName;

    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }
    
    // Change defaults according to XML file
    if(myParam_->Get("binaryFormat", ParamNode::PASS)->As<std::string>() == "yes")
    {
        ascii_ = false;
    }

    std::string endianness = "native";
    endianness = myParam_->Get("endianness", ParamNode::PASS)->As<std::string>();

    if(!ascii_) 
    {
      if(endianness == "little")
      {
        bigEndian_ = false;
      } else if (endianness == "big")
      {
        bigEndian_ = true;
      } else if (endianness == "native")
      {
        switch(HOST_ENDIAN_ORDER) 
        {  
        case LITTLE_ENDIAN_ORDER:
          bigEndian_ = false;
          break;
        case BIG_ENDIAN_ORDER:
          bigEndian_ = true;
          break;
        }        
      }
    }

    if(bigEndian_ && (HOST_ENDIAN_ORDER != BIG_ENDIAN_ORDER))
    {
      uiSwap_.swap = true;
      iSwap_.swap = true;
      dSwap_.swap = true;
    }
    else
    {
      uiSwap_.swap = false;
      iSwap_.swap = false;
      dSwap_.swap = false;
    }
    
    GmshHelper::InitElemNodeMap();

  }


  // **********************
  //   Default Destructor
  // **********************
  SimOutputGmsh::~SimOutputGmsh() {
    if(output_)
      output_->close();
    
    delete output_;
  }
  
  void SimOutputGmsh::Init( Grid* ptGrid, bool printGridOnly ) {
    
    ptGrid_ = ptGrid;

    if(printGridOnly)
      WriteGrid( printGridOnly );

  }


  void SimOutputGmsh::BeginMultiSequenceStep( UInt step,
                                             BasePDE::AnalysisType type,
                                             UInt numSteps )
  {
    currAnalysis_ = type;
    currMsStep_ = step;
    if ( step == 1 ) { 
      // first step must write the grid
      WriteGrid( false );
    }
    else if ( step == 2 ) { 
      /* if this is the second sequence step, one must not write the grid. 
         we leave this `append` file handle open until the end of the simu
         hence we need to open it only at the second sequence step. */
      // close old write handle
      if(output_) {
        output_->close();
        delete output_;
        output_ = NULL;
      }
      std::string gridFileName = fileName_ + ".msh";
      // Generate basename for output file
      fs::path filePath = dirName_ / gridFileName;
      
      if (ascii_) {
        output_ = new std::ofstream(filePath, std::ios_base::app);
      }
      else {
        output_ = new std::ofstream(filePath, std::ofstream::binary | std::ios_base::app);
      }
      if ( !output_ ) {
        EXCEPTION("Could not open file " << gridFileName
                  << " for appending Gmsh output");
      }
    } else {
      // relax, append file handle already opened
    }
  }
  
  void SimOutputGmsh::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory )
  {

    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGmsh) << "Registering output '" << actDof.resultName 
                          << "' with saveBegin " << saveBegin
                          << ", saveInc " << saveInc 
                          << ", saveEnd " << saveEnd;
  }


  //! Begin single analysis step
  void SimOutputGmsh::BeginStep( UInt stepNum, Double stepVal )
  {

    resultMap_.clear();
    
    actStep_ = stepNum;
    actStepVal_ = stepVal;
    if( currAnalysis_ == BasePDE::TRANSIENT ||
        currAnalysis_ == BasePDE::STATIC  ) { 
      actStep_ += stepNumOffset_;
    }

  }


  //! Add result to current step
  void SimOutputGmsh::AddResult( shared_ptr<BaseResult> sol )
  {
    ResultInfo & actDof = *(sol->GetResultInfo());
    
    LOG_DBG(simOutputGmsh) << "Adding result '" 
                           << actDof.resultName << "'";
      
    resultMap_[sol->GetResultInfo()->resultName].Push_back(sol);
  }


  //! End single analysis step
  void SimOutputGmsh::FinishStep( )
  {

    LOG_TRACE(simOutputGmsh) << "Starting to finish Step";
    
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
      LOG_DBG(simOutputGmsh) << "Writing result '" << title << "'";

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

  void SimOutputGmsh::
  WriteNodeElemDataTrans( const Vector<Double> & var,
                          const StdVector<std::string> & dofNames,
                          const std::string& name,
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {
    // assemble name for analysis step
    std::string analysisName = "transient";
    analysisName += "_" + lexical_cast<std::string>( currMsStep_ );

    std::string resultTag;

    // get number of entities
    UInt numEnt = var.GetSize() / dofNames.GetSize();
    if ( entityType == ResultInfo::NODE ) {
      resultTag = "NodeData";
    } else {
      resultTag = "ElementData";
    }

    // Write begin result tag
    (*output_) << "$" << resultTag << std::endl;

    // Write one string tag for result name
    (*output_) << 1 << std::endl;
    (*output_) << "\"" << name << "\"" << std::endl;

    // Write one real tag for time value
    (*output_) << 1 << std::endl;
    (*output_) << actStepVal_ << std::endl;

    // Write three integer tags for time step, number of dofs
    // and number of nodes/elems
    (*output_) << 3 << std::endl;
    (*output_) << (actStep_-1) << std::endl;
    UInt numDofs = dofNames.GetSize();
    (*output_) << numDofs << std::endl;
    (*output_) << numEnt << std::endl;

    for ( UInt iEnt = 0; iEnt < numEnt; iEnt++ ) {
      UInt offset = iEnt * numDofs;

      if(!ascii_) {
        UInt entNum = iEnt+1;
        entNum = uiSwap_.EndianSwapBytes(entNum);
        (*output_).write((char*) &entNum, sizeof(entNum));
      } else {
        (*output_) << (iEnt+1);
      }
      
      for ( UInt dof = 0; dof < numDofs; dof++ ) {
        if(!ascii_) {
          Double d = var[offset + dof];
          d = dSwap_.EndianSwapBytes(d);
          (*output_).write((char*) &d, sizeof(d));
        } else {
          (*output_) << " " << var[offset + dof];
          (*output_) << std::endl;        
        }
      }
    }
    
    (*output_) << "$End" << resultTag << std::endl;
  }

  void SimOutputGmsh::
  WriteNodeElemDataHarm( const Vector<Complex> & var,
                         const StdVector<std::string> & dofNames,
                         const std::string name,
                         const ResultInfo::EntryType entryType,
                         ResultInfo::EntityUnknownType entityType,
                         const Double freq,
                         const ComplexFormat outputFormat ) {
    UInt n = var.GetSize();
    Vector<Double> var1(n);
    Vector<Double> var2(n);
    
    // determine complete name of output quantities
    std::string outName1, outName2;
    if (outputFormat == REAL_IMAG) {
      outName1 = name + "-real";
      outName2 = name + "-imag";
    } else {
      outName1 = name + "-amp";
      outName2 = name + "-phase";
    }
    
    for ( UInt i = 0; i < n; i++ ) {
      if (outputFormat == REAL_IMAG) {
        var1[i] = var[i].real();
        var2[i] = var[i].imag();
      } else {
        var1[i] = std::abs(var[i]);
        var2[i] = CPhase(var[i]);
      }      
    }
    
    WriteNodeElemDataTrans( var1, dofNames, outName1,
                            entryType, entityType, freq );
    WriteNodeElemDataTrans( var2, dofNames, outName2,
                            entryType, entityType, freq );
  }
  
  //! End multisequence step
  void SimOutputGmsh::FinishMultiSequenceStep( ) {
    // set offset for step value and number to last values
    if( currAnalysis_ == BasePDE::TRANSIENT ||
          currAnalysis_ == BasePDE::STATIC ) {
      stepNumOffset_ = actStep_;
    }

  }


  //! Finalize the output
  void SimOutputGmsh::Finalize() {
    if(output_) {
      output_->close();
      delete output_;
      output_ = NULL;
    }
  }


  // **************
  //   WriteNodes
  // **************
  void SimOutputGmsh::WriteNodes( std::ofstream * gridFile ) {

    //get and write number of nodes on the level
    UInt numnodes = ptGrid_->GetNumNodes();
    (*output_) << "$Nodes" << std::endl;
    (*output_) << numnodes << std::endl;

    std::vector<GmshHelper::NodeStruct> nodes;
    if (!ascii_) {
      nodes.resize(numnodes);
    }

    // write x,y,z-coordinate
    for (UInt i = 1; i <= numnodes; i++) {
      Vector<Double> p;
      
      ptGrid_->GetNodeCoordinate3D(p, i);  

      if (ascii_) {
        (*output_) << i << " " << p[0] << " " << p[1] << " "
                   << p[2] << std::endl;
      }
       else {
         nodes[i-1].nodeId = uiSwap_.EndianSwapBytes(i);
         nodes[i-1].x = dSwap_.EndianSwapBytes(p[0]);
         nodes[i-1].y = dSwap_.EndianSwapBytes(p[1]);
         nodes[i-1].z = dSwap_.EndianSwapBytes(p[2]);

         (*output_).write((char*) &nodes[i-1].nodeId, sizeof(nodes[i-1].nodeId));
         (*output_).write((char*) &nodes[i-1].x, sizeof(nodes[i-1].x));
         (*output_).write((char*) &nodes[i-1].y, sizeof(nodes[i-1].y));
         (*output_).write((char*) &nodes[i-1].z, sizeof(nodes[i-1].z));

       }
    }

    if (!ascii_) {
      //  (*output_) << std::endl;
    }

    (*output_) << "$EndNodes" << std::endl;
  }


  // **************
  //   WriteCells
  // **************
  void SimOutputGmsh::WriteCells( std::ofstream * gridFile ) {
    // read information about number of elements 
    UInt numElems = ptGrid_->GetNumElems();
    UInt namedElemsOffset = numElems;
    numRegions_ = ptGrid_->GetNumRegions();
    
    // ================
    //  ELEMENT GROUPS
    // ================
    StdVector<std::string> elemNames;
    StdVector<Elem*> elemVec;
    ptGrid_->GetListElemNames(elemNames);

    for( UInt i = 0, n = elemNames.GetSize(); i < n; i++ ) {
      ptGrid_->GetElemsByName( elemVec, elemNames[i] );
      numElems += elemVec.GetSize();
    }

    // =============
    //  NODE GROUPS
    // =============
    StdVector<std::string> nodeNames;
    StdVector<UInt> nodeNumbers;

    ptGrid_->GetListNodeNames(nodeNames);

    for( UInt i = 0, n = nodeNames.GetSize(); i < n; i++ ) {
      // get named nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);
      numElems += nodeNumbers.GetSize();
    }

    
    
    (*output_) << "$Elements" << std::endl;
    (*output_) << numElems << std::endl;

    UInt numTags = 3;
    //    UInt elemRecSize = 1 + 3 + ptGrid_->GetMaxNumNodesPerElem();
    //    std::vector<UInt> elemRecords(numElems * elemRecSize );
    //    std::fill(elemRecords.begin(), elemRecords.end(), 0);

    std::vector<UInt> elemRecord(1 + numTags + ptGrid_->GetMaxNumNodesPerElem());

    std::map< Elem::FEType, std::vector<UInt> > elemRecords;

    // write elements
    for (UInt i = 0, n=ptGrid_->GetNumElems(); i < n; i++) {
      Elem::FEType feType;
      RegionIdType region;

      // Assign element number
      elemRecord[0] = i+1;
      // Retrieve type, region and connectivity
      ptGrid_->GetElemData(i+1, feType, region, &elemRecord[1+numTags]);
      // Assign region index
      // physical entity = elementary entity = partition = region + 1
      elemRecord[1] = elemRecord[2] = elemRecord[3] = region + 1;

      // Copy first part of element record to records vector for given type
      std::copy(elemRecord.begin(), elemRecord.begin() + 1 + numTags,
                std::back_inserter(elemRecords[feType]));
      
      typedef GmshHelper::ElemNodeMapType::right_map::const_iterator right_const_iterator;
      right_const_iterator right_iter = GmshHelper::elemNodeMap_[feType].right.begin();
      right_const_iterator iend = GmshHelper::elemNodeMap_[feType].right.end();

      std::vector<UInt> gmshConn(GmshHelper::elemNodeMap_[feType].size());
      
      // Iterate over element nodes and reorder connectivity for Gmsh
      for( ; right_iter != iend; ++right_iter)
      {
        gmshConn[right_iter->second] = elemRecord[1 + numTags + right_iter->first];
        
        // elemRecords[feType].push_back(elemRecord[1 + numTags + right_iter->second]);

        //        std::cout << elemRecord[1 + numTags + right_iter->second] << std::endl;
        
        // std::cout << right_iter->first << " --> " << right_iter->second << std::endl;
        // std::cout << (right_iter->first+1) << " <-- " << elemRecord[1 + numTags + right_iter->second] << std::endl;
      }

      std::copy(gmshConn.begin(), gmshConn.end(),
                std::back_inserter(elemRecords[feType]));

//      std::ostream_iterator< UInt > output( std::cout, " " );
//
//      std::copy(gmshConn.begin(), gmshConn.end(),
//                output);
//      std::cout << std::endl;

    }

    LOG_TRACE(simOutputGmsh) << "Writing named elements/nodes";

    // Write named elements
    UInt k = namedElemsOffset + 1;
    for( UInt i = 0, n = elemNames.GetSize(); i < n; i++ ) {
      ptGrid_->GetElemsByName( elemVec, elemNames[i] );
      numRegions_++;
      
      for( UInt j = 0, m = elemVec.GetSize(); j < m; j++ ) {
        Elem::FEType feType;
        RegionIdType region;

        // Assign element number
        elemRecord[0] = k++;
        // Retrieve type, region and connectivity
        ptGrid_->GetElemData(elemVec[j]->elemNum, feType, region, &elemRecord[1+numTags]);
        // Assign region index
        elemRecord[1] = elemRecord[2] = elemRecord[3] = numRegions_;
        
        // Copy first part of element record to records vector for given type
        std::copy(elemRecord.begin(), elemRecord.begin() + 1 + numTags,
                  std::back_inserter(elemRecords[feType]));
      
        typedef GmshHelper::ElemNodeMapType::right_map::const_iterator right_const_iterator;
        right_const_iterator right_iter = GmshHelper::elemNodeMap_[feType].right.begin();
        right_const_iterator iend = GmshHelper::elemNodeMap_[feType].right.end();

        std::vector<UInt> gmshConn(GmshHelper::elemNodeMap_[feType].size());
      
        // Iterate over element nodes and reorder connectivity for Gmsh
        for( ; right_iter != iend; ++right_iter)
        {
          gmshConn[right_iter->second] = elemRecord[1 + numTags + right_iter->first];
        }

        std::copy(gmshConn.begin(), gmshConn.end(),
                  std::back_inserter(elemRecords[feType]));
      }
    }

    for( UInt i = 0, n = nodeNames.GetSize(); i < n; i++ ) {
      // get named nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);

      numRegions_++;

      for( UInt j = 0, m = nodeNumbers.GetSize(); j < m; j++ ) {
        Elem::FEType feType = Elem::ET_POINT;

        // Assign element number
        elemRecord[0] = k++;
        elemRecord[1] = numRegions_;
        elemRecord[2] = numRegions_;
        elemRecord[3] = numRegions_;

        // Copy first part of element record to records vector for given type
        std::copy(elemRecord.begin(), elemRecord.begin() + 1 + numTags,
                  std::back_inserter(elemRecords[feType]));
      
        std::copy(&nodeNumbers[j], (&nodeNumbers[j])+1,
                  std::back_inserter(elemRecords[feType]));
      }
    }
    
    // Wedge, Hexa27, Pyra
    std::map< Elem::FEType, std::vector<UInt> >::iterator it, end;
    it = elemRecords.begin();
    end = elemRecords.end();
    
    for( ; it != end; it++ )
    {
      UInt type;
      UInt numElemNodes;
      UInt recSize;
      UInt numRecs;
      std::vector<UInt>& elemRecords = it->second;

      GmshHelper::FEType2ElemType(it->first, type, numElemNodes);
      recSize = 1 + numTags + numElemNodes;
      numRecs = elemRecords.size() / recSize;

      if(ascii_)
      {
        for (UInt i = 0; i < numRecs; i++) {
          (*output_) << elemRecords[i*recSize + 0] << " " << type << " " << numTags << " ";

          for (UInt j = 0; j < numTags; j++) {
            (*output_) << elemRecords[i*recSize + 1 + j] << " ";
//            std::cout << "tag " << elemRecords[i*recSize + 1 + j] << std::endl;
          }

          for (UInt j = 0; j < numElemNodes; j++) {
            (*output_) << elemRecords[i*recSize + 1 + numTags + j] << " ";            
//            std::cout << "elemNode " << elemRecords[i*recSize + 1 + numTags + j] << std::endl;
          }

          (*output_) << std::endl;
          
        }
      }
      else
      {
        UInt dummy;
        
        dummy = uiSwap_.EndianSwapBytes(type);
        (*output_).write((char*) &dummy, sizeof(dummy));
        dummy = uiSwap_.EndianSwapBytes(numRecs);
        (*output_).write((char*) &dummy, sizeof(dummy));
        dummy = uiSwap_.EndianSwapBytes(numTags);
        (*output_).write((char*) &dummy, sizeof(dummy));

        for(UInt i=0, n=elemRecords.size(); i<n; i++)
          elemRecords[i] = uiSwap_.EndianSwapBytes(elemRecords[i]);
          
        (*output_).write((char*) &elemRecords[0], elemRecords.size() * sizeof(elemRecords[0]));

        //        (*output_) << std::endl;
      }
    }

    (*output_) << "$EndElements" << std::endl;
  }

   void SimOutputGmsh::WriteRegions()
   {
     UInt numPhysicalNames = 0;

     StdVector<RegionIdType> volRegions;
     ptGrid_->GetVolRegionIds( volRegions );
     StdVector<RegionIdType> surfRegions;
     ptGrid_->GetSurfRegionIds( surfRegions );
     UInt dim = ptGrid_->GetDim();
     StdVector<std::string> nodeNames;
     ptGrid_->GetListNodeNames(nodeNames);
     StdVector<std::string> elemNames;
     ptGrid_->GetListElemNames(elemNames);

     numPhysicalNames = volRegions.GetSize() + surfRegions.GetSize();

     (*output_) << "$PhysicalNames" << std::endl;
     (*output_) << ( numPhysicalNames + 
                     nodeNames.GetSize() + elemNames.GetSize() ) << std::endl;

     for(UInt i=0; i<volRegions.GetSize(); i++) 
     {
       std::string regionName = ptGrid_->GetRegion().ToString(volRegions[i]);

       // Determine dimension of elements
       StdVector<Elem*> elems;
       ptGrid_->GetElems( elems, volRegions[i] );
       dim = 1000; // nonsense dimension
       for( UInt e = 0, m = elems.GetSize(); e < m; e++ ) {
         UInt edim = Elem::shapes[elems[e]->type].dim;
         
         if(dim != 1000) 
         {
           if(edim != dim)
             EXCEPTION("Physical entities must have same dimension.")
         }
         else
         {
           dim = edim;
         }
       }

       (*output_) << dim << " " << (volRegions[i]+1) << " \"" << regionName << "\"" << std::endl;
     }

     for(UInt i=0; i<surfRegions.GetSize(); i++) 
     {
       std::string regionName = ptGrid_->GetRegion().ToString(surfRegions[i]);

       // Determine dimension of elements
       StdVector<Elem*> elems;
       ptGrid_->GetElems( elems, surfRegions[i] );
       dim = 1000; // nonsense dimension
       for( UInt e = 0, m = elems.GetSize(); e < m; e++ ) {
         UInt edim = Elem::shapes[elems[e]->type].dim;
         
         if(dim != 1000) 
         {
           if(edim != dim)
             EXCEPTION("Physical entities must have same dimension.")
         }
         else
         {
           dim = edim;
         }
       }

       (*output_) << dim << " " << (surfRegions[i]+1) << " \"" << regionName << "\"" << std::endl;
     }

     // Named elems
     for( UInt i = 0, n = elemNames.GetSize(); i < n; i++ ) {
       // Determine dimension of elements
       StdVector<Elem*> elems;
       ptGrid_->GetElemsByName( elems, elemNames[i] );
       dim = 1000; // nonsense dimension
       for( UInt e = 0, m = elems.GetSize(); e < m; e++ ) {
         UInt edim = Elem::shapes[elems[e]->type].dim;
         
         if(dim != 1000) 
         {
           if(edim != dim)
             EXCEPTION("Physical entities must have same dimension.")
         }
         else
         {
           dim = edim;
         }
       }

       (*output_) << dim << " " << (++numPhysicalNames) << " \"" << elemNames[i] << "\"" << std::endl;
     }

     // =============
     //  NODE GROUPS
     // =============
     for( UInt i = 0, n = nodeNames.GetSize(); i < n; i++ ) {
       (*output_) << 0 << " " << (++numPhysicalNames) << " \"" << nodeNames[i] << "\"" << std::endl;
     }
     
     (*output_) << "$EndPhysicalNames" << std::endl;
  }
  
  void SimOutputGmsh::WriteGrid( bool printGridOnly ) {
    // fileName for grid file
    std::string gridFileName = fileName_ + ".msh";

    // -------------------------------------------------------------
    // Section for PrintGridOnly and for writing of first fixed grid
    // -------------------------------------------------------------

    output_ = OpenFile( gridFileName );

    // Write header
    (*output_) << "$MeshFormat" << std::endl;
    (*output_) << "2.1 " << ( ascii_ ? 0 : 1 ) << " " << sizeof(double) << std::endl;

    if(!ascii_)
    {
      UInt oneBinary = 1;
      
      oneBinary = uiSwap_.EndianSwapBytes(oneBinary);

      (*output_).write((char*) &oneBinary, sizeof(oneBinary));
      //      (*output_) << std::endl;
    }
    
    (*output_) << "$EndMeshFormat" << std::endl;

    WriteRegions();
    WriteNodes( output_ );
    WriteCells( output_ );
    
    if ( printGridOnly ){
      output_->close();
      delete output_;
      output_ = NULL;
    }

  }


  // ***********
  //  OpenFile
  // ***********
  std::ofstream * SimOutputGmsh::OpenFile( const std::string& name ) {
    std::string totalName;
    std::ofstream * outFile = NULL;

    // Generate basename for output file
    
    fs::path filePath = dirName_ / name;
    
    if (ascii_) {
      outFile = new std::ofstream(filePath);
    }
    else {
      outFile = new std::ofstream(filePath, std::ofstream::binary);
    }
    if ( !outFile ) {
      EXCEPTION("Could not open file " << totalName
                << " for writing Gmsh output");
    }

    return outFile;
  }

}
