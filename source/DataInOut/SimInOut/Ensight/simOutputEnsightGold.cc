// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdarg.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
#include "boost/date_time/posix_time/posix_time.hpp"

namespace fs = boost::filesystem;

#include "simOutputEnsightGold.hh"

namespace CoupledField {


  SimOutputEnsightGold::SimOutputEnsightGold(std::string fileName, ParamNode * inputNode) :
    SimOutput(fileName, inputNode),
    gridWritten_(false) {
    
    fileName_ = fileName;
    formatName_ = "ensight";
    dirName_ = "simoutput_ensight";
    
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );
    
    elemNames_[ET_LINE2] = "bar2";
    elemNames_[ET_LINE3] = "bar3";
    elemNames_[ET_TRIA3] = "tria3";
    elemNames_[ET_TRIA6] = "tria6";
    elemNames_[ET_QUAD4] = "quad4";
    elemNames_[ET_QUAD8] = "quad8";
    elemNames_[ET_QUAD9] = "quad8";
    elemNames_[ET_TET4] = "tetra4";
    elemNames_[ET_TET10] = "tetra10";
    elemNames_[ET_HEXA8] = "hexa8";
    elemNames_[ET_HEXA20] = "hexa20";
    elemNames_[ET_HEXA27] = "hexa20";
    elemNames_[ET_PYRA5] = "pyramid5";
    elemNames_[ET_PYRA13] = "pyramid5";
    elemNames_[ET_WEDGE6] = "penta6";
    elemNames_[ET_WEDGE15] = "penta15";
    
    actMSStep_ = 1;
    actStep_ = 1;

  }


  SimOutputEnsightGold::~SimOutputEnsightGold() {
  }


  void SimOutputEnsightGold::Init(Grid* ptGrid,
      bool printGridOnly)
  {
    ptGrid_ = ptGrid;
    printGridOnly_ = printGridOnly;
    std::cout << "In Init, printGridOnly_: " << printGridOnly_ << std::endl;
    try {
      pathSep_ = fs::path("/").native_directory_string();
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }    
    
    WriteGrid();    
  }
  
  void SimOutputEnsightGold::InitModule() {
    std::string fileName;
    std::ostringstream strBuffer;

    // concatenate case file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }    

    if(printGridOnly_)
    {
      strBuffer << dirName_ << pathSep_ << fileName_ << ".case";
      fileName = strBuffer.str();

      caseFile_.open( fileName.c_str(), std::ios::binary );
      if ( !caseFile_.good() ) {
        EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
                   << fileName << "!";); 
      }

      caseFile_ << "# Ensight Gold case file written by Simon" << std::endl; 
      caseFile_ << std::endl;
      caseFile_ << "FORMAT" << std::endl;
      caseFile_ << std::endl;
      caseFile_ << std::endl;
      caseFile_ << "type:  ensight gold" << std::endl;
      caseFile_ << std::endl;
      caseFile_ << std::endl;
      caseFile_ << "GEOMETRY" << std::endl;
      caseFile_ << std::endl;
      caseFile_ << std::endl;
      caseFile_ << "model: " << fileName_ << ".geo" << std::endl;
    }
    
    // concatenate geo file name
    strBuffer.str("");
    strBuffer << dirName_ << pathSep_ << fileName_ << ".geo";
    fileName = strBuffer.str();

    geoFile_.open( fileName.c_str(), std::ios::binary );
    if ( !geoFile_.good() ) {
      EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
          << fileName << "!";); 
    }
    
  }  
  
  void SimOutputEnsightGold::WriteGrid()
  {
    // ensure that grid gets only written once
    if(!gridWritten_)
      InitModule();
    else
      return;
    
    //clever little trick from Quake 2 to determine the endian
    //of the current system without depending on a preprocessor define

    uint8_t SwapTest[2] = { 1, 0 };
    bool bigEndianSystem = false;

    if( *(unsigned short *) SwapTest == 1 )
    {
      //little endian
      bigEndianSystem = false;
    }
    else
    {
      //big endian
      bigEndianSystem = true;
    }

    geoFile_ << "This is the 1st description line of the EnSight Gold geometry example" << std::endl;
    geoFile_ << "This is the 2nd description line of the EnSight Gold geometry example" << std::endl;
    geoFile_ << "node id given" << std::endl;
    geoFile_ << "element id given" << std::endl;
    
    // iterate over all regions
    StdVector<RegionIdType> regionIds;
    StdVector<UInt> nodeVec;
    StdVector<Double> xCoordVec, yCoordVec, zCoordVec;
    typedef std::ostream_iterator<UInt> uint_os_iter;
    typedef std::ostream_iterator<Double> double_os_iter;
    std::map<FEType, std::vector<UInt> > elemMap;
    std::map<FEType, std::vector<UInt> > elemIds;
    
    ptGrid_->GetRegionIds(regionIds);
    for ( UInt iReg = 0; iReg < regionIds.GetSize(); iReg++ ) {

      // get region name and elements
      std::string regionName = ptGrid_->RegionIdToName(regionIds[iReg]);
      ptGrid_->GetNodesByRegion( nodeVec, regionIds[iReg] );

      UInt i, j, nNodes;
      nNodes = nodeVec.GetSize();
      UInt nWidth = 10;
      
      regionIdToPart_[regionIds[iReg]] = (iReg+1);
      geoFile_ << "part" << std::endl;
      geoFile_ << std::setw(nWidth);
      geoFile_ << regionIdToPart_[regionIds[iReg]] << std::endl;
      geoFile_ << "Unstructured CFS++ grid." << std::endl;
      geoFile_ << "coordinates" << std::endl;
      geoFile_ << std::setw(nWidth);
      geoFile_ << nNodes << std::endl;
      
      // std::ios::floatfield will allow you to make the desired tabular form,
      // std::ios::fixed sets a fixed number format (dddd.dd), alternatives are
      // std::ios::scientific (that is d.dddddEdd) or 0 (general format)
      geoFile_.setf( std::ios::scientific, std::ios::floatfield );

      // the default precision is 6, let's make it 5
      geoFile_ << std::setprecision( 5 ); // number of digits after the '.'

      // the width of the columns that will be used (in nr of characters)
      
      for ( i = 0; i < nNodes; i++ ) {
        // note that the width needs to be set individually for each output (3)
        geoFile_ << std::setw( nWidth ) << nodeVec[i] << std::endl;
      }
      // std::copy (nodeVec.Begin (), nodeVec.End (), uint_os_iter (geoFile_, "\n"));
      
      xCoordVec.Resize(nNodes);
      yCoordVec.Resize(nNodes);
      zCoordVec.Resize(nNodes);
      for (UInt i = 0, nNodes=nodeVec.GetSize();
           i < nNodes;
           i++) {
        Point p;
        ptGrid_->GetNodeCoordinate(p, nodeVec[i]);
        
        xCoordVec[i] = p[0];
        yCoordVec[i] = p[1];
        zCoordVec[i] = p[2];
      }

      nWidth = 12;
      
      for ( i = 0; i < nNodes; i++ ) {
        geoFile_ << std::setw( nWidth ) << xCoordVec[i] << std::endl;
      }
      for ( i = 0; i < nNodes; i++ ) {
        geoFile_ << std::setw( nWidth ) << yCoordVec[i] << std::endl;
      }
      for ( i = 0; i < nNodes; i++ ) {
        geoFile_ << std::setw( nWidth ) << zCoordVec[i] << std::endl;
      }

      
      
      RegionElems2EnsightElems(regionIds[iReg], elemMap, elemIds);
      std::map<FEType, std::vector<UInt> >::iterator it, end, elIdsIt;
      
      it = elemMap.begin();
      end = elemMap.end();
      elIdsIt = elemIds.begin();
      
      for( ; it != end; it++, elIdsIt++ ) {
        FEType type = it->first;
        std::vector<UInt> elCon = it->second;
        std::vector<UInt> elIds = elIdsIt->second;
        UInt numNodes = NUM_ELEM_NODES[type];
        
        switch(type)
        {
        case ET_QUAD9:
          numNodes = 8;
          break;
        case ET_HEXA27:
          numNodes = 20;
          break;
        default:
          break;
        }
        
        UInt nElems = elIds.size();
      
      
        
        geoFile_ << elemNames_[type] << std::endl;
        nWidth = 10;
        geoFile_ << std::setw( nWidth ) << nElems << std::endl;
        
        for ( i = 0; i < nElems; i++ ) {
          geoFile_ << std::setw( nWidth ) << elIds[i] << std::endl;
        }

        for ( i = 0; i < nElems; i++ ) {
          geoFile_ << std::setw( nWidth ) << elCon[i*numNodes+0];
          for ( j = 1; j < numNodes; j++ ) {
            geoFile_ << " " << elCon[i*numNodes+j];
          }
          geoFile_ << std::endl;          
        }
      }
    }
    
    gridWritten_ = true;
    
//    geoFile_.close();
//    caseFile_.close();   

  }


  void SimOutputEnsightGold::RegisterResult( shared_ptr<BaseResult> sol,
                                             UInt saveBegin, UInt saveInc,
                                             UInt saveEnd,
                                             bool isHistory )
  {
    std::string resultName;
    TimeFileSetInfo tfsi;
    std::stringstream sstr;
    ResultInfo & actInfo = *(sol->GetResultInfo());

    //    std::cout << actInfo.resultName << std::endl;
    std::cout << "Registering output '" << actInfo.resultName
              << "' with saveBegin " << saveBegin
              << ", saveInc " << saveInc 
              << ", saveEnd " << saveEnd << std::endl;
    
    if(isHistory)
    {
      Warning("History results cannot be saved in Ensight format at this time!");
      return;
    }

    if((actInfo.entryType != ResultInfo::SCALAR) &&
       (actInfo.entryType != ResultInfo::VECTOR))
    {
      Warning("Only scalars and vectors can be saved for Ensight at this time!");
      return;
    }

    resultName = actInfo.resultName;
    tfsi.setNumber = resultSetInfo_.size()+1;
    tfsi.actStepIdx = 0;
    tfsi.stepNums.Resize(0);
    tfsi.stepVals.Resize(0);
    
//    for(UInt step=saveBegin; step<=saveEnd; step+=saveInc) {
//      tfsi.stepNums.Push_back(step);
//      tfsi.stepVals.Push_back(0.0);
//    }

    sstr.str("");
    sstr << resultName << "_ms" << actMSStep_ << "_s";
    tfsi.baseFileName = sstr.str();
    
    resultSetInfo_[resultName] = tfsi;
    
    
  }

  void SimOutputEnsightGold::BeginMultiSequenceStep( UInt step,
                                                     BasePDE::AnalysisType type,
                                                     UInt numSteps  )
  {
    std::cout << "In BeginMultiSequenceStep" << std::endl;

    std::string fileName;
    std::ostringstream strBuffer;

    if((type != BasePDE::HARMONIC) && (type != BasePDE::TRANSIENT))
    {
      Warning("Analysis Type not supported by Ensight format at this time!");
      return;
    }
    
    analysisType_ = type;
    
    actMSStep_ = step;

    // concatenate case file name
    try {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }    

    strBuffer << dirName_ << pathSep_ << fileName_ << "_ms"
              << step << ".case";
    fileName = strBuffer.str();

    caseFile_.open( fileName.c_str(), std::ios::binary );
    if ( !caseFile_.good() ) {
      EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
                 << fileName << "!";); 
    }

    caseFile_ << "# Ensight Gold case file written by Simon" << std::endl; 
    caseFile_ << std::endl;
    caseFile_ << "FORMAT" << std::endl;
    caseFile_ << std::endl;
    caseFile_ << std::endl;
    caseFile_ << "type:  ensight gold" << std::endl;
    caseFile_ << std::endl;
    caseFile_ << std::endl;
    caseFile_ << "GEOMETRY" << std::endl;
    caseFile_ << std::endl;
    caseFile_ << std::endl;
    caseFile_ << "model: " << fileName_ << ".geo" << std::endl;


    std::map<std::string, TimeFileSetInfo >::iterator it, end;
    it = resultSetInfo_.begin();
    end = resultSetInfo_.end();
    std::string resultName;
    std::stringstream sstr;
    
    // Iterate over all results and open result files.
    for( ; it != end; it++ ){
      resultName = it->first;
      sstr.str("");
      
      sstr << resultName << "_ms" << actMSStep_ << "_s";
      it->second.baseFileName = sstr.str();
//
//      sstr.str("");
//      sstr << dirName_ << pathSep_ << it->second.fileName;
//      it->second.resultFile = new std::ofstream( sstr.str().c_str(),
//                                                 std::ios::binary );
//      if ( !it->second.resultFile->good() ) {
//        EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
//                   << it->second.fileName << "!";); 
//      }
//      
//      (*it->second.resultFile) << "description line 1" << std::endl;
    }
    
  }
  
  void SimOutputEnsightGold::BeginStep( UInt stepNum, Double stepVal )
  {
    std::cout << "In BeginStep" << std::endl;
    actStep_ = stepNum;
    actStepVal_ = stepVal;
  }

  void SimOutputEnsightGold::AddResult( shared_ptr<BaseResult> sol )
  {
    std::cout << "In AddResult" << std::endl;   
    ResultInfo & actInfo = *(sol->GetResultInfo());
    std::string resultName = actInfo.resultName;
    UInt actStepIdx;
    
    // At the moment this writer only supports scalars and vectors
    if((actInfo.entryType != ResultInfo::SCALAR) &&
       (actInfo.entryType != ResultInfo::VECTOR))
    {
      return;
    }
   
    std::cout << actInfo.resultName << std::endl;

    // Add result to temporary vector 
    actStepIdx = resultSetInfo_[resultName].actStepIdx;
    
    if(actStepIdx == resultSetInfo_[resultName].stepNums.GetSize()) {
      resultSetInfo_[resultName].stepNums.Push_back(actStep_);
      resultSetInfo_[resultName].stepVals.Push_back(actStepVal_);
      resultSetInfo_[resultName].entityType = actInfo.definedOn;
      resultSetInfo_[resultName].entryType = actInfo.entryType;
    }
    
    resultSetInfo_[resultName].regionResults.Push_back(sol);
  }

  void SimOutputEnsightGold::FinishStep( )
  {
    std::cout << "In FinishStep" << std::endl;
    std::map<std::string, TimeFileSetInfo >::iterator it, end;
    it = resultSetInfo_.begin();
    end = resultSetInfo_.end();
    UInt nWidth;
    std::stringstream sstr;
    std::string fileName;
    
    // über alle ergebnisse iterieren
    // nur index hochzählen, wenn stepnum passt.
    for( ; it != end; it++ ){
      UInt actStepIdx = it->second.actStepIdx;
      if( it->second.stepNums[actStepIdx] != actStep_ )
        continue;
      
      it->second.actStepIdx++;
      
      sstr.str("");
      
      if(analysisType_ == BasePDE::TRANSIENT)
        sstr << it->second.baseFileName << actStep_ << ".dat";
      else
        sstr << it->second.baseFileName << actStep_ << "_real.dat";        
      fileName = sstr.str();
      
      
      sstr.str("");
      sstr << dirName_ << pathSep_ << fileName;
      std::ofstream resFile( sstr.str().c_str(), std::ios::binary );
      if ( !resFile.good() ) {
        EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
                   << fileName << "!";); 
      }
      
      resFile << "description line 1" << std::endl;
      
      // Iterate over all regions/parts and write results
      StdVector< shared_ptr<BaseResult> >::iterator resIt, resEnd;
      resIt = it->second.regionResults.Begin();
      resEnd = it->second.regionResults.End();
      
      if(analysisType_ == BasePDE::TRANSIENT)
      {
        for( ; resIt != resEnd; resIt++ ) {
          StdVector< shared_ptr<BaseResult> > resVec;
          Vector<Double> partSol;
          Vector<Double> solIt, solEnd;
          resVec.Push_back(*resIt);
          std::string regionName = (*resIt)->GetEntityList()->GetName(); 
          RegionIdType regionId = ptGrid_->RegionNameToId(regionName);

          FillGlobalVec<Double>(partSol, resVec, it->second.entityType );

          nWidth = 10;

          resFile << "part" << std::endl;
          resFile << std::setw(nWidth) << regionIdToPart_[regionId] << std::endl;
          resFile << "coordinates" << std::endl;

          nWidth = 12;
          resFile.setf( std::ios::scientific, std::ios::floatfield );
          resFile << std::setprecision( 5 ); // number of digits after the '.'

          for( UInt i=0, n=partSol.GetSize(); i<n; i++ )
          {
            resFile << std::setw(nWidth) << partSol[i] << std::endl;
          }
        }
      }
      
      if(analysisType_ == BasePDE::HARMONIC)
      {
        sstr.str("");
        sstr << it->second.baseFileName << actStep_ << "_imag.dat";        
        fileName = sstr.str();


        sstr.str("");
        sstr << dirName_ << pathSep_ << fileName;
        std::ofstream resFileImag( sstr.str().c_str(), std::ios::binary );
        if ( !resFileImag.good() ) {
          EXCEPTION ("SimOutputEnsightGold unable to open output file '" \
              << fileName << "!";); 
        }

        resFileImag << "description line 1" << std::endl;

        for( ; resIt != resEnd; resIt++ ) {
          StdVector< shared_ptr<BaseResult> > resVec;
          Vector<Complex> partSol;
          Vector<Complex> solIt, solEnd;
          resVec.Push_back(*resIt);
          std::string regionName = (*resIt)->GetEntityList()->GetName(); 
          RegionIdType regionId = ptGrid_->RegionNameToId(regionName);

          FillGlobalVec<Complex>(partSol, resVec, it->second.entityType );

          nWidth = 10;

          // real part
          resFile << "part" << std::endl;
          resFile << std::setw(nWidth) << regionIdToPart_[regionId] << std::endl;
          resFile << "coordinates" << std::endl;

          // imag part
          resFileImag << "part" << std::endl;
          resFileImag << std::setw(nWidth) << regionIdToPart_[regionId] << std::endl;
          resFileImag << "coordinates" << std::endl;

          nWidth = 12;
          resFile.setf( std::ios::scientific, std::ios::floatfield );
          resFile << std::setprecision( 5 ); // number of digits after the '.'
          resFileImag.setf( std::ios::scientific, std::ios::floatfield );
          resFileImag << std::setprecision( 5 ); // number of digits after the '.'

          
          for( UInt i=0, n=partSol.GetSize(); i<n; i++ )
          {
            resFile << std::setw(nWidth) << partSol[i].real() << std::endl;
            resFileImag << std::setw(nWidth) << partSol[i].imag() << std::endl;
          }
        }
      }

      resFile.close();  
    }
  }

  void SimOutputEnsightGold::FinishMultiSequenceStep( )
  {
    std::cout << "In FinishMultiSequenceStep" << std::endl;

    std::map<std::string, TimeFileSetInfo >::iterator it, end;
    it = resultSetInfo_.begin();
    end = resultSetInfo_.end();
    
    // Iterate over all results and close result files.
    for( ; it != end; it++ ){
//      it->second.resultFile->close();
//      delete it->second.resultFile;
    }
  }
  
  void SimOutputEnsightGold::Finalize()
  {
    std::cout << "In Finalize" << std::endl;

    geoFile_.close();
    caseFile_.close();   
  }


  // =========================================================================
  // MISCELLANEOUS METHODS
  // =========================================================================


  // ***************
  //   Type2ptElem
  // ***************
  UInt SimOutputEnsightGold::RegionElems2EnsightElems(const RegionIdType region,
      std::map<FEType, std::vector<UInt> > & elemMap,
      std::map<FEType, std::vector<UInt> > & elemIds) {

    StdVector<Elem*> elemVec;
    FEType type;
    UInt nElems, elemId;
    UInt i=0;
    StdVector<UInt>::iterator it, end;
    StdVector<UInt> nodeVec;
    std::map<UInt, UInt> globToLocNodeMap;
    
    ptGrid_->GetElems(elemVec,region);
    ptGrid_->GetNodesByRegion( nodeVec, region );
    nElems = elemVec.GetSize();

    elemMap.clear();
    elemIds.clear();
    
    it = nodeVec.Begin();
    end = nodeVec.End();
    for( i=1; it != end; it++, i++ )
      globToLocNodeMap[*it] = i; 
    
    for(i=0; i<nElems; i++)
    {
      type = elemVec[i]->ptElem->feType();
      elemId = elemVec[i]->elemNum;
      it = elemVec[i]->connect.Begin();
      end = elemVec[i]->connect.End();
      
      elemIds[type].push_back(elemId);

      switch(type)
      {
      case ET_LINE3:
        elemMap[type].push_back(globToLocNodeMap[*(it)]);        
        elemMap[type].push_back(globToLocNodeMap[*(it+2)]);        
        elemMap[type].push_back(globToLocNodeMap[*(it+1)]);        
        break;
      case ET_QUAD9:
        for( ; it != (it+8); it++ )
          elemMap[type].push_back(globToLocNodeMap[*it]);
        break;
      case ET_HEXA27:
        for( ; it != (it+20); it++ )
          elemMap[type].push_back(globToLocNodeMap[*it]);
        break;
      default:
        for( ; it != end; it++ )
          elemMap[type].push_back(globToLocNodeMap[*it]);
        
//        std::copy(elemVec[i]->connect.Begin(),
//                  elemVec[i]->connect.End(),
//                  std::back_inserter(elemMap[type]));
        continue;
      }

    }
    
    // This place should never be reached!
    return 0;
  }

}

