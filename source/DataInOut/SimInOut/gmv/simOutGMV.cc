// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <complex>
#include <ctime>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>
namespace fs = boost::filesystem;

#include <def_cfs_stats.hh>

#include <DataInOut/simInput.hh>
#include <DataInOut/Logging/cfslog.hh>
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"

#include "simOutGMV.hh"

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(simOutputGMV)
  DEFINE_LOG(simOutputGMV, "SimOutputGMV")

  #define CODE_NAME "CFS++"
  #define CODE_VER CFS_VERSION

   //*****************
  //   Constructor
  //*****************
  SimOutputGMV::SimOutputGMV( const std::string fileName,
                              ParamNode * outputNode ) :
    SimOutput( fileName, outputNode )
  {


    // Initialize variables
    formatName_ = "gmv";
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    std::ostringstream strBuffer;

    stepNumOffset_ = 0;
    stepValOffset_ = 0.0;
    dirName_ = "simoutput_gmv";
    fileName_ = fileName;

    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }
    
    // Does the grid change over time, or can we use a fixed grid
    fixedgrid_ = true;
    firstGridWritten_ = false;
    output = NULL;
    printUnit_ = false;
    ascii_ = false;
    charOutSize_ = 32;


    // Change defaults according to XML file
    if(myParam_->Get("binaryFormat", false)->AsString() == "no")
    {
        ascii_ = true;
    }

    if(myParam_->Get("fixedGrid", false)->AsString() == "no")
    {
        fixedgrid_ = false;
    }
    
    charOutSize_ = 32;


  }


  // **********************
  //   Default Destructor
  // **********************
  SimOutputGMV::~SimOutputGMV() {


  }
  
  void SimOutputGMV::Init( Grid* ptGrid, bool printGridOnly ) {
    
    ptGrid_ = ptGrid;

    if( printGridOnly ) 
      WriteGrid( true );

  }


  void SimOutputGMV::BeginMultiSequenceStep( UInt step,
                                             AnalysisType type,
                                             UInt numSteps )
  {

    currAnalysis_ = type;
    currMsStep_ = step;

  }
  
  void SimOutputGMV::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd,
                                     bool isHistory )
  {

    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGMV) << "Registering output '" << actDof.resultName 
                          << "' with saveBegin " << saveBegin
                          << ", saveInc " << saveInc 
                          << ", saveEnd " << saveEnd;
  }


  //! Begin single analysis step
  void SimOutputGMV::BeginStep( UInt stepNum, Double stepVal )
  {

    resultMap_.clear();
    
    actStep_ = stepNum;
    actStepVal_ = stepVal;
    if( currAnalysis_ == TRANSIENT ||
        currAnalysis_ == STATIC  ) { 
      actStep_ += stepNumOffset_;
      actStepVal_ += stepValOffset_;
    }

  }


  //! Add result to current step
  void SimOutputGMV::AddResult( shared_ptr<BaseResult> sol )
  {
    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGMV) << "Adding result '" 
                          << actDof.resultName << "'";
      
    resultMap_[sol->GetResultInfo()->resultName].Push_back(sol);
  }


  //! End single analysis step
  void SimOutputGMV::FinishStep( )
  {

    LOG_TRACE(simOutputGMV) << "Starting to finish Step";
    
    // Check, if any result at all is registered. 
    // If not leave!
    if( resultMap_.size() == 0 )
      return;

    // ----------------------
    // Open new file
    // ----------------------
    std::ostringstream strBuffer;
    std::string analysisString;
    Enum2String( currAnalysis_, analysisString );
    strBuffer <<  fileName_ << "-" << analysisString
              << "-" << currMsStep_ << ".gmv";
    if ( actStep_ < 10 ) strBuffer << "000";
    else if ( actStep_ < 100 ) strBuffer << "00";
    else if ( actStep_ < 1000 ) strBuffer << "0";
    else if ( actStep_ > 10000 ) {
      EXCEPTION("Number of gmv files exceeds 9999!");
    }
    
    strBuffer << actStep_;
    output = OpenFile( strBuffer.str() );

    // Print Grid
    WriteGrid( false );


    // -------------------------
    // Iterate over all results
    // -------------------------
    std::string title;
    
    // iterate over all result types
    ResultMapType::iterator it = resultMap_.begin();
    for( ; it != resultMap_.end(); it++ ) {
      
      // get result info object and results for current result type
      //ResultInfo & actDof = *(it->first);
      ResultInfo & actInfo = *(it->second[0]->GetResultInfo());
      const StdVector<shared_ptr<BaseResult> > actResults =
        it->second;

      title = actInfo.resultName;

      ResultInfo::EntryType entryType =  actInfo.entryType;
      ResultInfo::EntityUnknownType entityType = actInfo.definedOn;

      // check if result is defined on nodes or elements
      StdVector<std::string> & dofNames = actInfo.dofNames;  
      if(  actInfo.definedOn != ResultInfo::NODE &&
           actInfo.definedOn != ResultInfo::PFEM &&
           actInfo.definedOn != ResultInfo::ELEMENT &&
           actInfo.definedOn != ResultInfo::SURF_ELEM ) {
        Warning( "GMV can only write results on element and nodes",
                 __FILE__, __LINE__ );
        continue;
      }

      StdVector<std::string> myDofNames;
      myDofNames.Resize(dofNames.GetSize());
      
      for(UInt i=0; i < dofNames.GetSize(); i++)
      {
        myDofNames[i] = dofNames[i];
        if(printUnit_)
          myDofNames[i] += "_(" + actInfo.unit + ")";
      }

      LOG_DBG(simOutputGMV) << "Writing result '" << title << "'";
      
      if( actResults[0]->GetEntryType() == EntryType::DOUBLE ) {
        Vector<Double> gSol;
        FillGlobalVec<Double>(gSol, actResults, entityType );
        WriteNodeElemDataTrans( gSol, myDofNames, title, entryType,
                                entityType, actStepVal_ );        
      } else {
        Vector<Complex> gSol;
        FillGlobalVec<Complex>(gSol, actResults, entityType );
        WriteNodeElemDataHarm( gSol, myDofNames, title, entryType, entityType,
                               actStepVal_, actInfo.complexFormat );
      }
    }

  
    // print end of file and close file
    PrintFileEpilog( output, true );
    output->close();
    delete output;
  }

  //! End multisequence step
  void SimOutputGMV::FinishMultiSequenceStep( ) {
    // set offset for step value and number to last values
    if( currAnalysis_ == TRANSIENT ||
          currAnalysis_ == STATIC ) {
      stepNumOffset_ = actStep_;
      stepValOffset_ = actStepVal_;
    }

  }


  //! Finalize the output
  void SimOutputGMV::Finalize() {

    // nothing to do here
  }

  void SimOutputGMV::
  WriteNodeElemDataTrans( const Vector<Double> & var, 
                          const StdVector<std::string> & dofNames,
                          const std::string& name, 
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {
    
    // get number of entities
    UInt numEnt = 0;
    UInt loc;
    
    if ( entityType == ResultInfo::NODE ||
         entityType == ResultInfo::PFEM ) {
      loc = 1;
      numEnt = ptGrid_->GetNumNodes();
    } else {
      loc = 0;
      numEnt = ptGrid_->GetNumElems();
    }

    // write Result header
    if ( entryType == ResultInfo::SCALAR ) {
      WriteVariable(name, loc, var);
    } else if( entryType == ResultInfo::VECTOR ) {
      WriteVector(name, loc, var, dofNames);
    } else if( entryType == ResultInfo::TENSOR ) {
      WriteVector(name, loc, var, dofNames);
    }
  } 
    
  void SimOutputGMV::
  WriteNodeElemDataHarm( const Vector<Complex> & var, 
                         const StdVector<std::string> & dofNames,
                         const std::string name, 
                         const ResultInfo::EntryType entryType,
                         ResultInfo::EntityUnknownType entityType,
                         const Double freq, 
                         const ComplexFormat outputFormat )
  {

    // get number of entities
    UInt numEnt = 0;
    UInt loc;
    
    if ( entityType == ResultInfo::NODE ||
         entityType == ResultInfo::PFEM ) {
      loc = 1;
      numEnt = ptGrid_->GetNumNodes();
    } else {
      loc = 0;
      numEnt = ptGrid_->GetNumElems();
    }

    UInt numDofs = dofNames.GetSize();
    
    // write Result header
    if ( entryType == ResultInfo::SCALAR ) {
      Vector<Double> var1;
      Vector<Double> var2;
      std::string name1;
      std::string name2;
      
      var1.Resize(numEnt);
      var2.Resize(numEnt);
      
      if(outputFormat == REAL_IMAG)
      {
        for(UInt i=0; i<numEnt; i++)
        {
          var1[i] = var[i].real();
          var2[i] = var[i].imag();
        }
        name1 = name + "_real";
        name2 = name + "_imag";
      }
      else
      {
        for(UInt i=0; i<numEnt; i++)
        {
          var1[i] = std::abs(var[i]);
          var2[i] = CPhase(var[i]);
        }
        name1 = name + "_ampl";
        name2 = name + "_phase";
      }
      
      WriteVariable(name1, loc, var1);
      WriteVariable(name2, loc, var2);
    } else if( entryType == ResultInfo::VECTOR ||
               entryType == ResultInfo::TENSOR ) {
      Vector<Double> vars1;
      Vector<Double> vars2;
      StdVector< std::string > names1;
      StdVector< std::string > names2;

      vars1.Resize(numDofs*numEnt);
      vars2.Resize(numDofs*numEnt);
      names1.Resize(numDofs);
      names2.Resize(numDofs);
      
      if(outputFormat == REAL_IMAG)
      {
        for(UInt j=0; j<numDofs; j++)
        {
          for(UInt i=0; i<numEnt; i++)
          {
            vars1[i*numDofs+j] = var[i*numDofs+j].real();
            vars2[i*numDofs+j] = var[i*numDofs+j].imag();
          }
          names1[j] = dofNames[j] + "_real";
          names2[j] = dofNames[j] + "_imag";
        }
      }
      else
      {
        for(UInt j=0; j<numDofs; j++)
        {
          for(UInt i=0; i<numEnt; i++)
          {
            vars1[i*numDofs+j] = std::abs(var[i*numDofs+j]);
            vars2[i*numDofs+j] = CPhase(var[i*numDofs+j]);
          }
          names1[j] = dofNames[j] + "_ampl";
          names2[j] = dofNames[j] + "_phase";
        }
      }

      WriteVector(name, loc, vars1, names1);
      WriteVector(name, loc, vars2, names2);
    }
  }
  

  // ******************************
  //   WriteVariable
  // ******************************
  void SimOutputGMV::WriteVariable( const std::string name,
                                    const UInt dataType,
                                    const Vector<Double> & var) 
  {
    std::string outBuf = name;
    TruncateString(outBuf, charOutSize_);

    if (ascii_) {
      (*output) << std::endl;
      (*output) << "variable " << outBuf << " "
                << dataType << std::endl;

      UInt i;
      for (i=0; i<var.GetSize(); i++)
        (*output) << var[i] << " ";

      (*output) << std::endl << "endvars" << std::endl;
    }
    else {
      outBuf = "variable";
      TruncateString(outBuf, 8);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));

      outBuf = name;
      TruncateString(outBuf, charOutSize_);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));
      output->write((char*)&dataType,sizeof(UInt));

      const Double * ptvar=var.GetPointer();
      UInt len=var.GetSize();
      output->write((char*)ptvar,len * sizeof(double));

      outBuf = "endvars";
      TruncateString(outBuf, 8);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));
    }

  }

  // ******************************
  //   WriteVector
  // ******************************
  void SimOutputGMV::WriteVector( const std::string name,
                                  const UInt dataType,
                                  const Vector<Double> & var,
                                  const StdVector<std::string> & dofNames) 
  {
    UInt numDOFs = dofNames.GetSize();
    UInt numEntities = var.GetSize() / numDOFs;
    std::string outBuf = name;
    
    TruncateString(outBuf, charOutSize_);
    
    if (ascii_) {
      (*output) << std::endl;
      (*output) << "vectors " << outBuf << " "
                << dataType << " "
                << numDOFs << " "
                << 1 << std::endl;

      for(UInt i=0; i < numDOFs; i++)
      {
        (*output) << name << "_" << dofNames[i] << std::endl;
      }

      for(UInt i=0; i < numDOFs; i++)
      {
        for (UInt j=0; j<numEntities; j++)
          (*output) << var[j*numDOFs+i] << " ";

        (*output) << std::endl;
      }
      
      (*output) << "endvect" << std::endl;
    }
    else {
      UInt dofFlag = 1;

      outBuf = "vectors";
      TruncateString(outBuf, 8);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));

      outBuf = name;
      TruncateString(outBuf, charOutSize_);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));
      output->write((char*)&dataType,sizeof(UInt));
      output->write((char*)&numDOFs,sizeof(UInt));
      output->write((char*)&dofFlag,sizeof(UInt));

      for(UInt i=0; i < numDOFs; i++)
      {
        outBuf = name + "_" + dofNames[i];
        TruncateString(outBuf, charOutSize_);
        output->write(outBuf.c_str(),outBuf.length()*sizeof(char));
      }

      const Double * ptvar;
      Vector<Double> comp;
      
      comp.Resize(numEntities);

      for(UInt i=0; i < numDOFs; i++)
      {
        for (UInt j=0; j<numEntities; j++)
        {
          comp[j] = var[j*numDOFs+i];
        }

        ptvar=comp.GetPointer();
        output->write((char*)ptvar,numEntities * sizeof(double));
      }

      outBuf = "endvect";
      TruncateString(outBuf, 8);
      output->write(outBuf.c_str(),outBuf.length()*sizeof(char));
    }

  }


  // **************
  //   WriteNodes
  // **************
  void SimOutputGMV::WriteNodes( std::ofstream * gridFile ) {


    // write keyword
    (*gridFile) << "nodev   ";

    //get and write number of nodes on the level
    UInt numnodes = ptGrid_->GetNumNodes();
    if (ascii_)
      (*gridFile) << numnodes << std::endl;
    else
      gridFile->write((char*)&numnodes,sizeof(UInt));

    //get and write coodinates of nodes
    UInt i;

    // write x,y,z-coordinate
    for ( i = 1; i <= numnodes; i++) {
      Point p;
      
      ptGrid_->GetNodeCoordinate(p, i);  
        
      if (ascii_) {
        (*gridFile) << " " << p[0] << " " << p[1] << " "
                  << p[2] << std::endl;
      }
      else {
        gridFile->write((char*)&p[0],sizeof(Double));
        gridFile->write((char*)&p[1],sizeof(Double));
        gridFile->write((char*)&p[2],sizeof(Double));
      }
    }

  }


  // **************
  //   WriteCells
  // **************
  void SimOutputGMV::WriteCells( std::ofstream * gridFile ) {


    // write keyword
    (*gridFile) << "cells   ";

    // read information about number of elements 
    UInt numelem = ptGrid_->GetNumElems();

    if (ascii_)
      (*gridFile) << numelem << std::endl;
    else
      gridFile->write((char*)&numelem,sizeof(UInt));

    UInt i;
    StdVector<RegionIdType> regionIds;
    ptGrid_->GetRegionIds(regionIds);
    UInt numRegions = ptGrid_->GetNumRegions();

    std::vector<UInt> connect;
    UInt elemNum, numNodes;
    RegionIdType elemRegion;
    FEType  elemType;
    std::string gmvElemName;
    connect.resize(100);
    std::vector<UInt> elemRegions;
    StdVector<std::string> regionNames;
    
    ptGrid_->GetRegionNames(regionNames);
    
    
    elemRegions.resize(numelem);

    for ( i = 0; i < numelem; i++ ) {
      elemNum = i+1;
      ptGrid_->GetElemData(elemNum, elemType, elemRegion, &connect[0]);
      numNodes = NUM_ELEM_NODES[elemType];
      
      ElemType2GMVElemId(elemType, gmvElemName);

      elemRegions[i] = elemRegion+1;

      switch(elemType)
      {
      case ET_LINE3:
        connect[3] = connect[1];
        connect[1] = connect[2];
        connect[2] = connect[3];          
        break;
      default:
        break;
      }
    
      if (ascii_)
      {
        (*gridFile) << gmvElemName << " " << numNodes << std::endl;
        UInt j;
        for (j=0; j < numNodes; j++)
          (*gridFile) << " " << connect[j];
        (*gridFile) << std::endl;
      }
      else {
        (*gridFile) << gmvElemName;
        gridFile->write((char*)&numNodes,sizeof(UInt));
        gridFile->write((char*)&connect[0],numNodes * sizeof(UInt));
      }
    }

    // Write Regions -> Materials
    UInt aux = 0;
    if (ascii_)
      (*gridFile) << "material " << numRegions << " 0" << std::endl;
    else {
      (*gridFile) << "material";
      gridFile->write((char*)&numRegions,sizeof(UInt));
      gridFile->write((char*)&aux,sizeof(UInt));
    }

    for(i=0; i<numRegions; i++)
    {
      TruncateString(regionNames[i], charOutSize_);

      if (ascii_) {
        (*gridFile) << regionNames[i] << std::endl;
      } else {
        gridFile->write(&regionNames[i][0],charOutSize_*sizeof(char));
      }

    }

    // write for each element the according regionID
    for (i=0; i<numelem; i++) {

      if (ascii_)
        (*gridFile) << elemRegions[i] << " ";
      else {
        gridFile->write((char*)&elemRegions[i],sizeof(UInt));
      }
    }

    if (ascii_)
      (*gridFile) << std::endl;
  }

  void SimOutputGMV::ElemType2GMVElemId(FEType et, std::string & id)
  {
    switch(et)
    {
    case ET_LINE2:
      id       = "line    ";
      break;
    case ET_LINE3:
      id       = "3line   ";
      break;
    case ET_TRIA3:
      id       = "tri     ";
      break;
    case ET_TRIA6:
      id       = "6tri    ";
      break;
    case ET_QUAD4:
      id       = "quad    ";
      break;
    case ET_QUAD8:
      id       = "8quad   ";
      break;
    case ET_TET4:
      id       = "tet     ";
      break;
    case ET_TET10:
      id       = "ptet10  ";
      break;
    case ET_HEXA8:
      id       = "phex8   ";
      break;
    case ET_HEXA20:
      id       = "phex20  ";
      break;
    case ET_HEXA27:
      id       = "phex27  ";
      break;
    case ET_PYRA5:
      id       = "ppyrmd5 ";
      break;
    case ET_PYRA13:
      id       = "ppyrmd13";
      break;
    case ET_WEDGE6:
      id       = "pprism6 ";
      break;
    case ET_WEDGE15:
      id       = "pprism15";
      break;
    default:
      id       = "general ";
      break;
    }

  }

  // **************************************
  //   Write named entities (nodes, names)
  // **************************************
  void SimOutputGMV::WriteNamedEntities( std::ofstream * gridFile ) {


    StdVector<std::string> nodeNames, elemNames;
    StdVector<UInt> nodeNumbers;
    StdVector<Elem*> elems;
    UInt entType = 0;
    UInt numNamedEntities = 0;

    ptGrid_->GetListNodeNames(nodeNames);

    ptGrid_->GetListElemNames(elemNames);

    // Begin group section
    if (ascii_) {
      (*gridFile) << "groups" << std::endl;
    } else {
      (*gridFile) << "groups  ";
    }
    
    // write named nodes
    entType = 1;
    numNamedEntities = nodeNames.GetSize();
    for( UInt i = 0; i < numNamedEntities; i++ ) {
      
      // get nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);
      UInt numNodes = nodeNumbers.GetSize();

      // write name and number of nodes
      TruncateString(nodeNames[i], charOutSize_);

      if (ascii_) {
        (*gridFile) << nodeNames[i] <<  " " << entType << " " 
                  << numNodes;
      } else {
        gridFile->write(&nodeNames[i][0],charOutSize_*sizeof(char));
        gridFile->write((char*)&entType,sizeof(UInt));
        gridFile->write((char*)&numNodes,sizeof(UInt));
      }

      // write node numbers itself
      if (ascii_) {
        for( UInt iNode = 0; iNode < numNodes; iNode++) {
          (*gridFile) << " " << nodeNumbers[iNode];
        }
      } else {
        for( UInt iNode = 0; iNode < numNodes; iNode++) {
          gridFile->write((char*)&nodeNumbers[iNode],sizeof(UInt));
        }
      }
      
      if (ascii_) {
        (*gridFile) << std::endl;
      }
    }

    // write named elems
    entType = 0;
    numNamedEntities = elemNames.GetSize();
    for( UInt i = 0; i < numNamedEntities; i++ ) {
      
      // get nodes
      ptGrid_->GetElemsByName(elems, elemNames[i]);
      UInt numElems = elems.GetSize();

      // write name and number of nodes
      TruncateString(elemNames[i], charOutSize_);

      if (ascii_) {
        (*gridFile) << elemNames[i] <<  " " << entType << " " 
                  << numElems;
      } else {
        gridFile->write(&elemNames[i][0],charOutSize_*sizeof(char));
        gridFile->write((char*)&entType,sizeof(UInt));
        gridFile->write((char*)&numElems,sizeof(UInt));
      }

      // write elem numbers itself

      if (ascii_) {
        for( UInt iElem = 0; iElem < numElems; iElem++) {
          (*gridFile) << " " << elems[iElem]->elemNum;
        }
      } else {
        for( UInt iElem = 0; iElem < numElems; iElem++) {
          gridFile->write((char*)&elems[iElem]->elemNum,sizeof(UInt));
        }
      }
      
      if (ascii_) {
        (*gridFile) << std::endl;
      }
    }

    if (ascii_) {
      (*gridFile) << "endgrp" << std::endl;
    } else {
      (*gridFile) << "endgrp  ";
    }
      
    
  }
  

  void SimOutputGMV::WriteGrid( bool printGridOnly ) {


    // fileName for external grid file
    std::string extGridFileName = fileName_+ "_GRID.gmv";

    // -------------------------------------------------------------
    // Section for PrintGridOnly and for writing of first fixed grid
    // -------------------------------------------------------------
      
    if ( printGridOnly ||  
         ( fixedgrid_ == true &&
           firstGridWritten_ == false ) ){
           
      std::ofstream * gridFile  = OpenFile( extGridFileName );

      WriteNodes( gridFile );
      WriteCells( gridFile );
      WriteNamedEntities( gridFile );

      PrintFileEpilog( gridFile, false );
      gridFile->close();
      delete gridFile;

      firstGridWritten_ = true;

    }


    // if only grid is to be written ->leave
    if( printGridOnly )
      return;
      
    // ---------------------------
    // Section for new time / frequency step
    // ----------------------------
    
    
    if (fixedgrid_){
      if (ascii_)
        (*output) << "nodev fromfile \"" << extGridFileName
                  <<"\""<< std::endl;
      else 
          (*output) << "nodev   fromfile\"" << extGridFileName <<"\"";
      
      if (ascii_)
        (*output) << "cells fromfile \"" << extGridFileName 
                  <<"\""<< std::endl;
      else 
        (*output) << "cells   fromfile\"" << extGridFileName <<"\"";
      
      if (ascii_)
        (*output) << "groups fromfile \"" << extGridFileName
                  <<"\""<< std::endl;
      else 
        (*output) << "groups  fromfile\"" << extGridFileName <<"\"";

      if (ascii_)
        (*output) << "material fromfile \"" << extGridFileName
                  <<"\""<< std::endl;
      else 
        (*output) << "materialfromfile\"" << extGridFileName <<"\"";
      
    } else {
      WriteNodes( output );
      WriteCells( output );
      WriteNamedEntities( output );
    }
    
  }


  // ***********
  //  OpenFile
  // ***********
  std::ofstream * SimOutputGMV::OpenFile( const std::string& name ) {


    std::string totalName;
    std::ofstream * outFile = NULL;

    // Generate basename for output file
    totalName.append( dirName_ );
    std::string pathsep = fs::path("/").native_directory_string();
    totalName.append( pathsep );
    totalName.append( name );
    
    if (ascii_) {
      outFile = new std::ofstream(totalName.c_str());
    }
    else {
      outFile = new std::ofstream(totalName.c_str(), std::ofstream::binary);
    }
    if ( !outFile ) {
      EXCEPTION("Could not open file " << totalName
                << " for writing GMV output");
    }

    // Write header, problem time and step number
    if (ascii_) {
      (*outFile) << "gmvinput ascii" << std::endl;
    }
    else {
      (*outFile) << "gmvinput" << "iecxi4r8";
    }

    return outFile;
  }

  void SimOutputGMV::PrintFileEpilog( std::ofstream * outFile,
                                      bool printStepInfo ) {
    
    // In the end, add comment 
    struct tm* time_tm;
    time_t time_t;
    char buffer[64];
    std::string dummy;
    
    time(&time_t);
    time_tm = localtime(&time_t);
    strftime(buffer, sizeof(buffer), "%m/%d/%y", time_tm);
    
    // If file was already open, write end of variables
    // and the time and step number
    if (ascii_) {
      if( printStepInfo ) {
        (*outFile) << "probtime " << actStepVal_ << std::endl;
        (*outFile) << "cycleno " << actStep_ << std::endl;
      }
      (*outFile) << "codename " << CODE_NAME << std::endl;
      (*outFile) << "codever " << CODE_VER << std::endl;
      (*outFile) << "simdate " << buffer << std::endl;
      (*outFile) << "endgmv ";
      (*outFile).flush();
    }
    else {
      if( printStepInfo ) {
        (*outFile) << "probtime";
        outFile->write( (char*)&actStepVal_, sizeof(Double) );
        (*outFile) << "cycleno ";
        outFile->write( (char*)&actStep_, sizeof(UInt) );
      }
      (*outFile) << "codename";
      dummy = CODE_NAME;
      TruncateString(dummy, 8);
      outFile->write( dummy.c_str(), 8 );
      (*outFile) << "codever ";
      dummy = CODE_VER;
      TruncateString(dummy, 8);
      outFile->write( dummy.c_str(), 8 );
      (*outFile) << "simdate ";
      dummy = buffer;
      TruncateString(dummy, 8);
      outFile->write(dummy.c_str(), 8 );
      (*outFile) << "endgmv  ";
      (*outFile).flush();
    }
  }

  void SimOutputGMV::TruncateString(std::string & str,
                                    UInt maxLen)
  {
    UInt len = str.length();
    std::ostringstream strBuffer;

    if(len > maxLen)
    {
      strBuffer.clear();
      strBuffer << str;
      len = maxLen;
    }
    else
    { 
      strBuffer.clear();
      strBuffer << str;
      for(UInt i=0; i<maxLen-len; i++)
        strBuffer << '\0';
    }

    str = strBuffer.str();
    if(ascii_)
      str.resize(len);
    else
      str.resize(maxLen);
  }

  
}
