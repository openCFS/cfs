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

    ENTER_FCN( "SimOutputGMV::SimOutputGMV" );

    // Initialize variables
    formatName_ = "gmv";
    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    std::ostringstream strBuffer;

    dirName_ = "simoutput_gmv";
    fileName_ = fileName;

    try 
    {
      fs::create_directory( dirName_ );
    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }
    
    currStep_ = 0;
    lastStep_ = 0;
    lastTime_ = 0.0;
    currTime_ = 0.0;

    // Does the grid change over time, or can we use a fixed grid
    fixedgrid_ = true;
    printGridOnly_ = commandLine->GetPrintGrid();
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

    ENTER_FCN( "SimOutputGMV::~SimOutputGMV" );

  }

  void SimOutputGMV::BeginMultiSequenceStep( UInt step,
                                             AnalysisType type,
                                             UInt numSteps )
  {

  }
  
  void SimOutputGMV::RegisterResult( shared_ptr<BaseResult> sol,
                                     UInt saveBegin, UInt saveInc,
                                     UInt saveEnd )
  {
    ENTER_FCN( "SimOutputGMV::RegisterResult" );

    std::string type;
    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGMV) << "Registering output '" << actDof.resultName
                          << "' with saveBegin " << saveBegin
                          << ", saveInc " << saveInc 
                          << ", saveEnd " << saveEnd;
  }


  //! Begin single analysis step
  void SimOutputGMV::BeginStep( UInt stepNum, Double stepVal )
  {
    ENTER_FCN( "SimOutputGMV::BeginStep" );
    // std::cout << "BeginStep: stepNum: " << stepNum  << " stepVal: " << stepVal << std::endl;
    resultMap_.clear();
    
    currStep_ = stepNum;
    currTime_ = stepVal;
  }


  //! Add result to current step
  void SimOutputGMV::AddResult( shared_ptr<BaseResult> sol )
  {
    ENTER_FCN( "SimOutputGMV::AddResult" );
    std::string type;
    ResultInfo & actDof = *(sol->GetResultInfo());

    LOG_DBG(simOutputGMV) << "Adding result '" 
                          << actDof.resultName  << "'";
      
    resultMap_[sol->GetResultInfo()->resultName].Push_back(sol);
  }


  //! End single analysis step
  void SimOutputGMV::FinishStep( )
  {
    ENTER_FCN( "SimOutputGMV::FinishStep" );

    LOG_TRACE(simOutputGMV) << "Starting to finish Step";

    // Check, if any result at all is registered. 
    // If not leave!
    if( resultMap_.size() == 0 )
      return;

    

    WriteGrid();
    
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
  }

  //! End multisequence step
  void SimOutputGMV::FinishMultiSequenceStep( )
  {
    ENTER_FCN( "SimOutputGMV::FinishMultiSequenceStep" );
  }


  //! Finalize the output
  void SimOutputGMV::Finalize()
  {
    struct tm* time;
    time_t time_t;
    char buffer[64];

    memset(&time_t, sizeof(time_t), 0);
    memset(buffer, sizeof(buffer), 0);
    time = localtime(&time_t);
    strftime(buffer, sizeof(buffer), "%m/%d/%y", time);

    // The last gmv-file is closed by the destructor
    if ( output != NULL ) {
      if ( printGridOnly_ == false ) {
        if (ascii_) {
          (*output) << "probtime " << lastTime_ << std::endl;
          (*output) << "cycleno " << lastStep_ << std::endl;
          (*output) << "codename " << CODE_NAME << std::endl;
          (*output) << "codever " << CODE_VER << std::endl;
          (*output) << "simdate " << buffer << std::endl;
          (*output) << "endgmv ";
          (*output).flush();
        }
        else {
          (*output) << "probtime";
          output->write( (char*)&lastTime_, sizeof(Double) );
          (*output) << "cycleno ";
          output->write( (char*)&lastStep_, sizeof(UInt) );
          (*output) << "endgmv  ";
          (*output).flush();
        }
      }
      delete output;
    }
  }


  void SimOutputGMV::
  WriteNodeElemDataTrans( const Vector<Double> & var, 
                          const StdVector<std::string> & dofNames,
                          const std::string& name, 
                          ResultInfo::EntryType entryType,
                          ResultInfo::EntityUnknownType entityType,
                          Double time ) {
    ENTER_FCN ( "SimOutputGMV::WriteNodeElemDataTrans" );
    
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
    ENTER_FCN ( "SimOutputGMV::WriteNodeElemDataHarm" );

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
            vars1[i*numDofs+j] = var[i*numDofs*2+j].real();
            vars2[i*numDofs+j] = var[i*numDofs*2+j].imag();
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
            vars1[i*numDofs+j] = std::abs(var[i*numDofs*2+j]);
            vars2[i*numDofs+j] = CPhase(var[i*numDofs*2+j]);
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
  void SimOutputGMV::WriteNodes() {

    ENTER_FCN( "SimOutputGMV::WriteNodes" );

    // write keyword
    (*output) << "nodev   ";

    //get and write number of nodes on the level
    UInt numnodes = ptGrid_->GetNumNodes();
    if (ascii_)
      (*output) << numnodes << std::endl;
    else
      output->write((char*)&numnodes,sizeof(UInt));

    //get and write coodinates of nodes
    UInt i;

    // write x,y,z-coordinate
    for ( i = 1; i <= numnodes; i++) {
      Point p;
      
      ptGrid_->GetNodeCoordinate(p, i);  
        
      if (ascii_) {
        (*output) << " " << p[0] << " " << p[1] << " "
                  << p[2] << std::endl;
      }
      else {
        output->write((char*)&p[0],sizeof(Double));
        output->write((char*)&p[1],sizeof(Double));
        output->write((char*)&p[2],sizeof(Double));
      }
    }

  }


  // **************
  //   WriteCells
  // **************
  void SimOutputGMV::WriteCells() {

    ENTER_FCN( "SimOutputGMV::WriteCells" );

    // write keyword
    (*output) << "cells   ";

    // read information about number of elements 
    UInt numelem = ptGrid_->GetNumElems();

    if (ascii_)
      (*output) << numelem << std::endl;
    else
      output->write((char*)&numelem,sizeof(UInt));

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
    std::vector<std::string> regionNames;
    
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
        (*output) << gmvElemName << " " << numNodes << std::endl;
        UInt j;
        for (j=0; j < numNodes; j++)
          (*output) << " " << connect[j];
        (*output) << std::endl;
      }
      else {
        (*output) << gmvElemName;
        output->write((char*)&numNodes,sizeof(UInt));
        output->write((char*)&connect[0],numNodes * sizeof(UInt));
      }
    }

    // Write Regions -> Materials
    UInt aux = 0;
    if (ascii_)
      (*output) << "material " << numRegions << " 0" << std::endl;
    else {
      (*output) << "material";
      output->write((char*)&numRegions,sizeof(UInt));
      output->write((char*)&aux,sizeof(UInt));
    }

    for(i=0; i<numRegions; i++)
    {
      TruncateString(regionNames[i], charOutSize_);

      if (ascii_) {
        (*output) << regionNames[i] << std::endl;
      } else {
        output->write(&regionNames[i][0],charOutSize_*sizeof(char));
      }

    }

    // write for each element the according regionID
    for (i=0; i<numelem; i++) {

      if (ascii_)
        (*output) << elemRegions[i] << " ";
      else {
        output->write((char*)&elemRegions[i],sizeof(UInt));
      }
    }

    if (ascii_)
      (*output) << std::endl;
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
  void SimOutputGMV::WriteNamedEntities() {

    ENTER_FCN( "SimOutputGMV::WriteNamedEntities" );

    StdVector<std::string> nodeNames, elemNames;
    StdVector<UInt> nodeNumbers;
    StdVector<Elem*> elems;
    UInt entType = 0;
    UInt numNamedEntities = 0;

    ptGrid_->GetListNodeNames(nodeNames);

    ptGrid_->GetListElemNames(elemNames);

    // check if there are any named nodes / elems in the grid
    if( ( nodeNames.GetSize() == 0 ) && 
        ( elemNames.GetSize() == 0 ) )
    {
      return;
    }

    // Begin group section
    if (ascii_) {
      (*output) << "groups" << std::endl;
    } else {
      (*output) << "groups  ";
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
        (*output) << nodeNames[i] <<  " " << entType << " " 
                  << numNodes;
      } else {
        output->write(&nodeNames[i][0],charOutSize_*sizeof(char));
        output->write((char*)&entType,sizeof(UInt));
        output->write((char*)&numNodes,sizeof(UInt));
      }

      // write node numbers itself
      if (ascii_) {
        for( UInt iNode = 0; iNode < numNodes; iNode++) {
          (*output) << " " << nodeNumbers[iNode];
        }
      } else {
        for( UInt iNode = 0; iNode < numNodes; iNode++) {
          output->write((char*)&nodeNumbers[iNode],sizeof(UInt));
        }
      }
      
      if (ascii_) {
        (*output) << std::endl;
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
        (*output) << elemNames[i] <<  " " << entType << " " 
                  << numElems;
      } else {
        output->write(&elemNames[i][0],charOutSize_*sizeof(char));
        output->write((char*)&entType,sizeof(UInt));
        output->write((char*)&numElems,sizeof(UInt));
      }

      // write elem numbers itself

      if (ascii_) {
        for( UInt iElem = 0; iElem < numElems; iElem++) {
          (*output) << " " << elems[iElem]->elemNum;
        }
      } else {
        for( UInt iElem = 0; iElem < numElems; iElem++) {
          output->write((char*)&elems[iElem]->elemNum,sizeof(UInt));
        }
      }
      
      if (ascii_) {
        (*output) << std::endl;
      }
    }

    if (ascii_) {
      (*output) << "endgrp" << std::endl;
    } else {
      (*output) << "endgrp  ";
    }
      
    
  }
  

  void SimOutputGMV::WriteGrid() {

    ENTER_FCN( "SimOutputGMV::WriteGrid" );
    
    struct tm* time_tm;
    time_t time_t;
    char buffer[64];
    std::string dummy;

    time(&time_t);
    time_tm = localtime(&time_t);
    strftime(buffer, sizeof(buffer), "%m/%d/%y", time_tm);

    // ---------------------------
    // Section for PrintGridOnly
    // ----------------------------
    if (printGridOnly_ == true) {
      if(fixedgrid_ == true) {
        OpenFile(-1);
      }
      else {
        OpenFile(0);
      }

      WriteNodes();
      WriteCells();
      WriteNamedEntities();
      
      if (ascii_)
      {
        (*output) << "codename " << CODE_NAME << std::endl;
        (*output) << "codever " << CODE_VER << std::endl;
        (*output) << "simdate " << buffer << std::endl;
        (*output) << "endgmv" << std::endl;
        (*output).flush();
      }
      else 
      {
        (*output) << "codename";
        dummy = CODE_NAME;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "codever ";
        dummy = CODE_VER;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "simdate ";
        dummy = buffer;
        TruncateString(dummy, 8);
        output->write(dummy.c_str(), 8 );
        (*output) << "endgmv  ";
        (*output).flush();
      }
      return;
    }


    // ---------------------------
    // Section for first 
    // fixedGrid step
    // ----------------------------
    if (fixedgrid_ == true &&
        firstGridWritten_ == false) {
      OpenFile(-1);
      WriteNodes();
      WriteCells();
      WriteNamedEntities();
      
      if (ascii_)
      {
        (*output) << "codename " << CODE_NAME << std::endl;
        (*output) << "codever " << CODE_VER << std::endl;
        (*output) << "simdate " << buffer << std::endl;
        (*output) << "endgmv" << std::endl;
        (*output).flush();
      }
      else 
      {
        (*output) << "codename";
        dummy = CODE_NAME;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "codever ";
        dummy = CODE_VER;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "simdate ";
        dummy = buffer;
        TruncateString(dummy, 8);
        output->write(dummy.c_str(), 8 );
        (*output) << "endgmv  ";
        (*output).flush();
      }
        
      firstGridWritten_ = true;
        
      return;
    }
  
    // ---------------------------
    // Section for new time step
    // ----------------------------
    
    // Write Only if time step has changed
    // since last write of grid
    if (currStep_ != lastStep_) { 

      OpenFile(currStep_);

      if (fixedgrid_){
        if (ascii_)
          (*output) << "nodev fromfile \"" << nameGridFile_ 
                    <<"\""<< std::endl;
        else 
          (*output) << "nodev   fromfile\"" << nameGridFile_ <<"\"";
            
        if (ascii_)
          (*output) << "cells fromfile \"" << nameGridFile_ 
                    <<"\""<< std::endl;
        else 
          (*output) << "cells   fromfile\"" << nameGridFile_ <<"\"";
      
        if (ascii_)
          (*output) << "material fromfile \"" << nameGridFile_
                    <<"\""<< std::endl;
        else 
          (*output) << "materialfromfile\"" << nameGridFile_ <<"\"";
        
        // Since there exist no 'fromfile'-construct for groups, we have
        // to write it out all the time
        //        WriteNamedEntities();
      } else {
        WriteNodes();
        WriteCells();
        WriteNamedEntities();
      }
    }
  }


  // ***********
  //  OpenFile
  // ***********
  void SimOutputGMV::OpenFile( const Integer num ) {

    ENTER_FCN( "SimOutputGMV::OpenFile") ;

    std::string filename;
    std::ostringstream strBuffer;

    // Generate basename for output file
    filename.append( dirName_ );
    std::string pathsep = fs::path("/").native_directory_string();
    filename.append( pathsep );
    filename.append( fileName_ );

    // In the case of a fixed grid we write the grid description to a
    // separate file
    if ( num == -1 ) {
      filename.append( "_GRID.gmv" );
      nameGridFile_ = fileName_ + "_GRID.gmv";
    }

    // Normal output file
    else {
      strBuffer.clear();
      strBuffer << filename << ".gmv";
      if ( num < 10 ) strBuffer << "000";
      else if ( num < 100 ) strBuffer << "00";
      else if ( num < 1000 ) strBuffer << "0";
      else if ( num > 10000 ) {
        EXCEPTION("Number of gmv files exceeds 9999!");
      }
      strBuffer << num;
      filename = strBuffer.str();
    }

    struct tm* time_tm;
    time_t time_t;
    char buffer[64];
    std::string dummy;
    
    time(&time_t);
    time_tm = localtime(&time_t);
    strftime(buffer, sizeof(buffer), "%m/%d/%y", time_tm);

    // If file was already open, write end of variables
    // and the time and step number
    if (output && num > 1) {
      if (ascii_) {
        (*output) << "probtime " << lastTime_ << std::endl;
        (*output) << "cycleno " << lastStep_ << std::endl;
        (*output) << "codename " << CODE_NAME << std::endl;
        (*output) << "codever " << CODE_VER << std::endl;
        (*output) << "simdate " << buffer << std::endl;
        (*output) << "endgmv ";
        (*output).flush();
      }
      else {
        (*output) << "probtime";
        output->write( (char*)&lastTime_, sizeof(Double) );
        (*output) << "cycleno ";
        output->write( (char*)&lastStep_, sizeof(UInt) );
        (*output) << "codename";
        dummy = CODE_NAME;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "codever ";
        dummy = CODE_VER;
        TruncateString(dummy, 8);
        output->write( dummy.c_str(), 8 );
        (*output) << "simdate ";
        dummy = buffer;
        TruncateString(dummy, 8);
        output->write(dummy.c_str(), 8 );
        (*output) << "endgmv  ";
        (*output).flush();
      }
    }
    delete output;

    // Update time stepping information
    lastTime_ = currTime_;
    lastStep_ = currStep_;

    // check what kind of data for input
    std::string typedata;
  
    if (ascii_) {
      output = new std::ofstream(filename.c_str());
    }
    else {
      output = new std::ofstream(filename.c_str(), std::ofstream::binary);
    }
    if ( !output ) {
      EXCEPTION("Could not open file " << filename
                << " for writing GMV output");
    }

    // Write header, problem time and step number
    if (ascii_) {
      (*output) << "gmvinput ascii" << std::endl;
    }
    else {
      (*output) << "gmvinput" << "iecxi4r8";
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
