#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdio.h>
#include <complex>

#include "outGMV.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "General/environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {


  //*****************
    //   Constructor
    //*****************
    WriteResultsGMV::WriteResultsGMV( const Char *const filename)
      : WriteResults(filename) {

      ENTER_FCN( "WriteResultsGMV::WriteResultsGMV" );

      UInt nameLength = std::strlen(filename);
      namedir_ = new Char[ nameLength + 4 + 1 ]; 
      std::strcpy( namedir_, filename );
      std::strcat( namedir_, "_gmv" );

      Char *command = new Char[ nameLength + 4 + 9 + 1 ];
      std::sprintf( command, "mkdir -p %s", namedir_ );
      std::system( command );
      delete[] command;

      currStep_ = 0;
      lastStep_ = 0;
      lastTime_ = 0.0;
      currTime_ = 0.0;

      firstGridWritten_ = FALSE;
      output = NULL;

      ascii_ = TRUE;
      fixedgrid_ = TRUE;

      // Set allowed lenght of characters for names
      charOutSize_ = 32;
      strBuffer_ = new Char[charOutSize_];

      // Output format can be either ascii (default) or binary
      ascii_ = !params->IsSet( "binaryFormat", "gmv" );
 
      // Does the grid change over time, or can we use a fixed grid
      fixedgrid_ = params->IsSet( "fixedGrid", "gmv" );
    }


  // **********************
  //   Default Destructor
  // **********************
  WriteResultsGMV::~WriteResultsGMV() {

    ENTER_FCN( "WriteResultsGMV::~WriteResultsGMV" );

    // The last gmv-file is closed by the destructor
    if ( output != NULL ) {
      if ( PrintGridOnly == FALSE ) {
        if (ascii_) {
          (*output) << "\nendvars\n" ;
          (*output) << "probtime " << lastTime_ << std::endl;
          (*output) << "cycleno " << lastStep_ << std::endl;
          (*output) << "endgmv ";
        }
        else {
          (*output) << "endvars ";
          (*output) << "probtime";
          output->write( (char*)&lastTime_, sizeof(Double) );
          (*output) << "cycleno ";
          output->write( (char*)&lastStep_, sizeof(UInt) );
          (*output) << "endgmv  ";
        }
      }
      delete output;
    }
    
    delete [] namedir_;
    
    delete [] strBuffer_;
  }


  // **************
  //   WriteNodes
  // **************
  void WriteResultsGMV::WriteNodes() {

    ENTER_FCN( "WriteResultsGMV::WriteNodes" );

    // write keyword
    (*output) << "nodev   ";

    //get and write number of nodes on the level
    UInt numnodes=ptGrid_->GetNumNodes();
    if (ascii_)
      (*output) << numnodes << std::endl;
    else
      output->write((char*)&numnodes,sizeof(UInt));

    UInt dim=ptGrid_->GetDim();

    //get and write coodinates of nodes
    UInt i;

    if ( dim == 2 ) {

      Point<2> point;

      // write x,y,z-coordinate
      for ( i = 1; i <= numnodes; i++ ) {

        ptGrid_->GetNodeCoordinate(point,i);

        if (ascii_) {
          (*output) << " " << point[0] << " " << point[1] << " "
                    << 0 << std::endl;
        }
        else {
          Double z=0.;
          output->write((char*)&point[0],sizeof(Double));
          output->write((char*)&point[1],sizeof(Double));
          output->write((char*)&z,sizeof(Double));
        }
      }
    }
    else {

      Point<3> point;

      // write x,y,z-coordinate
      for ( i = 1; i <= numnodes; i++ ) {

        ptGrid_->GetNodeCoordinate(point,i);
        
        if (ascii_) {
          (*output) << " " << point[0] << " " << point[1] << " "
                    << point[2] << std::endl;
        }
        else {
          output->write((char*)&point[0],sizeof(Double));
          output->write((char*)&point[1],sizeof(Double));
          output->write((char*)&point[2],sizeof(Double));
        }
      }
    }
  }


  // **************
  //   WriteCells
  // **************
  void WriteResultsGMV::WriteCells() {

    ENTER_FCN( "WriteResultsGMV::WriteCells" );

    // write keyword
    (*output) << "cells   ";

    if (!ptGrid_) {
      Error( "ptGrid_ is not initialized", __FILE__, __LINE__ );
    }

    // read information about number of elements 
    UInt numelem =ptGrid_->GetNumElems();

    if (ascii_)
      (*output) << numelem << std::endl;
    else
      output->write((char*)&numelem,sizeof(UInt));

    UInt i, dim;
    const Elem * ptElem = NULL;

    for ( i = 0; i < numelem; i++ ) {

      ptElem = ptGrid_->GetElem(i+1);
      
      StdVector<UInt> const & connect = ptElem->connect;
      dim = ptElem->ptElem->GetDim();

      if ( dim == 1 )
        switch ( connect.GetSize() ) {
        case 2:
          if (ascii_)
            (*output) << "line 2" << std::endl;
          else {
            (*output) << "line    ";
            UInt nn=2;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        case 3:
          if (ascii_)
            (*output) << "3line 3" << std::endl;
          else {
            (*output) << "3line   ";
            UInt nn=3;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        default:
          Error("This type of element is not implemented",
                __FILE__, __LINE__);
          
        } 
      else if ( dim == 2 ) {
        switch ( connect.GetSize() ) {
          //
          //            }
        case 3: 
          if (ascii_)
            (*output) << "tri 3" << std::endl;
          else {
            (*output) << "tri     ";
            UInt nn=3;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        case 4:
          if (ascii_)
            (*output) << "quad 4" << std::endl;
          else {
            (*output) << "quad    ";
            UInt nn=4;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        case 6: 
          if (ascii_)
            (*output) << "6tri 6" << std::endl;
          else {
            (*output) << "6tri    ";
            UInt nn=6;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        case 8:
          if (ascii_)
            (*output) << "8quad 8" << std::endl;
          else {
            (*output) << "8quad   ";
            UInt nn=8;
            output->write((char*)&nn,sizeof(UInt));
          }
          break;
        default:
          Error("This type of element is not implemented",
                __FILE__, __LINE__);
        }
      }
      else
        {
          switch (connect.GetSize())
            {
            case 4:
              if (ascii_)
                (*output) << "tet 4" << std::endl;
              else {
                (*output) << "tet     ";
                UInt nn=4;
                output->write((char*)&nn,sizeof(UInt));
              }
              break;
            case 5:
              if (ascii_)
                (*output) << "ppyrmd5 5" << std::endl;
              else {
                (*output) << "ppyrmd5 ";
                UInt nn=5;
                output->write((char*)&nn,sizeof(UInt));
              }
              break;
            case 6:
              if (ascii_)
                (*output) << "pprism6 6" << std::endl;
              else {
                (*output) << "pprism6 ";
                UInt nn=6;
                output->write((char*)&nn,sizeof(UInt));
              }
              break;
            case 10:
              if (ascii_)
                (*output) << "ptet10 10" << std::endl;
              else {
                (*output) << "ptet10  ";
                Integer nn=10;
                output->write((char*)&nn,sizeof(Integer));
              }
              break;
            case 13:
              if (ascii_)
                (*output) << "ppyrmd13 13" << std::endl;
              else {
                (*output) << "ppyrmd13";
                Integer nn=13;
                output->write((char*)&nn,sizeof(Integer));
              }
              break;
            case 15:
              if (ascii_) {
                (*output) << "pprism15 15" << std::endl;
                // (*output) << "pprism15 16" << std::endl;
              }
              else {
                (*output) << "pprism15";
                UInt nn=15;
                // UInt nn=16;
                output->write((char*)&nn,sizeof(UInt));
              }

              // Note: Due to a bug up to the current version of GMV (v3.5)
              // we have to add a 16. node for the wedge element by simply
              // copying the last node
              // connect.Push_back(connect[connect.GetSize()-1]);
              break;

            case 8:
              if (ascii_)
                (*output) << "phex8 8" << std::endl;
              else 
                {
                  (*output) << "phex8   ";
                  UInt nn=8;
                  output->write((char*)&nn,sizeof(UInt));
                }
              break;
            case 20:
              if (ascii_)
                (*output) << "phex20 20" << std::endl;
              else 
                {
                  (*output) << "phex20  ";
                  UInt nn=20;
                  output->write((char*)&nn,sizeof(UInt));
                }
              break;
            default:
              std::cout << connect.GetSize() << std::endl;
              Error("This type of element is not implemented",
                    __FILE__, __LINE__);
            }
        }

      if (ascii_) 
        {
          UInt j;
          for (j=0; j< connect.GetSize(); j++)
            (*output) << " " << connect[j] ;
          (*output) << std::endl;
        }
      else 
        {
          UInt * ptcon=connect.GetPointer();
          UInt len=connect.GetSize();
          output->write((char*)ptcon,len * sizeof(UInt));
        }

    }
  }


  // ******************
  //   WriteMaterials
  // ******************
  void WriteResultsGMV::WriteMaterials() {
    
    ENTER_FCN( "WriteResultsGMV::WriteMaterials");

    StdVector<UInt> regionID;
    StdVector<RegionIdType> subdoms;
    StdVector<Elem*> elemSD;
    std::string regionName;
    UInt aux;
      
    ptGrid_->GetRegionIds(subdoms);

    regionID.Resize(ptGrid_->GetNumElems());

    if (ascii_)
      (*output) << "material " << subdoms.GetSize() << " 0" << std::endl;
    else {
      (*output) << "material";
      aux = subdoms.GetSize();
      output->write((char*)&aux,sizeof(UInt));
      aux = 0;
      output->write((char*)&aux,sizeof(UInt));
    }
  
    // loop over all subdomains
    for (UInt iSD=0; iSD<subdoms.GetSize(); iSD++) {
      regionName = ptGrid_->RegionIdToName( subdoms[iSD] );

      TruncateString(regionName,"",strBuffer_);

      if (ascii_) {
        (*output) << strBuffer_ << std::endl;
      } else {
        output->write(strBuffer_,charOutSize_*sizeof(char));
      }
        
      

      ptGrid_->GetElems(elemSD,subdoms[iSD]);

      // loop over all elemtns
      for (UInt iElem=0; iElem<elemSD.GetSize(); iElem++) {
        regionID[elemSD[iElem]->elemNum -1 ] = iSD+1;
      }

    }

    // write for each element the according regionID
    for (UInt i=0; i<regionID.GetSize(); i++) {

      if (ascii_)
        (*output) << regionID[i] << " ";
      else {
        aux = regionID[i];
        output->write((char*)&aux,sizeof(UInt));
      }
    }
  
  
    if (ascii_)
      (*output) << std::endl;
  
  }


  // **************************************
  //   Write named entities (nodes, names)
  // **************************************
  void WriteResultsGMV::WriteNamedEntities() {

    ENTER_FCN( "WriteResultsGMV::WriteNamedEntities" );

    StdVector<std::string> nodeNames, elemNames;
    StdVector<UInt> nodeNumbers, elemNumbers;
    UInt entType = 0;
    UInt aux = 0;

    ptGrid_->GetListNodeNames(nodeNames);
    ptGrid_->GetListNodeNames(elemNames);
    
    // Check if there are any named nodes / elems in the grid
    if( nodeNames.GetSize() == 0 &&
        elemNames.GetSize() == 0 ) {
      return;
    }
    
    // Begin group section
    if (ascii_) {
      (*output) << "groups" << std::endl;
    } else {
      (*output) << "groups  ";
    }
    
    // write names nodes
    entType = 1;
      
    for( UInt i = 0; i < nodeNames.GetSize(); i++ ) {
      
      // get nodes
      ptGrid_->GetNodesByName(nodeNumbers, nodeNames[i]);

      // write name and number of nodes

      TruncateString(nodeNames[i],"",strBuffer_);
      if (ascii_) {
        (*output) << strBuffer_ <<  " " << entType << " " 
                  << nodeNumbers.GetSize();
      } else {
        output->write(strBuffer_,charOutSize_*sizeof(char));
        output->write((char*)&entType,sizeof(UInt));
        aux = nodeNumbers.GetSize();
        output->write((char*)&aux,sizeof(UInt));
      }

      // write node numbers itself
      
      if (ascii_) {
        for( UInt iNode = 0; iNode < nodeNumbers.GetSize(); iNode++) {
          (*output) << " " << nodeNumbers[iNode];
        }
      } else {
        for( UInt iNode = 0; iNode < nodeNumbers.GetSize(); iNode++) {
          aux = nodeNumbers[iNode];
          output->write((char*)&aux,sizeof(UInt));
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
  


  // ******************************
  //   WriteNodeVariableTransient
  // ******************************
  void WriteResultsGMV::WriteNodeVariableTransient( const Vector<Double> & var,
                                                    const std::string name,
                                                    const UInt dataType ) {

    ENTER_FCN("WriteResultsGMV::WriteNodeVariableTransient");
    
    if (ascii_) {
      (*output) << std::endl;
    }

    TruncateString(name,"",strBuffer_);
    if (ascii_) {
      (*output) << strBuffer_ << " " << dataType << std::endl;
    }
    else {

      output->write(strBuffer_,charOutSize_*sizeof(char));
      output->write((Char*)&dataType,sizeof(UInt));
    }

    if (ascii_) {
      UInt i;
      for (i=0; i<var.GetSize(); i++)
        (*output) << var[i] << " ";
    }
    else {
      Double * ptvar=var.GetPointer();
      UInt len=var.GetSize();
      output->write((char*)ptvar,len * sizeof(Double));
    }
  }


  // *****************************
  //   WriteNodeVariableHarmonic
  // *****************************
  void WriteResultsGMV::
  WriteNodeVariableHarmonic( const Vector<Complex> & var, 
                             const std::string name, 
                             const UInt dataType,
                             const ComplexFormat outputFormat ) {

    ENTER_FCN( "void WriteResultsGMV::WriteNodeVariableHarmonic" );
    
    UInt i;
    Double val;
    std::string auxName, suffix;
  
    if (outputFormat == REAL_IMAG) {
  
      // **********************
      // REAL-IMAGINARY FORMAT
      // **********************
      
      // --- Real Part ---
      if (ascii_) (*output) << std::endl;
    
      auxName = name;
      suffix = "-Real";
      TruncateString(auxName,suffix,strBuffer_);
      
      if (ascii_)
        (*output) << strBuffer_ << " "<< dataType << std::endl;
      else {
        output->write(strBuffer_,charOutSize_*sizeof(char));
        output->write((Char*)&dataType,sizeof(UInt));
      }
    
      if (ascii_) {
        UInt i;
        for (i=0; i<var.GetSize(); i++)
          (*output) << var[i].real() << " ";
      }
      else {
        for (i=0; i<var.GetSize(); i++){
          val = var[i].real();
          output->write((char*)(&val),sizeof(Double));
        }
      }

      // --- Imaginary Part ---
      if (ascii_) (*output) << std::endl;
      auxName = name;
      suffix = "-Imag";
      TruncateString(auxName,suffix,strBuffer_);
       
      if (ascii_)
        (*output) << strBuffer_ << " " << dataType << std::endl;
      else 
        {
          output->write(strBuffer_,charOutSize_*sizeof(char));
          output->write((Char*)&dataType,sizeof(UInt));
        }
    
      if (ascii_) 
        {
          UInt i;
          for (i=0; i<var.GetSize(); i++)
            (*output) << var[i].imag();
        }
      else {
        for (i=0; i<var.GetSize(); i++){
          val = var[i].imag();
          output->write((char*)(&val),sizeof(Double));
        }
      }
    
    }


    // **********************
    // AMPLITUDE-PHASE FORMAT
    // ********************** 
    else if ( outputFormat == AMPLITUDE_PHASE ) {
    
      // --- Amplitude ---
      if (ascii_) (*output) << std::endl;
      auxName = name;
      suffix = "-Ampl";
      TruncateString(auxName,suffix,strBuffer_);
      
      if (ascii_)
        (*output) << strBuffer_ << " " << dataType << std::endl;
      else 
        {
          output->write(strBuffer_,charOutSize_*sizeof(char));
          output->write((Char*)&dataType,sizeof(UInt));
        }
    
      if (ascii_) 
        {
          UInt i;
          for (i=0; i<var.GetSize(); i++)
            (*output) << std::abs(var[i]) << " ";
        }
      else 
        {
          for (i=0; i<var.GetSize(); i++){
            val = std::abs(var[i]);
            output->write((char*)(&val),sizeof(Double));
          }
        }

      // --- Phase ---
      if (ascii_) (*output) << std::endl;
      auxName = name;
      suffix = "-Phase";
      TruncateString(auxName,suffix,strBuffer_);
            
      if (ascii_)
        (*output) << strBuffer_ << " " << dataType << std::endl;
      else 
        {
          output->write(strBuffer_,charOutSize_*sizeof(char));
          output->write((Char*)&dataType,sizeof(UInt));
        }
    
      if (ascii_) 
        {
          for (i=0; i<var.GetSize(); i++) {
            if (abs(var[i].imag()) > 1e-16)
              (*output) << std::arg(var[i])*180/PI << " ";
            else
              (*output) << "0.0 ";
          }
        }
      else {
        for (i=0; i<var.GetSize(); i++){
          if (abs(var[i].imag()) > 1e-16)
            val = std::arg(var[i])*180/PI;
          else
            val = 0.0;
          output->write((char*)(&val),sizeof(Double));
        }
      }
    
    }
  }


  // *****************
  //   WriteVelocity
  // *****************
  void WriteResultsGMV::WriteVelocity( const Vector<Double> *var,
                                       const std::string name,
                                       const UInt dataType ) {

    ENTER_FCN( "WriteResultsGMV::WriteVelocity" );
    
    (*output) << "velocity";
    if (ascii_)
      (*output) << " " << dataType << std::endl;
    else {
      output->write((char*)&dataType,sizeof(UInt));
    }

    if (ascii_) {
      UInt i,j;
      for (i=0; i<3; i++) {
        for (j=0; j<var[i].GetSize(); j++)
          (*output) << var[i][j] << " ";
        (*output) << std::endl;
      }
    }
    else {
      UInt i;
      for(i=0; i<3; i++) {
        Double * ptvar=var[i].GetPointer();
        UInt len=var[i].GetSize();
        output->write((char*)ptvar,len * sizeof(Double));
      }
    }
  }



  void WriteResultsGMV::WriteGrid() {

    ENTER_FCN( "WriteResultsGMV::WriteVelocity" );
    
    // ---------------------------
    // Section for PrintGridOnly
    // ----------------------------
    if (PrintGridOnly == TRUE) {
      if(fixedgrid_ == TRUE) {
        OpenFile(-1);
      }
      else {
        OpenFile(0);
      }

      WriteNodes();
      WriteCells();
      WriteMaterials();
      WriteNamedEntities();
      
      if (ascii_)
        (*output) << "\nendgmv";
      else 
        (*output) << "endgmv  ";
      
      return;
    }

    // ---------------------------
    // Section for first 
    // fixedGrid step
    // ----------------------------
    if (fixedgrid_ == TRUE &&
        firstGridWritten_ == FALSE) {
      OpenFile(-1);
      WriteNodes();
      WriteCells();
      WriteMaterials();
      WriteNamedEntities();
      
      if (ascii_)
        (*output) << "\nendgmv";
      else 
        (*output) << "endgmv  ";
      firstGridWritten_ = TRUE;
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
        WriteNamedEntities();
      } else {
        WriteNodes();
        WriteCells();
        WriteMaterials();
        WriteNamedEntities();
      }
      if (ascii_)
        (*output) << std::endl << "variable" << std::endl;
      else
        (*output) << "variable";
    }
  }


  void WriteResultsGMV::
  WriteNodeSolutionTransient( const NodeStoreSol<Double> &sol, 
                              const UInt step, const Double time ) {

    ENTER_FCN( "WriteResultsGMV::WriteNodeSolutionTransient" );

    UInt iDof;
    Vector<Double> solhelp;
    std::string outString;
    StdVector<SolutionType> solTypes;
  
    currTime_ = time;
    currStep_ = step;
    WriteGrid();  

    UInt type=1; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    // Get all solutiontypes
    sol.GetSolutionTypes(solTypes);
    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<sol.GetNumSolutions(); iSol++) {

      // Iterate over all degrees of freedom
      for (iDof=0; iDof< sol.GetDof(solTypes[iSol]); iDof++) {
        if (sol.GetDof(solTypes[iSol]) > 1) {
          char nrStr[10];
          sprintf(nrStr,"%i",iDof+1);
          Enum2String( solTypes[iSol], outString );
          outString += nrStr;
        }
        else {
          Enum2String( solTypes[iSol], outString );
        }
      
        // Get the solution for one solutiontype and one dof
        sol.GetGlobalSolVectorSingleDof(solTypes[iSol],iDof,solhelp);

        // Write node solution to file
        WriteNodeVariableTransient(solhelp, outString , type);
      }
    }
  }


  void WriteResultsGMV::
  WriteElemSolutionTransient( const ElemStoreSol<Double> &data, 
                              const UInt step, const Double time ) {

    ENTER_FCN ( "WriteResultsGMV::WriteElemSolutionTransient" );

    UInt type=0; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    currTime_ = time;
    currStep_ = step;
    UInt i = 0;
    Vector<Double> solhelp;
    StdVector<SolutionType> solType;

    WriteGrid(); 
    data.GetSolutionTypes(solType);

    // Iterate over all degrees of freedom 
    for (i=0; i<data.GetDof(); i++) {
      char nrStr[10];
      sprintf(nrStr,"%i",i+1);
      std::string sumString;
      Enum2String( solType[0], sumString );
      sumString += nrStr;

      // Get global solution vector for on dof
      data.GetGlobalSolVectorSingleDof(i,solhelp);
   
      // Write global solution to file
      WriteNodeVariableTransient(solhelp, sumString , type);
    }
  }

  void WriteResultsGMV::
  WriteNodeSolutionHarmonic(const NodeStoreSol<Complex> & sol, 
                            const UInt step,
                            const Double frequency, 
                            const ComplexFormat format) {

    ENTER_FCN( "WriteResultsGMV::WriteNodeSolutionHarmonic" );

    UInt iDof;
    Vector<Complex> solhelp;
    std::string outString, errMsg;
    StdVector<SolutionType> solTypes;
  
    currTime_ = frequency;
    currStep_ = step;
    WriteGrid();  
  
    UInt type=1; // 0 - for cell 
    // 1 - for node
    // 2 - for face data
  
    // Get all solutiontypes
    sol.GetSolutionTypes(solTypes);
  
    // Iterate over all solutiontypes
    for (UInt iSol=0; iSol<sol.GetNumSolutions(); iSol++) {
    
      // GMV can not visualize tensor data
      if (sol.GetDof(solTypes[iSol]) > 3){
        (*warning) << "OutGMV::WriteElemSolutionTransient: GMV can only "
                   << "visualize 3 dimensional data.\n Your solution has "
                   << GenStr(sol.GetDof(solTypes[iSol]))
                   << " degrees of freedom!";
        Warning(__FILE__, __LINE__);
        return;
      }
        
      // Iterate over all degrees of freedom
      for (iDof=0; iDof< sol.GetDof(solTypes[iSol]); iDof++) {
        if (sol.GetDof(solTypes[iSol]) > 1) {
          char nrStr[10];
          sprintf(nrStr,"%i",iDof+1);
          Enum2String( solTypes[iSol], outString );
          outString += nrStr;
        }
        else {
          Enum2String( solTypes[iSol], outString );
        }
      
        // Get the solution for one solutiontype and one dof
        sol.GetGlobalSolVectorSingleDof(solTypes[iSol],iDof,solhelp);

        // Write node solution to file
        WriteNodeVariableHarmonic(solhelp, outString , type, format);
      }
    }
  
  }

  void WriteResultsGMV::
  WriteElemSolutionHarmonic(const ElemStoreSol<Complex> & sol, 
                            const UInt step,
                            const Double frequency,
                            const ComplexFormat format) {
    
    ENTER_FCN( "WriteResultsGMV::WriteElemSolutionHarmonic" );

    UInt type=0; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    currTime_ = frequency;
    currStep_ = step;
    std::string errMsg;
    UInt i = 0;
    Vector<Complex> solhelp;
    StdVector<SolutionType> solType;
 
    WriteGrid(); 
    sol.GetSolutionTypes(solType);

    // GMV can not visualize tensor data
    if (sol.GetDof() > 3){
      (*warning) << "OutGMV::WriteElemSolutionTransient: GMV can only "
                 << "visualize 3 dimensional data.\n Your solution has "
                 << GenStr(sol.GetDof())
                 << " degrees of freedom!";
      Warning(__FILE__, __LINE__);
      return;
    }
 
    // Iterate over all degrees of freedom 
    for (i=0; i<sol.GetDof(); i++) {
      char nrStr[10];
      sprintf(nrStr,"%i",i+1);
      std::string sumString;
      Enum2String( solType[0], sumString );
      sumString += nrStr;

      // Get global solution vector for on dof
      sol.GetGlobalSolVectorSingleDof(i,solhelp);
   
      // Write global solution to file
      WriteNodeVariableHarmonic(solhelp, sumString, type, format);
    }
  }


  // ***********
  //  OpenFile
  // ***********
  void WriteResultsGMV::OpenFile( const Integer num ) {

    ENTER_FCN( "WriteResultsGMV::OpenFile") ;

    std::string filename;

    // Generate basename for output file
    filename.append( namedir_ );
    filename.append( "/" );
    filename.append( namefile_ );

    // In the case of a fixed grid we write the grid description to a
    // separate file
    if ( num == -1 ) {
      filename.append( "_GRID.gmv" );
      nameGridFile_ = namefile_ + "_GRID.gmv";
    }

    // Normal output file
    else {
      filename.append( ".gmv" );
      if ( num < 10 ) filename.append( "000" );
      else if ( num < 100 ) filename.append( "00" );
      else if ( num < 1000 ) filename.append( "0" );
      else if ( num > 10000 ) {
        Info->Error( "Number of gmv file exceeds 9999!",
                     __FILE__, __LINE__ );
      }
      filename.append( GenStr( num ) );
    }

    // If file was already open, write end of variables
    // and the time and step number
    if (output && num > 1) {
      if (ascii_) {
        (*output) << "\nendvars\n" ;
        (*output) << "probtime " << lastTime_ << std::endl;
        (*output) << "cycleno " << lastStep_ << std::endl;
        (*output) << "endgmv";
      }
      else {
        (*output) << "endvars ";
        (*output) << "probtime";
        output->write( (char*)&lastTime_, sizeof(Double) );
        (*output) << "cycleno ";
        output->write( (char*)&lastStep_, sizeof(UInt) );
        (*output) << "endgmv  ";
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
      output = new std::ofstream(filename.c_str(),std::ofstream::binary);
    }
    if ( !output ) {
      (*error) << "Could not open file '" << filename
               << "' for writing GMV output";
      Error( __FILE__, __LINE__ );
    }

    // Write header, problem time and step number
    if (ascii_) {
      (*output) << "gmvinput ascii" << std::endl;
    }
    else {
      (*output) << "gmvinput" << "iecxi4r8";
    }
  }


  // ********
  //   Init
  // ********
  void WriteResultsGMV::Init( Grid *aptgrid ) {

    ENTER_FCN( "WriteResultsGMV::Init" );

    ptGrid_ = aptgrid;

    // Initialize history files
    InitHistoryFiles();

  }

  void WriteResultsGMV::TruncateString(const std::string name, 
                                       const std::string suffix, 
                                       char * result) {

    ENTER_FCN( "WriteResultsGMV::TruncateString" );

    UInt totalSize, insert_pos;

    // initalize an empty string with #00
    std::string aux(charOutSize_,'\0');

    // insert the name
    aux.insert(0,name,0,name.length());
    
    // if a 'suffix' was given, copy it at the end of the
    // buffer. if the size of  'name+suffix' is too long,
    // adapt the insert position, so that some characters of
    // 'name' get overwritten.
    
    if( suffix.length() > 0 ) {
      totalSize = name.length() + suffix.length();

      if( totalSize > charOutSize_ ) {
        insert_pos = charOutSize_ - suffix.length();
      } else {
        insert_pos = name.length();
      }
      aux.insert(insert_pos,suffix,0,suffix.length());
    }
    
    // copy the result string into the given char * buffer
    aux.copy(result,charOutSize_);
  }

} // end of namespace
