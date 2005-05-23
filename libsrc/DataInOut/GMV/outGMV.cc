#include <iostream>
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

  ENTER_FCN( "WriteResultsGMV :: WriteResultsGMV" );

  Integer nameLength = std::strlen(filename);
  namedir_ = new Char[ nameLength + 4 + 1 ];
  std::strcpy( namedir_, filename );
  std::strcat( namedir_, "_gmv" );

  Char *command = new Char[ nameLength + 4 + 9 + 1 ];
  std::sprintf( command, "mkdir -p %s", namedir_ );
  std::system( command );
  delete[] command;

  currStep_ = -1;
  lastStep_ = -1;
  lastTime_ = 0.0;
  currTime_ = 0.0;

  firstGridWritten_ = FALSE;
  output = NULL;
  ptgrid = NULL;

  ascii_ = TRUE;
  fixedgrid_ = TRUE;

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
        output->write( (char*)&lastStep_, sizeof(Integer) );
        (*output) << "endgmv  ";
      }
      delete output;
    }

    delete [] namedir_;
  }


  // **************
  //   WriteNodes
  // **************
  void WriteResultsGMV::WriteNodes() {

    // write keyword
    (*output) << "nodev   ";

    //get and write number of nodes on the level
    Integer numnodes=ptgrid->GetNumNodes();
    if (ascii_)
      (*output) << numnodes << std::endl;
    else
      output->write((char*)&numnodes,sizeof(Integer));

    Integer dim=ptgrid->GetDim();

    //get and write coodinates of nodes
    Integer i;

    if ( dim == 2 ) {

      Point<2> point;

      // write x,y,z-coordinate
      for ( i = 1; i <= numnodes; i++ ) {

        ptgrid->GetNodeCoordinate(point,i);

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

        ptgrid->GetNodeCoordinate(point,i);
        
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

    if (!ptgrid) {
      Error( "ptgrid is not initialized", __FILE__, __LINE__ );
    }

    // read information about number of elements 
    Integer numelem; 
    numelem=ptgrid->GetNumVolElems();

    if (ascii_)
      (*output) << numelem << std::endl;
    else
      output->write((char*)&numelem,sizeof(Integer));

    StdVector<Integer> connect;

    Integer dim=ptgrid->GetDim();
  
    Integer i;
    for ( i = 0; i < numelem; i++ ) {

      ptgrid->GetElemNodes(connect, i+1);

      if ( dim == 2 ) {
        switch ( connect.GetSize() ) {
          //          case 2:
          //            if (ascii_)
          //              (*output) << "line 2" << std::endl;
          //            else {
          //              (*output) << "line    ";
          //              Integer nn=2;
          //              output->write((char*)&nn,sizeof(Integer));
          //            }
        case 3: 
          if (ascii_)
            (*output) << "tri 3" << std::endl;
          else {
            (*output) << "tri     ";
            Integer nn=3;
            output->write((char*)&nn,sizeof(Integer));
          }
          break;
        case 4:
          if (ascii_)
            (*output) << "quad 4" << std::endl;
          else {
            (*output) << "quad    ";
            Integer nn=4;
            output->write((char*)&nn,sizeof(Integer));
          }
          break;
        case 6: 
          if (ascii_)
            (*output) << "6tri 6" << std::endl;
          else {
            (*output) << "6tri    ";
            Integer nn=6;
            output->write((char*)&nn,sizeof(Integer));
          }
          break;
        case 8:
          if (ascii_)
            (*output) << "8quad 8" << std::endl;
          else {
            (*output) << "8quad   ";
            Integer nn=8;
            output->write((char*)&nn,sizeof(Integer));
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
                Integer nn=4;
                output->write((char*)&nn,sizeof(Integer));
              }
              break;
            case 5:
              if (ascii_)
                (*output) << "ppyrmd5 5" << std::endl;
              else {
                (*output) << "ppyrmd5 ";
                Integer nn=5;
                output->write((char*)&nn,sizeof(Integer));
              }
              break;
            case 6:
              if (ascii_)
                (*output) << "pprism6 6" << std::endl;
              else {
                (*output) << "pprism6 ";
                Integer nn=6;
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
                Integer nn=15;
                // Integer nn=16;
                output->write((char*)&nn,sizeof(Integer));
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
                  Integer nn=8;
                  output->write((char*)&nn,sizeof(Integer));
                }
              break;
            case 20:
              if (ascii_)
                (*output) << "phex20 20" << std::endl;
              else 
                {
                  (*output) << "phex20  ";
                  Integer nn=20;
                  output->write((char*)&nn,sizeof(Integer));
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
          Integer j;
          for (j=0; j< connect.GetSize(); j++)
            (*output) << " " << connect[j] ;
          (*output) << std::endl;
        }
      else 
        {
          Integer * ptcon=connect.GetPointer();
          Integer len=connect.GetSize();
          output->write((char*)ptcon,len * sizeof(Integer));
        }

    }
  }


  // ******************
  //   WriteMaterials
  // ******************
  void WriteResultsGMV::WriteMaterials() {

    StdVector<Integer> regionID;
    StdVector<RegionIdType> subdoms;
    StdVector<Elem*> elemSD;
    std::string regionName;
    Integer aux;
    Char * str = NULL;
    if (! ascii_)
      str =new Char[8];
      
    ptgrid->GetVolRegionIds(subdoms);

    regionID.Resize(ptgrid->GetNumVolElems());

    if (ascii_)
      (*output) << "material " << subdoms.GetSize() << " 0" << std::endl;
    else {
      (*output) << "material";
      aux = subdoms.GetSize();
      output->write((char*)&aux,sizeof(Integer));
      aux = 0;
      output->write((char*)&aux,sizeof(Integer));
    }
  
    // loop over all subdomains
    for (Integer iSD=0; iSD<subdoms.GetSize(); iSD++) {
      regionName = ptgrid->RegionIdToName( subdoms[iSD] );
      if (ascii_)
        (*output) << regionName << std::endl;
      else {
        to8Char(regionName,str);
        (*output) << str;
      }
      

      ptgrid->GetVolElems(elemSD,subdoms[iSD]);

      // loop over all elemtns
      for (Integer iElem=0; iElem<elemSD.GetSize(); iElem++) 
        regionID[elemSD[iElem]->elemNum -1 ] = iSD+1;
    }

    // write for each element the according regionID
    for (Integer i=0; i<regionID.GetSize(); i++) {

      if (ascii_)
        (*output) << regionID[i] << " ";
      else {
        aux = regionID[i];
        output->write((char*)&aux,sizeof(Integer));
      }
    }
  
  
    if (ascii_)
      (*output) << std::endl;
  
    if (str)
      delete[] str;
  }


  // ******************************
  //   WriteNodeVariableTransient
  // ******************************
  void WriteResultsGMV::WriteNodeVariableTransient( const Vector<Double> var, 
                                                    const std::string name, 
                                                    const Integer dataType ) {
  
    if (ascii_) {
      (*output) << std::endl;
    }

    if (ascii_) {
      (*output) << name << " " << dataType << std::endl;
    }
    else {
      Char * str=new Char[8];  
      to8Char(name,str);
      (*output) << str;
      output->write((Char*)&dataType,sizeof(Integer));
      delete [] str;
    }

    if (ascii_) {
      Integer i;
      for (i=0; i<var.GetSize(); i++)
        (*output) << var[i] << " ";
    }
    else {
      Double * ptvar=var.GetPointer();
      Integer len=var.GetSize();
      output->write((char*)ptvar,len * sizeof(Double));
    }
  }


  // *****************************
  //   WriteNodeVariableHarmonic
  // *****************************
  void WriteResultsGMV::
  WriteNodeVariableHarmonic( const Vector<Complex> var, 
                             const std::string name, 
                             const Integer dataType,
                             const ComplexFormat outputFormat ) {
    Integer i;
    Double val;
  
    if (outputFormat == REAL_IMAG) {
  
      // **********************
      // REAL-IMAGINARY FORMAT
      // **********************
      
      // --- Real Part ---
      if (ascii_) (*output) << std::endl;
    
    
      if (ascii_)
        (*output) << name << "-Real " << dataType << std::endl;
      else {
        Char * str=new Char[8];  
        to8Char(name,str);
        str[7] = 'R';
        (*output) << str;
        output->write((Char*)&dataType,sizeof(Integer));
        delete [] str;
      }
    
      if (ascii_) {
        Integer i;
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
    
    
      if (ascii_)
        (*output) << name << "-Imag " << dataType << std::endl;
      else 
        {
          Char * str=new Char[8];  
          to8Char(name,str);
          str[7] = 'I';
          (*output) << str;
          output->write((Char*)&dataType,sizeof(Integer));
          delete [] str;
        }
    
      if (ascii_) 
        {
          Integer i;
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
    
    
      if (ascii_)
        (*output) << name << "-Amplitude " << dataType << std::endl;
      else 
        {
          Char * str=new Char[8];  
          to8Char(name,str);
          str[7] = 'A';
          (*output) << str;
          output->write((Char*)&dataType,sizeof(Integer));
          delete [] str;
        }
    
      if (ascii_) 
        {
          Integer i;
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
    
    
      if (ascii_)
        (*output) << name << "-Phase " << dataType << std::endl;
      else 
        {
          Char * str=new Char[8];  
          to8Char(name,str);
          str[7] = 'P';
          (*output) << str;
          output->write((Char*)&dataType,sizeof(Integer));
          delete [] str;
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
                                       const Integer dataType ) {

    (*output) << "velocity";
    if (ascii_)
      (*output) << " " << dataType << std::endl;
    else {
      output->write((char*)&dataType,sizeof(Integer));
    }

    if (ascii_) {
      Integer i,j;
      for (i=0; i<3; i++) {
        for (j=0; j<var[i].GetSize(); j++)
          (*output) << var[i][j] << " ";
        (*output) << std::endl;
      }
    }
    else {
      Integer i;
      for(i=0; i<3; i++) {
        Double * ptvar=var[i].GetPointer();
        Integer len=var[i].GetSize();
        output->write((char*)ptvar,len * sizeof(Double));
      }
    }
  }



  void WriteResultsGMV::WriteGrid()
  {

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
      if (ascii_)
        (*output) << "\nendgmv";
      else 
        (*output) << "endgmv  ";
      delete output;
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
  
    // Write Only if time time step has changed
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
      } else {
        WriteNodes();
        WriteCells();
        WriteMaterials();
      }
      if (ascii_)
        (*output) << std::endl << "variable" << std::endl;
      else
        (*output) << "variable";
    }
  }





  void WriteResultsGMV::
  WriteNodeSolutionTransient( const NodeStoreSol<Double> &sol, 
                              const Integer step, const Double time ) {
    ENTER_FCN( "WriteResultsGMV::WriteNodeSolutionTransient" );

    Integer iDof,i,j;
    Integer nrDofs = 0;
    Double help;
    Vector<Double> solhelp;
    std::string outString, errMsg;
    StdVector<SolutionType> solTypes;
  
    currTime_ = time;
    currStep_ = step;
    WriteGrid();  

    Integer type=1; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    // Get all solutiontypes
    sol.GetSolutionTypes(solTypes);

    // Iterate over all solutiontypes
    for (Integer iSol=0; iSol<sol.GetNumSolutions(); iSol++) {

      // GMV can not visualize tensor data
      if (sol.GetDof(solTypes[iSol]) > 3){
        errMsg  = "OutGMV::WriteNodeSolutionTransient: GMV can only ";
        errMsg += "visualize 3 dimensional data.\n Your solution has ";
        errMsg += Info->GenStr(sol.GetDof(solTypes[iSol]));
        errMsg += " degrees of freedom!";
        Error(errMsg.c_str(), __FILE__, __LINE__);
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
        WriteNodeVariableTransient(solhelp, outString , type);
      }
    }
  }


  void WriteResultsGMV::
  WriteElemSolutionTransient( const ElemStoreSol<Double> &data, 
                              const Integer step, const Double time ) {

    ENTER_FCN ( "WriteResultsGMV::WriteElemSolutionTransient" );

    Integer type=0; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    currTime_ = time;
    currStep_ = step;
    std::string errMsg;
    Integer i = 0;
    Vector<Double> solhelp;
    StdVector<SolutionType> solType;

    WriteGrid(); 
    data.GetSolutionTypes(solType);

    // GMV can not visualize tensor data
    if (data.GetDof() > 3){
      errMsg  = "OutGMV::WriteElemSolutionTransient: GMV can only ";
      errMsg += "visualize 3 dimensional data.\n Your solution has ";
      errMsg += Info->GenStr(data.GetDof());
      errMsg += " degrees of freedom!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
 
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
                            const Integer step,
                            const Double frequency, 
                            const ComplexFormat format)
  {
    ENTER_FCN( "WriteResultsGMV::WriteNodeSolutionHarmonic" );

    Integer iDof,i,j;
    Integer nrDofs = 0;
    Double help;
    Vector<Complex> solhelp;
    std::string outString, errMsg;
    StdVector<SolutionType> solTypes;
  
    currTime_ = frequency;
    currStep_ = step;
    WriteGrid();  
  
    Integer type=1; // 0 - for cell 
    // 1 - for node
    // 2 - for face data
  
    // Get all solutiontypes
    sol.GetSolutionTypes(solTypes);
  
    // Iterate over all solutiontypes
    for (Integer iSol=0; iSol<sol.GetNumSolutions(); iSol++) {
    
      // GMV can not visualize tensor data
      if (sol.GetDof(solTypes[iSol]) > 3){
        errMsg  = "OutGMV::WirteNodeSolutionHarmonic: GMV can only ";
        errMsg += "visualize 3 dimensional data.\n Your solution has ";
        errMsg += Info->GenStr(sol.GetDof(solTypes[iSol]));
        errMsg += " degrees of freedom!";
        Error(errMsg.c_str(), __FILE__, __LINE__);
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
                            const Integer step,
                            const Double frequency,
                            const ComplexFormat format)
  {
    ENTER_FCN( "WriteResultsGMV::WriteElemSolutionHarmonic" );

    Integer type=0; // 0 - for cell 
    // 1 - for node
    // 2 - for face data

    currTime_ = frequency;
    currStep_ = step;
    std::string errMsg;
    Integer i = 0;
    Vector<Complex> solhelp;
    StdVector<SolutionType> solType;
 
    WriteGrid(); 
    sol.GetSolutionTypes(solType);

    // GMV can not visualize tensor data
    if (sol.GetDof() > 3){
      errMsg  = "OutGMV::WriteElemSolutionHarmonic: GMV can only ";
      errMsg += "visualize 3 dimensional data.\n Your solution has ";
      errMsg += Info->GenStr(sol.GetDof());
      errMsg += " degrees of freedom!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
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
      filename.append( Info->GenStr( num ) );
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
        output->write( (char*)&lastStep_, sizeof(Integer) );
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
      (*output) << "gmvinput" << "ieeei4r8";
    }
  }


  // ********
  //   Init
  // ********
  void WriteResultsGMV::Init( Grid *aptgrid ) {
    ptgrid = aptgrid;

    // Initialize history files
    InitHistoryFiles();

  }

  void WriteResultsGMV::to8Char(const std::string name, char * result)
  {
    std::string aux;
    Integer i;

    aux="        "; 
    if (name.size()> 8) {
      for (i=0; i<8; i++)
        aux[i]=name[i];
    }
    else {
      for (i=0; i<name.size(); i++)
        aux[i]=name[i];
    }
    

    strcpy(result,aux.c_str());

  }

} // end of namespace
