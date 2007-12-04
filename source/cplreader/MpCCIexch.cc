#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <math.h>
#include <unistd.h>

#include <boost/tokenizer.hpp>
#include <boost/filesystem/operations.hpp>
namespace bfs=boost::filesystem;

#include <cplreaderdefs.hh>
#include "settings.hh"

#ifndef CPLREADER_STANDALONE
#include "integlib.h"
#include "integlib/elemIntegr.hh"
#endif

#include "MpCCIexch.hh"


namespace CoupledField
{

  MpCCIexch::MpCCIexch(int argc, char *argv[], FileReader * ptFileReader)
  {
    ptFileReader_ = ptFileReader;

    if(!ptFileReader_)
      EXCEPTION("Invalid pointer to file reader!");

    ptFileReader_->Init();
  }
  
  MpCCIexch::~MpCCIexch()
  {
  }

    
  void MpCCIexch::PutExchangeGrid2MpCCI()
  {
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

    Settings& settings = Settings::Instance();
    uint32_t dim = settings.GetInt("dim");
    int nodeOffset = 0;
    int elemOffset = 0;
    int numPartitions = ptFileReader_->GetNumPartitions();
    std::vector<UInt>::iterator it, it2, it3, end;
    int localNodeNum = 1;
    std::set<UInt> regionNodeSet;
    UInt maxNumElemNodes = 0;
    std::string coordCfgFileName, coordDatFileName, connCfgFileName;
    std::string connDatFileName, typesCfgFileName, typesDatFileName;
    std::string fileName, importMeshCmd;
    std::ostringstream sstr;
    std::ofstream outFile;

    ClearMeshTempFiles(coordCfgFileName, coordDatFileName, connCfgFileName,
                       connDatFileName, typesCfgFileName, typesDatFileName,
                       importMeshCmd);
    
    // First read everything into internal buffers
    for (int actPart=0; actPart<numPartitions; actPart++)
    {
      ptFileReader_->ReadNodalCoords(NodalCoords_[actPart], actPart);
      ptFileReader_->ReadTopology(Topology_[actPart],
                                  numNodesPerElem_[actPart],
                                  elemTypes_[actPart],
                                  actPart);

      // Determine the maximum number of element nodes
      it = elemTypes_[actPart].begin();
      end = elemTypes_[actPart].end();
      for( ; it != end; it++ )
      {
        maxNumElemNodes = NUM_ELEM_NODES[*it] > maxNumElemNodes ? NUM_ELEM_NODES[*it] : maxNumElemNodes;
      }
    }      

    for (int actPart=0; actPart<numPartitions; actPart++)
    {
      // Put all nodes in a partition into a set to get an ordered list
      it = Topology_[actPart].begin ();
      end = Topology_[actPart].end ();
      std::set<UInt>::iterator regSetIt, regSetEnd;
      regionNodeSet.clear();
      regionNodeSet.insert(it, end);
      regSetIt = regionNodeSet.begin(); 
      regSetEnd = regionNodeSet.end(); 

      for( ; regSetIt != regSetEnd; regSetIt++ )
      {
        nodesInPartition_[actPart][*regSetIt] = localNodeNum;
        localNodeNum++;
      }
      int nNodes = nodesInPartition_[actPart].size();
      int nElems = ptFileReader_->GetNumElems(actPart);

      // Write node coordinates
      outFile.open(coordDatFileName.c_str(),
                   std::ofstream::out | std::ofstream::app);

      regSetIt = regionNodeSet.begin(); 
      regSetEnd = regionNodeSet.end(); 
      for( ; regSetIt != regSetEnd; regSetIt++ )
      {
        int idx = *regSetIt - 1;
        outFile << NodalCoords_[actPart][idx*3+0] << " "
                << NodalCoords_[actPart][idx*3+1] << " "
                << NodalCoords_[actPart][idx*3+2] << std::endl;
      }
      outFile.close();

      // Replace old node numbers in connectivity with new ones 
      it = elemTypes_[actPart].begin();
      end = elemTypes_[actPart].end();
      it2 = Topology_[actPart].begin();
      std::vector<UInt> newTopo;
      newTopo.resize(maxNumElemNodes * elemTypes_[actPart].size());
      it3 = newTopo.begin();
      for( ; it != end; it++ )
      {
        int numElemNodes = NUM_ELEM_NODES[*it];

        for(int i=0; i<numElemNodes; i++)
        {
          *(it3+i) = nodesInPartition_[actPart][*(it2+i)];
        }

        it3 += maxNumElemNodes;
        it2 += numElemNodes;
      }        

      // Write out element connectivity
      outFile.open(connDatFileName.c_str(),
                   std::ofstream::out | std::ofstream::app);
      typedef std::ostream_iterator<int> int_os_iter;
      std::copy (newTopo.begin (), newTopo.end (), int_os_iter (outFile, "\n"));
      outFile.close();


      // Write out element types
      outFile.open(typesDatFileName.c_str(),
                   std::ofstream::out | std::ofstream::app);
      std::copy (elemTypes_[actPart].begin (), elemTypes_[actPart].end (), int_os_iter (outFile, "\n"));
      outFile.close();


      // Write out nodes of region
      sstr.str("");
      sstr << settings.GetString("name") << "_data/RegionNodes_partition" << (actPart+1) << ".dat";
      fileName = sstr.str();
      outFile.open(fileName.c_str(),
                   std::ofstream::out | std::ofstream::trunc);

      regSetIt = regionNodeSet.begin(); 
      regSetEnd = regionNodeSet.end(); 
      for( ; regSetIt != regSetEnd; regSetIt++ )
      {
        int nodeNum = *regSetIt;

        outFile << nodesInPartition_[actPart][nodeNum] << std::endl;
      }
      regionNodeSet.clear();
      outFile.close();

      // Write out elements of region
      sstr.str("");
      sstr << settings.GetString("name") << "_data/RegionElems_partition" << (actPart+1) << ".dat";
      fileName = sstr.str();
      outFile.open(fileName.c_str(),
                   std::ofstream::out | std::ofstream::trunc);

      for(int i=1; i<=nElems; i++)
      {
        outFile << (elemOffset + i) << std::endl;
      }
      outFile.close();

      // Write out dimension of region      
      sstr.str("");
      sstr << settings.GetString("name") << "_data/RegionDim_partition" << (actPart+1) << ".dat";
      fileName = sstr.str();
      outFile.open(fileName.c_str(),
                   std::ofstream::out | std::ofstream::trunc);
      outFile << dim << std::endl;
      outFile.close();

      elemOffset += nElems;
      nodeOffset += nNodes;
    }
    
    outFile.open(coordCfgFileName.c_str(),
                 std::ofstream::out | std::ofstream::trunc);
    outFile << "PATH /Mesh/Nodes/Coordinates" << std::endl;
    outFile << "INPUT-CLASS TEXTFP" << std::endl;
    outFile << "INPUT-SIZE 64" << std::endl;
    outFile << "RANK 2" << std::endl;
    outFile << "DIMENSION-SIZES " << nodeOffset << " " << 3 << std::endl;
    outFile << "OUTPUT-CLASS FP" << std::endl;;
    outFile << "OUTPUT-SIZE 64" << std::endl;;
    outFile << "OUTPUT-ARCHITECTURE IEEE" << std::endl;;
    outFile << "OUTPUT-BYTE-ORDER LE" << std::endl;
    outFile << "CHUNKED-DIMENSION-SIZES "
            << (nodeOffset < 1000 ? nodeOffset : 1000)
            << " 3" << std::endl;
    outFile << "COMPRESSION-TYPE GZIP" << std::endl;
    outFile << "COMPRESSION-PARAM 9" << std::endl;

    outFile.close();

    outFile.open(connCfgFileName.c_str(),
                 std::ofstream::out | std::ofstream::trunc);
    outFile << "PATH /Mesh/Elements/Connectivity" << std::endl;
    outFile << "INPUT-CLASS TEXTIN" << std::endl;
    outFile << "INPUT-SIZE 32" << std::endl;
    outFile << "RANK 2" << std::endl;
    outFile << "DIMENSION-SIZES " << elemOffset << " "
            << maxNumElemNodes << std::endl;
    outFile << "OUTPUT-CLASS UIN" << std::endl;
    outFile << "OUTPUT-SIZE 32" << std::endl;
    outFile << "OUTPUT-ARCHITECTURE STD" << std::endl;
    outFile << "OUTPUT-BYTE-ORDER LE" << std::endl;
    outFile << "CHUNKED-DIMENSION-SIZES " 
            << (elemOffset < 1000 ? elemOffset : 1000)
            << " " << maxNumElemNodes << std::endl;
    outFile << "COMPRESSION-TYPE GZIP" << std::endl;
    outFile << "COMPRESSION-PARAM 9" << std::endl;
    outFile.close();

    outFile.open(typesCfgFileName.c_str(),
                 std::ofstream::out | std::ofstream::trunc);
    outFile << "PATH /Mesh/Elements/Types" << std::endl;
    outFile << "INPUT-CLASS TEXTIN" << std::endl;
    outFile << "INPUT-SIZE 32" << std::endl;
    outFile << "RANK 1" << std::endl;
    outFile << "DIMENSION-SIZES " << elemOffset << std::endl;
    outFile << "OUTPUT-CLASS UIN" << std::endl;
    outFile << "OUTPUT-SIZE 32" << std::endl;
    outFile << "OUTPUT-ARCHITECTURE STD" << std::endl;
    outFile << "OUTPUT-BYTE-ORDER LE" << std::endl;
    outFile << "CHUNKED-DIMENSION-SIZES "
            << (elemOffset < 1000 ? elemOffset : 1000) << std::endl;
    outFile << "COMPRESSION-TYPE GZIP" << std::endl;
    outFile << "COMPRESSION-PARAM 9" << std::endl;
    outFile.close();

    if(std::system(importMeshCmd.c_str()))
    {
      EXCEPTION("Import of mesh into HDF5 file failed!");
    }

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

  }

  void MpCCIexch::Couple()
  {
    Settings& settings = Settings::Instance();
    bool calcSrc = settings.GetInt("calcSrc");
    bool good = true;
    UInt counter = 0;
    UInt numFiles = ptFileReader_->GetNumFiles();
    UInt numPartitions = ptFileReader_->GetNumPartitions();
    std::ostringstream sstr;
    std::string writeResultCmd;
    std::string writeAuxCmd;        
    std::string stepNumsFileName;
    std::string stepValuesFileName;
    std::string resultScriptFileName;
    std::string resultDatFileName;
    std::vector<double> acouSrc;
    std::vector<double> velocity;
    std::vector<double> flowdata;
    std::vector<double> tmpDatLHSrc;
    std::vector<double> tmpDatFluidVel;
    std::vector<double> tmpDatFluidPres;
    UInt nodeOffset;
    UInt actPart;
    double stepVal = counter*settings.GetDouble("timeStep");
    UInt numNodesInPartition;
    std::ofstream outFile;
    std::vector<std::string> outputFields;
    
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Time Data Conversion Phase!" << std::endl;
    std::cout << "                        ";
    std::cout << "Number of transient files: " << numFiles << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

    ClearDatasetTempFiles(stepNumsFileName, stepValuesFileName,
                          resultScriptFileName, resultDatFileName,
                          writeResultCmd, writeAuxCmd);

    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");
    Tok t(settings.GetString("outputfields"), sep);
    Tok::iterator it, end;
    it = t.begin();
    end = t.end();
    
    outputFields.resize(3);
    for( ; it != end; it++)
    {
      if(*it == "acouRhsLoad")
        outputFields[0] = "acouRhsLoad"; 
      if(*it == "fluidMechVelocity")
        outputFields[1] = "fluidMechVelocity"; 
      if(*it == "fluidMechPressure")
        outputFields[2] = "fluidMechPressure"; 
    }
    
    
    while ( counter < numFiles ) 
    {
      nodeOffset = 0;
      stepVal = counter*settings.GetDouble("timeStep");
      std::cout << "----------------------------------------"
                << "----------------------------------------"
                << std::endl;
      std::cout << "                        "
                << "Step Number: " << (counter+1) << "; ";
      std::cout << " Time Step: " << stepVal << std::endl;
      std::cout << "----------------------------------------"
                << "----------------------------------------"
                << std::endl;


      for (actPart = 0; actPart<numPartitions; actPart++)
      {
        int numNodesPart = ptFileReader_->GetNumNodes(actPart);
        numNodesInPartition = nodesInPartition_[actPart].size();
        //              int numResults = ptFileReader_->GetNumResults();

        // Read acou source, velocity and pressure. 
        ptFileReader_->ReadNodalValues(flowdata,
                                       actPart,
                                       counter);


        // by convention the flowdata array has to have
        // the following structure:
        // flowdata [0]      [1]  [2]  [3]  [4]   [5]   [6]   [7]      [8] ... 
        //          acousrc1 vx1  vy1  vz1  usr11 usr21 usr31 acousrc2 vx2 ...

#ifndef CPLREADER_STANDALONE
        if(calcSrc)
        {
          for (int inode=0; inode<numNodesPart; inode++)
          {
            flowdata[inode*7] = 0.0;
          }

          CalculateAcouSrcs(actPart, flowdata);
        }
#endif

        // Prepare vectors with Lighthill source, velocity and pressure. 
        std::map<UInt, UInt>::iterator it, end;
        it = nodesInPartition_[actPart].begin();
        end = nodesInPartition_[actPart].end();

        tmpDatLHSrc.resize(numNodesInPartition);
        tmpDatFluidVel.resize(numNodesInPartition*3);
        tmpDatFluidPres.resize(numNodesInPartition);
        int newNodeIdx = 0;
        for( ; it != end; it++ ) {
          int oldNodeIdx = it->first-1;
          tmpDatLHSrc[newNodeIdx] = flowdata[oldNodeIdx*7+0];

          tmpDatFluidVel[newNodeIdx*3+0] = flowdata[oldNodeIdx*7+1];
          tmpDatFluidVel[newNodeIdx*3+1] = flowdata[oldNodeIdx*7+2];
          tmpDatFluidVel[newNodeIdx*3+2] = flowdata[oldNodeIdx*7+3];

          tmpDatFluidPres[newNodeIdx] = flowdata[oldNodeIdx*7+4];

          newNodeIdx++;                           
        }

        typedef std::ostream_iterator<double> double_os_iter;

        // Write Lighthill source term              
        if(outputFields[0] == "acouRhsLoad") 
        {
          outFile.open(resultScriptFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);
          outFile << "RESULT_TYPE=Mesh" << std::endl;
          outFile << "MULTISTEP=1" << std::endl;
          outFile << "ANALYSIS_TYPE=transient" << std::endl;
          outFile << "RESULT_NAME=acouRhsLoad" << std::endl;
          outFile << "RESULT_REGION=partition" << (actPart+1) << std::endl;
          outFile << "DESC_DOFNAMES=-" << std::endl;
          outFile << "DESC_DEFINEDON=1" << std::endl;
          outFile << "DESC_ENTRYTYPE=1" << std::endl;
          outFile << "DESC_NUMDOFS=1" << std::endl;
          outFile << "DESC_REGIONS=partition" << (actPart+1) << std::endl;
          outFile << "DESC_UNIT=\"kg m^-3 s^-2\"" << std::endl;
          outFile << "STEP_NUM=" << (counter+1) << std::endl;
          outFile << "STEP_VALUE=" << stepVal << std::endl;
          outFile << "RESULT_ON=Nodes" << std::endl;
          outFile << "COMPLEX_TYPE=Real" << std::endl;
          outFile << "RESULT_DATAFILE=Result.dat" << std::endl;
          outFile << "RESULT_RANK=1" << std::endl;
          outFile << "RESULT_DIMENSION_SIZES=" << numNodesInPartition << std::endl;
          outFile.close(); outFile.clear();

          outFile.open(resultDatFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);

          std::copy (tmpDatLHSrc.begin (), tmpDatLHSrc.end (), double_os_iter (outFile, "\n"));
          outFile.close(); outFile.clear();

          if(std::system(writeResultCmd.c_str()))
          {
            EXCEPTION("Import of Lighthill dataset into HDF5 file failed!");
          }
        }

        // Write fluid velocities
        if(outputFields[1] == "fluidMechVelocity") 
        {
          outFile.open(resultScriptFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);
          outFile << "RESULT_TYPE=Mesh" << std::endl;
          outFile << "MULTISTEP=1" << std::endl;
          outFile << "ANALYSIS_TYPE=transient" << std::endl;
          outFile << "RESULT_NAME=fluidMechVelocity" << std::endl;
          outFile << "RESULT_REGION=partition" << (actPart+1) << std::endl;
          outFile << "DESC_DOFNAMES=\"x;y;z\"" << std::endl;
          outFile << "DESC_DEFINEDON=1" << std::endl;
          outFile << "DESC_ENTRYTYPE=3" << std::endl;
          outFile << "DESC_NUMDOFS=3" << std::endl;
          outFile << "DESC_REGIONS=partition" << (actPart+1) << std::endl;
          outFile << "DESC_UNIT=\"m s^-1\"" << std::endl;
          outFile << "STEP_NUM=" << (counter+1) << std::endl;
          outFile << "STEP_VALUE=" << stepVal << std::endl;
          outFile << "RESULT_ON=Nodes" << std::endl;
          outFile << "COMPLEX_TYPE=Real" << std::endl;
          outFile << "RESULT_DATAFILE=Result.dat" << std::endl;
          outFile << "RESULT_RANK=2" << std::endl;
          outFile << "RESULT_DIMENSION_SIZES=\"" << numNodesInPartition 
          << " 3\"" << std::endl;
          outFile.close(); outFile.clear();

          outFile.open(resultDatFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);

          std::copy (tmpDatFluidVel.begin (), tmpDatFluidVel.end (),
              double_os_iter (outFile, "\n"));
          outFile.close(); outFile.clear();

          if(std::system(writeResultCmd.c_str()))
          {
            EXCEPTION("Import of velocity into HDF5 file failed");
          }
        }

        // Write fluid pressure
        if(outputFields[2] == "fluidMechPressure") 
        {
          outFile.open(resultScriptFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);
          outFile << "RESULT_TYPE=Mesh" << std::endl;
          outFile << "MULTISTEP=1" << std::endl;
          outFile << "ANALYSIS_TYPE=transient" << std::endl;
          outFile << "RESULT_NAME=fluidMechPressure" << std::endl;
          outFile << "RESULT_REGION=partition" << (actPart+1) << std::endl;
          outFile << "DESC_DOFNAMES=-" << std::endl;
          outFile << "DESC_DEFINEDON=1" << std::endl;
          outFile << "DESC_ENTRYTYPE=1" << std::endl;
          outFile << "DESC_NUMDOFS=1" << std::endl;
          outFile << "DESC_REGIONS=partition" << (actPart+1) << std::endl;
          outFile << "DESC_UNIT=Pa" << std::endl;
          outFile << "STEP_NUM=" << (counter+1) << std::endl;
          outFile << "STEP_VALUE=" << stepVal << std::endl;
          outFile << "RESULT_ON=Nodes" << std::endl;
          outFile << "COMPLEX_TYPE=Real" << std::endl;
          outFile << "RESULT_DATAFILE=Result.dat" << std::endl;
          outFile << "RESULT_RANK=1" << std::endl;
          outFile << "RESULT_DIMENSION_SIZES=" << numNodesInPartition << std::endl;
          outFile.close(); outFile.clear();

          outFile.open(resultDatFileName.c_str(),
              std::ofstream::out | std::ofstream::trunc);

          std::copy (tmpDatFluidPres.begin (), tmpDatFluidPres.end (),
              double_os_iter (outFile, "\n"));
          outFile.close(); outFile.clear();

          if(std::system(writeResultCmd.c_str()))
          {
            EXCEPTION("Import of pressure into HDF5 file failed");
          }
        }
        
        nodeOffset += numNodesPart;

      }//end of for

      outFile.open(stepNumsFileName.c_str(),
                   std::ofstream::out | std::ofstream::app);
      outFile << (counter+1) << std::endl;
      outFile.close(); outFile.clear();

      outFile.open(stepValuesFileName.c_str(),
                   std::ofstream::out | std::ofstream::app);
      outFile << stepVal << std::endl;
      outFile.close(); outFile.clear();

      counter++;
    }//end of while

    if(std::system(writeAuxCmd.c_str()))
    {
      EXCEPTION("Import of additional data into HDF5 file failed");
    }

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Time Data Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
  }

  void MpCCIexch::CalculateAcouSrcs(const int partitionIdx,
                                    std::vector<double>& flowdata)
  {
#ifndef CPLREADER_STANDALONE
    Settings& settings = Settings::Instance();
    int nElems = ptFileReader_->GetNumElems(partitionIdx);

    Matrix<Double> coordMatQuad;
    Matrix<Double> nodaldTijdxjQuad;
    Matrix<Double> nodalVelQuad;

    coordMatQuad.Resize(2, 4);
    nodaldTijdxjQuad.Resize(2, 4);
    nodalVelQuad.Resize(2, 4);
    Vector<Double> elemVecQuad;

    Matrix<Double> coordMatHex;
    Matrix<Double> nodaldTijdxjHex;
    Matrix<Double> nodalVelHex;

    coordMatHex.Resize(3, 8);
    nodaldTijdxjHex.Resize(3, 8);
    nodalVelHex.Resize(3, 8);
    Vector<Double> elemVecHex;

    int k=0;

    std::map<UInt, ElemIntegr *> ptElemIntegr;


    for( int i=0; i<nElems; i++)
    {
      switch(elemTypes_[partitionIdx][i])
      {
      case ET_LINE2:
        break;

      case ET_TRIA3:
        break;

      case ET_QUAD4:
        if(!ptElemIntegr[ET_QUAD4])
          ptElemIntegr[ET_QUAD4]=new ElemIntegr(ET_QUAD4);

        for( int n=0; n<numNodesPerElem_[partitionIdx][i]; n++)
        {
          UInt idx = (Topology_[partitionIdx][k+n]-1)*7;

          coordMatQuad[0][n] = NodalCoords_[partitionIdx][(Topology_[partitionIdx][k+n]-1)*3+0];
          coordMatQuad[1][n] = NodalCoords_[partitionIdx][(Topology_[partitionIdx][k+n]-1)*3+1];

          nodalVelQuad[0][n] = flowdata[idx+1];
          nodalVelQuad[1][n] = flowdata[idx+2];
        }

        try {
          ptElemIntegr[ET_QUAD4]->PerformIntegration( coordMatQuad, nodaldTijdxjQuad,
                                                      nodalVelQuad, elemVecQuad);
        } catch (CoupledField::Exception &ex)
        {
          std::cerr << "Warning: An Exception occurred during source term "
                    << "computation!" << std::endl;
          std::cerr << "Corner coords: " << coordMatQuad << std::endl;

          if(settings.GetInt("verbose"))
            std::cerr << ex.what();

          std::cerr << "Setting contribution to acousrc to zero!" << std::endl;
          elemVecQuad.Resize(4);
          elemVecQuad.Init(0.0);
        }

        for( int n=0; n<numNodesPerElem_[partitionIdx][i]; n++)
        {
          UInt idx = (Topology_[partitionIdx][k+n]-1)*7;

          flowdata[idx+0] += elemVecQuad[n];
        }
        break;

      case ET_QUAD8:
        break;

      case ET_TET4:
        break;

      case ET_WEDGE6:
        break;

      case ET_HEXA8:
        if(!ptElemIntegr[ET_HEXA8])
          ptElemIntegr[ET_HEXA8]=new ElemIntegr(ET_HEXA8);

        for( int n=0; n<numNodesPerElem_[partitionIdx][i]; n++)
        {
          int idx = (Topology_[partitionIdx][k+n]-1)*7;
          int idx2 = (Topology_[partitionIdx][k+n]-1)*3;

          coordMatHex[0][n] = NodalCoords_[partitionIdx][idx2+0];
          coordMatHex[1][n] = NodalCoords_[partitionIdx][idx2+1];
          coordMatHex[2][n] = NodalCoords_[partitionIdx][idx2+2];

          nodalVelHex[0][n] = flowdata[idx+1];
          nodalVelHex[1][n] = flowdata[idx+2];
          nodalVelHex[2][n] = flowdata[idx+3];
        }

        try {
          ptElemIntegr[ET_HEXA8]->PerformIntegration( coordMatHex, nodaldTijdxjHex,
                                                      nodalVelHex, elemVecHex);
        } catch (CoupledField::Exception &ex)
        {
          std::cerr << "Warning: An Exception occurred during source term "
                    << "computation!" << std::endl;
          std::cerr << "Corner coords: " << coordMatHex << std::endl;

          if(settings.GetInt("verbose"))
            std::cerr << ex.what();

          std::cerr << "Setting contribution to acousrc to zero!" << std::endl;
          elemVecHex.Resize(8);
          elemVecHex.Init(0.0);
        }

        for( int n=0; n<numNodesPerElem_[partitionIdx][i]; n++)
        {
          int idx = (Topology_[partitionIdx][k+n]-1)*7;

          flowdata[idx+0] += elemVecHex[n];
        }
        break;

      case ET_PYRA5:
        break;
      }

      k += numNodesPerElem_[partitionIdx][i];
    }
#endif
  }

  void MpCCIexch::ClearMeshTempFiles(std::string& coordCfgFileName,
                                     std::string& coordDatFileName,
                                     std::string& connCfgFileName,
                                     std::string& connDatFileName,
                                     std::string& typesCfgFileName,
                                     std::string& typesDatFileName,
                                     std::string& importMeshCmd)
  {
    Settings& settings = Settings::Instance();
    std::ostringstream sstr;

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Coordinates.cfg";
    coordCfgFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Coordinates.dat";
    coordDatFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Connectivity.cfg";
    connCfgFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Connectivity.dat";
    connDatFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Types.cfg";
    typesCfgFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Types.dat";
    typesDatFileName = sstr.str();
    bfs::remove(bfs::path(sstr.str()));
    
    sstr.str("");
    sstr << settings.GetString("name") << "_data/"
         << settings.GetString("name") << ".h5";
    bfs::remove(bfs::path(sstr.str()));

    sstr.str("");
    sstr << settings.GetString("baseExeDir") << "/hdf5_import_mesh.sh "
         << settings.GetString("name") << "_data "
         << settings.GetString("name") << ".h5";
    importMeshCmd = sstr.str();
  }
  
  void MpCCIexch::ClearDatasetTempFiles(std::string& stepNumsFileName,
                                        std::string& stepValuesFileName,
                                        std::string& resultScriptFileName,
                                        std::string& resultDatFileName,
                                        std::string& writeResultCmd,
                                        std::string& writeAuxCmd)
  {
    Settings& settings = Settings::Instance();
    std::ostringstream sstr;

    sstr.str("");
    sstr << settings.GetString("name") << "_data/StepNumbers.dat";
    stepNumsFileName = sstr.str();
    bfs::remove(bfs::path(stepNumsFileName));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/StepValues.dat";
    stepValuesFileName = sstr.str();
    bfs::remove(bfs::path(stepValuesFileName));

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Result.sh";
    resultScriptFileName = sstr.str();

    sstr.str("");
    sstr << settings.GetString("name") << "_data/Result.dat";
    resultDatFileName = sstr.str();

    sstr.str("");
    sstr << settings.GetString("baseExeDir") << "/hdf5_import_dataset.sh "
         << settings.GetString("name") << "_data "
         << settings.GetString("name") << ".h5 "
         << settings.GetString("type") << " "
         << settings.GetString("floatDataset") << " "
         << settings.GetInt("extfiles"); 
    
    writeResultCmd = sstr.str();

    sstr.str("");
    sstr << settings.GetString("baseExeDir") << "/hdf5_aux.sh "
         << settings.GetString("name") << "_data "
         << settings.GetString("name") << ".h5 "
         << settings.GetString("type") << " "
         << settings.GetString("deltemp");
    writeAuxCmd = sstr.str();
    
    sstr.str("");
    sstr << "rm -f " << settings.GetString("name") << "_data/TimeStep_*;\n";
    if(std::system(sstr.str().c_str()))
    {
      EXCEPTION("Having trouble to clean up HDF5 time step files.");
    }
  }
} // end of namespace
