#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <math.h>
#include <unistd.h>

#include<boost/filesystem/operations.hpp>
namespace bfs=boost::filesystem;

#include <General/environment.hh>

#include "settings.hh"
#include "integlib.h"
#include "integlib/elemIntegr.hh"

#include "MpCCIexch.hh"


#if 0

#if (MPCCI_REL == 303)
#include <cci.h>
#else
#include <mpcci.h>
#endif

#endif

namespace CoupledField
{

    MpCCIexch::MpCCIexch(int argc, char *argv[], FileReader * ptFileReader)
    {
        Settings& settings = Settings::Instance();
        char buf[1024];

        ptFileReader_ = ptFileReader;
        
        if(!ptFileReader_)
          EXCEPTION("Invalid pointer to file reader!");

        ptFileReader_->Init();
#if 0
        quantityId1_ = 7;
        quantityDim1_ = 1;

        quantityId2_ = 8;
        quantityDim2_ = 3;

        //Define specific coupling description for cplreader side
        meshId_    = 1;
        GlobalDim_ = 3;

        CCI_Init_with_id_string(&argc, &argv, "MESHIO");
#endif
    }

    MpCCIexch::~MpCCIexch()
    {
#if 0
        CCI_Finalize();
#endif
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
      std::vector<int>::iterator it, it2, it3, end;
      int localNodeNum = 1;
      std::set<int> regionNodeSet;
      int maxNumElemNodes = 0;

      std::ofstream outFile;
      std::ostringstream sstr;
      std::string fileName;

      sstr.str("");
      sstr << settings.GetString("name") << "_data/Coordinates.cfg";
      std::string coordCfgFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));

      sstr.str("");
      sstr << settings.GetString("name") << "_data/Coordinates.dat";
      std::string coordDatFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));

      sstr.str("");
      sstr << settings.GetString("name") << "_data/Connectivity.cfg";
      std::string connCfgFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));
      
      sstr.str("");
      sstr << settings.GetString("name") << "_data/Connectivity.dat";
      std::string connDatFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));

      sstr.str("");
      sstr << settings.GetString("name") << "_data/Types.cfg";
      std::string typesCfgFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));
      
      sstr.str("");
      sstr << settings.GetString("name") << "_data/Types.dat";
      std::string typesDatFileName = sstr.str();
      bfs::remove(bfs::path(sstr.str()));

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
        it = Topology_[actPart].begin ();
        end = Topology_[actPart].end ();
        std::set<int>::iterator regSetIt, regSetEnd;
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
        //            std::copy (NodalCoords_.begin (), NodalCoords_.end (), double_os_iter (outFile, "\n"));

        //            outFile.write((char*)&Topology_[0], sizeof(int)*nElems*numNodesPerElem_[0]) << std::endl;
        outFile.close();
        
        it = elemTypes_[actPart].begin();
        end = elemTypes_[actPart].end();
        it2 = Topology_[actPart].begin();
        std::vector<int> newTopo;
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
        
        outFile.open(connDatFileName.c_str(),
            std::ofstream::out | std::ofstream::app);

        /*
        it = Topology_[actPart].begin ();
        end = Topology_[actPart].end ();
        for( ; it != end; it++ )
        {
          *it = nodesInPartition_[actPart][*it];
          outFile << nodesInPartition_[actPart][*it] << std::endl;
        }
*/
        typedef std::ostream_iterator<int> int_os_iter;

        std::copy (newTopo.begin (), newTopo.end (), int_os_iter (outFile, "\n"));
        //            outFile.write((char*)&Topology_[0], sizeof(int)*nElems*numNodesPerElem_[0]) << std::endl;
        outFile.close();


        outFile.open(typesDatFileName.c_str(),
            std::ofstream::out | std::ofstream::app);
        
        std::copy (elemTypes_[actPart].begin (), elemTypes_[actPart].end (), int_os_iter (outFile, "\n"));
        //            outFile.write((char*)&elemTypes_[0], sizeof(int)*nElems) << std::endl;
        outFile.close();


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

//        std::copy (regionNodeSet.begin (), regionNodeSet.end (), int_os_iter (outFile, "\n"));
        outFile.close();

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

        sstr.str("");
        sstr << settings.GetString("name") << "_data/RegionDim_partition" << (actPart+1) << ".dat";
        fileName = sstr.str();
        outFile.open(fileName.c_str(),
            std::ofstream::out | std::ofstream::trunc);

        outFile << dim << std::endl;
        outFile.close();

        elemOffset += nElems;
        nodeOffset += nNodes;
        
#if 0
        MESHIO_DEBUG( "Define partition" );
        // MpCCI part
        CCI_Def_partition(meshId_, actPart+1);
        std::cout << "Define nodes" << std::endl;
        //define the nodes
        CCI_Def_nodes(meshId_, actPart+1, ptFileReader_->GetDim(),
            ptFileReader_->GetNumNodes(), 0, NULL,
            CCI_DOUBLE, ptFileReader_->GetNodePtr());


        MESHIO_DEBUG( "Define elements" );
        //define the elements
        CCI_Def_elems(meshId_, actPart+1, nElems, 0, NULL, 
            nElems, &elemTypes[0], &numNodesPerElem[0],
            &Topology[0]);
#endif

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
      outFile.close();

      outFile.open(connCfgFileName.c_str(),
          std::ofstream::out | std::ofstream::trunc);
      outFile << "PATH /Mesh/Elements/Connectivity" << std::endl;
      outFile << "INPUT-CLASS TEXTIN" << std::endl;
      outFile << "INPUT-SIZE 32" << std::endl;
      outFile << "RANK 2" << std::endl;
      outFile << "DIMENSION-SIZES " << elemOffset << " "<< maxNumElemNodes << std::endl;
      outFile << "OUTPUT-CLASS UIN" << std::endl;
      outFile << "OUTPUT-SIZE 32" << std::endl;
      outFile << "OUTPUT-ARCHITECTURE STD" << std::endl;
      outFile << "OUTPUT-BYTE-ORDER LE" << std::endl;
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
      outFile.close();
      
      
      sstr.str("");
      sstr << settings.GetString("name") << "_data/"
      << settings.GetString("name") << ".h5";
      bfs::remove(bfs::path(sstr.str()));

      sstr.str("");
      sstr << settings.GetString("baseExeDir") << "/hdf5_import_mesh.sh "
      << settings.GetString("name") << "_data "
      << settings.GetString("name") << ".h5";

      if(std::system(sstr.str().c_str()))
      {
        exit(1);
      }

      //Close the definition phase; contact detection.
      //i take part in the coupling
#if 0        
      MESHIO_DEBUG( "Closing MpCCI Setup" );
      MpCCIprocess_ = 1;
      CCI_Close_setup(MpCCIprocess_);
#endif

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
#if 0
        int globalConvergence = CCI_CONTINUE;
        int myConvergence     = CCI_CONTINUE;
#endif
        int counter = 0;
        int numFiles = ptFileReader_->GetNumFiles();
        int numPartitions = ptFileReader_->GetNumPartitions();
        std::ostringstream sstr;
        
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

        sstr.str("");
        sstr << settings.GetString("name") << "_data/StepNumbers.dat";
        bfs::remove(bfs::path(sstr.str()));

        sstr.str("");
        sstr << settings.GetString("name") << "_data/StepValues.dat";
        bfs::remove(bfs::path(sstr.str()));
        

        while ( counter < numFiles ) 
        {
            int nodeOffset = 0;
            int actPart=0;
            double stepVal = counter*settings.GetDouble("timeStep");
            std::cout << "----------------------------------------"
                      << "----------------------------------------"
                      << std::endl;
            std::cout << "                        "
                      << "Step Number: " << (counter+1) << "; ";
            std::cout << " Time Step: " << stepVal << std::endl;
            std::cout << "----------------------------------------"
                      << "----------------------------------------"
                      << std::endl;
            

            for (; actPart<numPartitions; actPart++)
            {
              int numNodesPart = ptFileReader_->GetNumNodes(actPart);
//              int numResults = ptFileReader_->GetNumResults();

              // acoustic source
              std::vector<double> acouSrc;
              acouSrc.resize(numNodesPart); 

              //*3 due to velx, vely, velz (for 2D=0)
              std::vector<double> velocity;
              velocity.resize(3*numNodesPart); 

              std::vector<double> flowdata_;
              ptFileReader_->ReadNodalValues(flowdata_,
                                             actPart,
                                             counter);


              // by convention the flowdata array has to have
              // the following structure:
              // flowdata [0]      [1]  [2]  [3]  [4]   [5]   [6]   [7]      [8] ... 
              //          acousrc1 vx1  vy1  vz1  usr11 usr21 usr31 acousrc2 vx2 ...

              // Putting separately acouSrc (scalar) and Vel (vector)

              int k = 0;
              for (int inode=0; inode<numNodesPart; inode++)
              {
                acouSrc[inode] = flowdata_[inode*7];

                // acouSrc[inode] = 0;

                velocity[k] = flowdata_[inode*7+1];
                velocity[k+1] = flowdata_[inode*7+2];
                velocity[k+2] = flowdata_[inode*7+3];


                k = k+3;
              }

              if(calcSrc)
                CalculateAcouSrcs(actPart, velocity, acouSrc);

              std::vector<double> tmpDat;
              std::map<int, int>::iterator it, end;
              it = nodesInPartition_[actPart].begin();
              end = nodesInPartition_[actPart].end();
              tmpDat.resize(nodesInPartition_[actPart].size());
              int newNodeIdx = 0;
              for( ; it != end; it++ ) {
                int oldNodeIdx = it->first-1;
                tmpDat[newNodeIdx] = acouSrc[oldNodeIdx];
                newNodeIdx++;                           
              }
              
              std::ofstream outFile;
              typedef std::ostream_iterator<double> double_os_iter;
              std::string fileName;
                
              sstr.str("");
              sstr << settings.GetString("name") << "_data/Result.sh";
              fileName = sstr.str();
              outFile.open(fileName.c_str(),
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
              outFile << "RESULT_DIMENSION_SIZES=" << tmpDat.size() << std::endl;
              outFile.close();

              sstr.str("");
              sstr << settings.GetString("name") << "_data/Result.dat";
              fileName = sstr.str();
              outFile.open(fileName.c_str(),
                           std::ofstream::out | std::ofstream::trunc);

              std::copy (tmpDat.begin (), tmpDat.end (), double_os_iter (outFile, "\n"));
              outFile.close();

              sstr.str("");
              sstr << settings.GetString("baseExeDir") << "/hdf5_import_dataset.sh "
                   << settings.GetString("name") << "_data "
                   << settings.GetString("name") << ".h5";
              if(std::system(sstr.str().c_str()))
                exit(1);

              sstr.str("");
              sstr << settings.GetString("name") << "_data/StepNumbers.dat";
              fileName = sstr.str();
              outFile.open(fileName.c_str(),
                           std::ofstream::out | std::ofstream::app);
              outFile << (counter+1) << std::endl;
              outFile.close();

              sstr.str("");
              sstr << settings.GetString("name") << "_data/StepValues.dat";
              fileName = sstr.str();
              outFile.open(fileName.c_str(),
                           std::ofstream::out | std::ofstream::app);
              outFile << stepVal << std::endl;
              outFile.close();

           #if 0
              // acouSrc scalar
              CCI_Put_nodes( meshId_, actPart+1, quantityId1_,
                             quantityDim1_, numNodesPart, 0, NULL,
                             CCI_DOUBLE, &acouSrc[0]);
                    
              // Vx, Vy, Vz
              CCI_Put_nodes( meshId_, actPart+1, quantityId2_,
                             quantityDim2_, numNodesPart, 0, NULL,
                             CCI_DOUBLE, &velocity[0]);
           #endif

              nodeOffset += numNodesPart;
                

            }//end of for

         #if 0
            //Send values via MpCCI
            int nQuantityIds=0;
            int quantityIds[] = {0, 0};
            int nLocalMeshIds=1;
            int localMeshIds[] = {1};
            int remoteCodeId = 2;
                
            int comm = CCI_COMM_RCODE[remoteCodeId];

            quantityIds[nQuantityIds] = quantityId1_;
            nQuantityIds++;
            quantityIds[nQuantityIds] = quantityId2_;
            nQuantityIds++;
                
            CCI_Send(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, comm);
                
            CCI_Check_convergence(myConvergence,&globalConvergence,CCI_ANY_CODE);
                
            std::cout<<"CCI_CONTINUE_VALUE"<<CCI_CONTINUE<<std::endl;
         #endif
            counter++;
        }//end of while
     
        sstr.str("");
        sstr << settings.GetString("baseExeDir") << "/hdf5_aux.sh "
             << settings.GetString("name") << "_data "
             << settings.GetString("name") << ".h5 "
             << settings.GetString("type") << " "
             << settings.GetString("deltemp");
        if(std::system(sstr.str().c_str()))
          exit(1);

        std::cout << "========================================"
                  << "========================================"
                  << std::endl;
        std::cout << "                        "
                  << "Leaving Time Data Conversion Phase!" << std::endl;
        std::cout << "========================================"
                  << "========================================"
                  << std::endl;
    }

#if 0
    int MpCCIexch::ElemTypes2MpCCI(elemType et)
    {
        switch(et)
        {
        case ET_LINE2:
            return CCI_ELEM_LINE;
            break;

        case ET_TRIA3:
            return CCI_ELEM_TRIANGLE;
            break;

        case ET_QUAD4:
            return CCI_ELEM_QUAD;
            break;

        case ET_QUAD8:
            return CCI_ELEM_QUAD8;
            break;

        case ET_TET4:
            return CCI_ELEM_TETRAHEDRON;
            break;
                
        case ET_WEDGE6:
            return CCI_ELEM_PRISM;
            break;

        case ET_HEXA8:
            return CCI_ELEM_HEXAHEDRON;
            break;

        case ET_PYRA5:
            return CCI_ELEM_PYRAMID;
            break;
        }
    }
#endif
    
    void MpCCIexch::CalculateAcouSrcs(const int partitionIdx,
                                      const std::vector<double>& velocity,
                                      // const std::vector<double>& pressure,
                                      std::vector<double>& acouSrc)
    {
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
                    coordMatQuad[0][n] = NodalCoords_[partitionIdx][(Topology_[partitionIdx][k+n]-1)*3+0];
                    coordMatQuad[1][n] = NodalCoords_[partitionIdx][(Topology_[partitionIdx][k+n]-1)*3+1];
                    
                    nodalVelQuad[0][n] = velocity[(Topology_[partitionIdx][k+n]-1)*3+0];
                    nodalVelQuad[1][n] = velocity[(Topology_[partitionIdx][k+n]-1)*3+1];
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
                    acouSrc[Topology_[partitionIdx][k+n]-1] += elemVecQuad[n];
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
                int idx = (Topology_[partitionIdx][k+n]-1)*3;
                
                coordMatHex[0][n] = NodalCoords_[partitionIdx][idx+0];
                coordMatHex[1][n] = NodalCoords_[partitionIdx][idx+1];
                coordMatHex[2][n] = NodalCoords_[partitionIdx][idx+2];
                
                nodalVelHex[0][n] = velocity[idx+0];
                nodalVelHex[1][n] = velocity[idx+1];
                nodalVelHex[2][n] = velocity[idx+2];
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
                acouSrc[Topology_[partitionIdx][k+n]-1] += elemVecHex[n];
              }
              break;
              
            case ET_PYRA5:
                break;
            }

            k += numNodesPerElem_[partitionIdx][i];
        }
    }

} // end of namespace
