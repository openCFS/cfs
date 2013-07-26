/* 
 * File:   MeshFile.cpp
 * Author: Mace
 * 
 * Created on 31. Mai 2013, 23:56
 */

#include "MeshFile.hh"
#include "DataFile.hh"

#include <vector>
#include <string>
#include <iostream>

#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmiotypes.h"

namespace CCM {
  
  MeshFile::MeshFile(std::string fileName):DataFile(fileName) {
    CollectChildEntities(root_, topologies_, kCCMIOTopology);
  }
  
  MeshFile::MeshFile(const MeshFile& orig):DataFile(orig.fileName_) {
  }
  
  MeshFile::~MeshFile() {
  }
  
  void MeshFile::ReadFaces(std::vector<FaceData>* faceDatas) {
    std::cout << "Reading Faces:" << std::endl;
    FaceData faceData;
    for (uint iTopology = 0; iTopology < topologies_.size(); iTopology++) {
      CCMIOID topologyID = topologies_[iTopology];
      //faceData.topolgyID = faceData.topolgyID;
      OpenLabel(topologyID, faceData.topologyLabel);
      std::cout << "  Topology: " << faceData.topologyLabel << " (id: " << faceData.topolgyID << ")" << std::endl;
      
      std::vector<CCMIOID> facesIDs;
      
      // reading internal faces
      CollectChildEntities(topologyID, facesIDs, kCCMIOInternalFaces);
      faceData.isBoundary = false;
      for (uint iFaces = 0; iFaces < facesIDs.size(); iFaces++) {
        CCMIOID facesID = facesIDs[iFaces];
        faceData.id = facesID.id;
        OpenLabel(facesID, faceData.label);
        CCMIOEntitySize(&err_, facesID, (CCMIOSize_t*) &faceData.size, NULL);
        std::cout << "    Faces: " << faceData.label << " (id: " << faceData.id << ", size: " << faceData.size << ")" << std::endl;
        
        CCMIOID mapID;
        CCMIOReadFaces(&err_, facesID, kCCMIOInternalFaces, &mapID, (CCMIOSize_t*) &faceData.faceVerticesSize, NULL, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceData.mapID = mapID.id;
        OpenLabel(mapID, faceData.mapLabel);
        std::cout << "      Map: " << faceData.mapLabel << " (id: " << faceData.mapID << ")" << std::endl;
        std::cout << "      Numer of Edges: " << (faceData.faceVerticesSize - faceData.size) << std::endl;
        faceData.faceVertices = new int[faceData.faceVerticesSize];
        CCMIOReadFaces(&err_, facesID, kCCMIOInternalFaces, NULL, (CCMIOSize_t*) 0, faceData.faceVertices, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceData.faceCells = new int[faceData.size * 2];
        CCMIOReadFaceCells(&err_, facesID, kCCMIOInternalFaces, faceData.faceCells, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceDatas->push_back(faceData);
        std::cout << "      Completed" << faceData.mapLabel << std::endl;
      }
      
      // reading boundary faces
      faceData.isBoundary = true;
      CollectChildEntities(topologyID, facesIDs, kCCMIOBoundaryFaces);
      for (uint iFaces = 0; iFaces < facesIDs.size(); iFaces++) {
        CCMIOID facesID = facesIDs[iFaces];
        faceData.id = facesID.id;
        OpenLabel(facesID, faceData.label);
        CCMIOEntitySize(&err_, facesID, (CCMIOSize_t*) &faceData.size, NULL);
        std::cout << "    " << faceData.label << " (id: " << faceData.id << ", size: " << faceData.size << ")" << std::endl;
        
        CCMIOID mapID;
        CCMIOReadFaces(&err_, facesID, kCCMIOBoundaryFaces, &mapID, (CCMIOSize_t*) &faceData.faceVerticesSize, NULL, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceData.mapID = mapID.id;
        OpenLabel(mapID, faceData.mapLabel);
        std::cout << "      Map: " << faceData.mapLabel << " (id: " << faceData.mapID << ")" << std::endl;
        std::cout << "      Numer of Edges: " << (faceData.faceVerticesSize - faceData.size) << std::endl;
        faceData.faceVertices = new int[faceData.faceVerticesSize];
        CCMIOReadFaces(&err_, facesID, kCCMIOBoundaryFaces, NULL, NULL, faceData.faceVertices, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceData.faceCells = new int[faceData.size];
        CCMIOReadFaceCells(&err_, facesID, kCCMIOBoundaryFaces, faceData.faceCells, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        faceDatas->push_back(faceData);
        std::cout << "      Completed" << faceData.mapLabel << std::endl;
      }
    }
    std::cout << "Reading Faces Completed" << std::endl << std::endl;
  }
  
  void MeshFile::ReadVertices(std::vector<VertexData>* vertexDatas) {
    std::cout << "Reading Vertices:" << std::endl;
    VertexData vertexData;
    std::vector<CCMIOID> verticesIDs;
    CollectChildEntities(root_, verticesIDs, kCCMIOVertices);
    for (uint iVertices=0; iVertices < verticesIDs.size(); iVertices++) {
      CCMIOID verticesID = verticesIDs[iVertices];
      vertexData.id = verticesID.id;
      OpenLabel(verticesID, vertexData.label);
      CCMIOEntitySize(&err_, verticesID, (CCMIOSize_t*) &vertexData.size, NULL);
      std::cout << "    " << vertexData.label << " (id: " << vertexData.id << ", size: " << vertexData.size << ")" << std::endl;
      
      CCMIOID mapID;
      vertexData.coordinates = new double[vertexData.size*3];
      CCMIOReadVerticesd(&err_, verticesID, (CCMIOSize_t*) &vertexData.dim, &vertexData.scale, &mapID, 
              vertexData.coordinates, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
      vertexData.mapID = mapID.id;
      OpenLabel(mapID, vertexData.mapLabel);
      std::cout << "        Map: " << vertexData.mapLabel << " (id: " << vertexData.mapID << ")" << std::endl;
      std::cout << "        Dim: " << vertexData.dim << std::endl;
      std::cout << "        Scale Factor " << vertexData.scale << std::endl;
      vertexDatas->push_back(vertexData);
    }
    std::cout << "Reading Vertices Completed" << std::endl << std::endl;
  }
  
  void MeshFile::ReadCells(std::vector<CellData>* cellDatas) {
    std::cout << "Reading Cells:" << std::endl;
    CellData cellData;
    for (uint iTopology = 0; iTopology < topologies_.size(); iTopology++) {
      CCMIOID topologyID = topologies_[iTopology];
      cellData.topolgyID = cellData.topolgyID;
      OpenLabel(topologyID, cellData.topologyLabel);
      std::cout << "  Topology: " << cellData.topologyLabel << " (id: " << cellData.topolgyID << ")" << std::endl;
      
      std::vector<CCMIOID> cellsIDs;
      
      // reading internal faces
      CollectChildEntities(topologyID, cellsIDs, kCCMIOCells);
      for (uint iCells = 0; iCells < cellsIDs.size(); iCells++) {
        CCMIOID cellsID = cellsIDs[iCells];
        cellData.id = cellsID.id;
        OpenLabel(cellsID, cellData.label);
        CCMIOEntitySize(&err_, cellsID, (CCMIOSize_t*) &cellData.size, NULL);
        std::cout << "    Cells: " << cellData.label << " (id: " << cellData.id << ", size: " << cellData.size << ")" << std::endl;
        
        CCMIOID mapID;
        CCMIOReadCells(&err_, cellsID, &mapID, NULL, (CCMIOIndex_t) kCCMIOStart, (CCMIOIndex_t) kCCMIOEnd);
        cellData.mapID = mapID.id;
        OpenLabel(mapID, cellData.mapLabel);
        std::cout << "      Map: " << cellData.mapLabel << " (id: " << cellData.mapID << ")" << std::endl;
        cellDatas->push_back(cellData);
      }
    }
    std::cout << "Reading Cells Completed" << std::endl << std::endl;
  }


}

