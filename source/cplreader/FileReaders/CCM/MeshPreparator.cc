/* 
 * File:   MeshPreparator.cpp
 * Author: Mace
 * 
 * Created on 2. Juni 2013, 15:15
 */

#include "MeshPreparator.hh"

#include <string>
#include <vector>
#include <iostream>

#include "MeshFile.hh"
#include "Utilities.hh"
#include "Meshes.hh"
#include "divUtil.h"
#include "PreAllocatedDataAcceptor.hh"
#include "ValueMapCatcher.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace CCM {
  
  MeshPreparator::MeshPreparator() {
  }
  
  MeshPreparator::MeshPreparator(const MeshPreparator& orig) {
  }
  
  MeshPreparator::~MeshPreparator() {
  }
  
  void MeshPreparator::PrepareMesh(MeshFile* meshFile, DualizableMesh& mesh, std::vector<ConsecutiveMap>& valueConMaps) {
    meshFile_ = meshFile;
    
    ReadMesh();
    
    
    RemapRawData();
    ClearMaps();
    
    CreateVertices();
    ClearVertices();
    
    DetermineFaceCounts();
    CreateFaceToVertex();
    CreateFaceToCell();
    ClearFaces();
    
    CreateCellToFace();
    
    MarkGhosts();
    CreateGhostFaces();
    
    ReadCellPositions();
    
    for (uint i=0; i < cellConMaps_.size(); i++) {
      valueConMaps.push_back(cellConMaps_[i]);
    }
    
    mesh = m_;
    meshFile_ = NULL;
  }
  
  void MeshPreparator::ReadMesh() {
    std::cout << "READING MESH FILE" << std::endl << std::endl;
    meshFile_->ReadMaps(&maps_);
    meshFile_->ReadFaces(&faces_);
    meshFile_->ReadVertices(&vertices_);
  }
  
  
  
  void MeshPreparator::RemapRawData() {
    std::cout << "Remapping raw data consecutively" << std::endl;
    
    std::cout << "  Creating consecutive cell maps" << std::endl;
    std::vector<Map> cellMaps;
    std::vector<Map> boundaryFaceMaps;
    ValueMapCatcher vmc(maps_, cellMaps, boundaryFaceMaps);
    meshFile_->SetVerbose(false);
    meshFile_->ReadFields(NULL, &vmc);
    
    uint numCells = 0;
    uint numRealCells = 0;
    uint numGhostCells = 0;
    uint maxCellValue = 0;
    for (uint i=0; i < cellMaps.size(); i++) {
      Map map = cellMaps[i];
      maxCellValue = map.maxValue > maxCellValue ? map.maxValue : maxCellValue;
      ConsecutiveMap conMap;
      conMap.id = map.id;
      conMap.label = map.label;
      conMap.size = map.size;
      conMap.startValue = numCells + 1;
      numCells += cellMaps[i].size;
      cellConMaps_.push_back(conMap);
      std::cout << "    " << conMap.startValue << " - " << (conMap.startValue + conMap.size - 1) << ": " << cellMaps[i].label << std::endl;
    }
    numRealCells = numCells;
    cellConMapsForFaces_.resize(faces_.size());
    for (uint i=0; i < boundaryFaceMaps.size(); i++) {
      Map map = boundaryFaceMaps[i];
      ConsecutiveMap conMap;
      conMap.id = map.id;
      conMap.label = map.label;
      conMap.size = map.size;
      conMap.startValue = numCells + 1;
      numCells += map.size;
      numGhostCells += map.size;
      cellConMaps_.push_back(conMap);
      for (uint j=0; j < faces_.size(); j++) {
        if (faces_[j].mapID == map.id) {
          cellConMapsForFaces_[j] = conMap;
        }
      }
      std::cout << boundaryFaceMaps.size() << "    " << conMap.startValue << " - " << (conMap.startValue + conMap.size - 1) << ": " << boundaryFaceMaps[i].label << std::endl;
    }
    m_.cellCount = numCells;
    m_.cellGhostCount = numGhostCells;
    std::cout << "    Real cells:  " << numRealCells << std::endl;
    std::cout << "    Ghost cells: " << numGhostCells << std::endl;
    std::cout << "    All cells:   " << numCells << std::endl;
    
    std::cout << "  Creating cell remaps" << std::endl;
    int* rawCellCellMap = new int[maxCellValue + 1];
    int cell = 1;
    for (uint i=0; i < cellMaps.size(); i++) {
      Map map = cellMaps[i];
      for (uint j=0; j < map.size; j++) {
        rawCellCellMap[map.value[j]] = cell;
        cell++;
      }
    }
    
    std::cout << "  Applying cell remaps to raw Face -> Cell Data" << std::endl;
    for (uint i=0; i < faces_.size(); i++) {
      uint size = faces_[i].isBoundary ? faces_[i].size : faces_[i].size * 2;
      int* faceCells = faces_[i].faceCells;
      for (uint j=0; j< size; j++) {
        faceCells[j] = rawCellCellMap[faceCells[j]];
      }
    }
    delete[] rawCellCellMap;
    
    std::cout << "  Creating consecutive vertex maps" << std::endl;
    uint numVertices = 0;
    uint maxVertexValue = 0;
    for (uint i=0; i < vertices_.size(); i++) {
      Map map = GetMap(&vertices_[i]);
      maxVertexValue = map.maxValue > maxVertexValue ? map.maxValue : maxVertexValue;
      ConsecutiveMap conMap;
      conMap.id = map.id;
      conMap.label = map.label;
      conMap.size = map.size;
      conMap.startValue = numVertices + 1;
      numVertices += vertices_[i].size;
      vertexConMaps_.push_back(conMap);
      std::cout << "    " << conMap.startValue << " - " << (conMap.startValue + conMap.size) << ": " << vertices_[i].label << std::endl;
    }
    std::cout << "  Creating vertex remaps" << std::endl;
    int* rawVertexVertexMap = new int[maxVertexValue + 1];
    int vertex = 1;
    for (uint i=0; i < vertices_.size(); i++) {
      Map map = GetMap(&vertices_[i]);
      for (uint j=0; j < map.size; j++) {
        rawVertexVertexMap[map.value[j]] = vertex;
        vertex++;
      }
    }
    m_.vertexCount = numVertices;
    std::cout << "    Vertices:    " << numVertices << std::endl;
    std::cout << "  Applying cell remaps to raw Face -> Vertex Data" << std::endl;
    for (uint i=0; i < faces_.size(); i++) {
      uint size = faces_[i].faceVerticesSize;
      int* faceVertices = faces_[i].faceVertices;
      uint index = 0;
      for (uint j=0; j < size; j++) {
        if (j == index) {
          index += faceVertices[index] + 1;
        } else {
          faceVertices[j] = rawVertexVertexMap[faceVertices[j]];
        }
      }
    }
    delete[] rawVertexVertexMap;
  }
  
  Map MeshPreparator::GetMap(MappedData* mappedData) {
    Map map;
    for (uint i=0; i < maps_.size(); i++) {
      if (maps_[i].id == mappedData->mapID) {
        return maps_[i];
      }
    }
    return map;
  }

  void MeshPreparator::ClearMaps() {
    std::cout << "Clearing raw map data" << std::endl;
    while (maps_.size() > 0) {
      Map map = maps_.back();
      maps_.pop_back();
      int* array = map.value;
      map.value = NULL;
      delete[] array;
    }
  }
  
  
  
  void MeshPreparator::CreateVertices() {
    std::cout << "Creating Vertex Coordinates" << std::endl;
    m_.vertexPosition = CreateZeroArray<double>(3*(m_.vertexCount + 1));
    uint actPos = 3;
    for (uint i=0; i < vertices_.size(); i++) {
      VertexData vd = vertices_[i];
      double* data = &m_.vertexPosition[actPos];
      uint size = vd.size * 3;
      for (uint j=0; j < size; j++) {
        data[j] = vd.coordinates[j];
      }
      actPos += size;
    }
  }
  
  void MeshPreparator::ClearVertices() {
    std::cout << "Clearing raw vertex data" << std::endl;
    for (uint i=0; i < vertices_.size(); i++) {
      VertexData vd = vertices_[i];
      double* array = vd.coordinates;
      vd.coordinates = NULL;
      delete[] array;
    }
    vertices_.clear();
  }
  
  
  
  void MeshPreparator::DetermineFaceCounts() {
    uint numRealFaces = 0;
    uint numBoundaryEdges = 0;
    uint numEdges = 0;
    std::cout << "Determining number of Faces" << std::endl;
    for (uint i=0; i < faces_.size(); i++) {
      FaceData fd = faces_[i];
      numRealFaces += fd.size;
      numEdges += (fd.faceVerticesSize - fd.size);
      if (fd.isBoundary) {
        numBoundaryEdges += (fd.faceVerticesSize - fd.size);
      }
    }
    std::cout << "    Real Faces:  " << numRealFaces << std::endl;
    uint numGhostFaces = numBoundaryEdges / 2;
    std::cout << "    Ghost Faces: " << numGhostFaces << std::endl;
    uint numFaces = numRealFaces + numGhostFaces;
    std::cout << "    All Faces:   " << numFaces << std::endl;
    std::cout << "    Real Face -> Vertices:  " << numEdges << std::endl;
    std::cout << "    Ghost Face -> Vertices: " << (numGhostFaces * 3) << std::endl;
    numEdges += numGhostFaces * 3;
    std::cout << "    All Face -> Vertices:   " << numEdges << std::endl;
    m_.faceCount = numFaces;
    m_.faceGhostCount = numGhostFaces;
    m_.faceVertexLength = numEdges;
    m_.cellFaceLength = numFaces * 2;
  }

  void MeshPreparator::CreateFaceToVertex() {
    std::cout << "Creating Face -> Vertex Data (Except Ghost Faces)" << std::endl;
    m_.faceVertexCount = CreateZeroArray<uint>(m_.faceCount + 1);
    m_.faceVertexIndex = CreateZeroArray<uint>(m_.faceCount + 1);
    m_.faceVertex = CreateZeroArray<uint>(m_.faceVertexLength + 1);
    uint faceNumber = 1;
    uint faceVertexIndex = 0;
    for (uint i=0; i < faces_.size(); i++) {
      FaceData fd = faces_[i];
      uint size = fd.size;
      uint iFaceVertexIndex = 0;
      for (uint j=0; j< size; j++) {
        uint numVerts = fd.faceVertices[iFaceVertexIndex];
        iFaceVertexIndex++;
        m_.faceVertexCount[faceNumber] = numVerts;
        m_.faceVertexIndex[faceNumber] = faceVertexIndex;
        for (uint k=0; k < numVerts; k++) {
          m_.faceVertex[faceVertexIndex] = fd.faceVertices[iFaceVertexIndex];
          iFaceVertexIndex++;
          faceVertexIndex++;
        }
        faceNumber++;
      }
    }
    // caring about ghost faces;
    while (faceNumber <= m_.faceCount) {
      m_.faceVertexCount[faceNumber] = 3;
      m_.faceVertexIndex[faceNumber] = faceVertexIndex;
      faceVertexIndex += 3;
      faceNumber++;
    }
  }
  
  void MeshPreparator::CreateFaceToCell() {
    std::cout << "Creating Face -> Cell Data (Except Ghost Faces), Counting Faces per Cell" << std::endl;
    m_.faceCell = CreateZeroArray<uint>(2*(m_.faceCount + 1));
    m_.cellFaceCount = CreateZeroArray<uint>(m_.cellCount + 1);
    
    uint faceCellIndex = 2;
    for (uint i=0; i < faces_.size(); i++) {
      FaceData fd = faces_[i];
      uint size = fd.size;
      if (!fd.isBoundary) {
        size *= 2;
        for (uint j=0; j < size; j++) {
          m_.cellFaceCount[fd.faceCells[j]]++;
          m_.faceCell[faceCellIndex] = fd.faceCells[j];
          faceCellIndex++;
        }
      } else {
        uint ghostCell = cellConMapsForFaces_[i].startValue;
        for (uint j=0; j < size; j++) {
          m_.cellFaceCount[fd.faceCells[j]]++;
          m_.faceCell[faceCellIndex] = fd.faceCells[j];
          m_.cellFaceCount[ghostCell] += 1 + m_.faceVertexCount[faceCellIndex / 2];
          faceCellIndex++;
          m_.faceCell[faceCellIndex] = ghostCell;
          faceCellIndex++;
          ghostCell++;
        }
      }
    }
  }

  void MeshPreparator::ClearFaces() {
    std::cout << "Clearing raw face data" << std::endl;
    for (uint i=0; i < faces_.size(); i++) {
      FaceData fd = faces_[i];
      int* array = fd.faceVertices;
      fd.faceVertices = NULL;
      delete[] array;
      array = fd.faceCells;
      fd.faceCells = NULL;
      delete[] array;
    }
    faces_.clear();
  }
  
  
  
  void MeshPreparator::CreateCellToFace() {
    std::cout << "Creating Cell -> Face Data" << std::endl;
    m_.cellFaceIndex = CreateZeroArray<uint>(m_.cellCount + 1);
    m_.cellFace = CreateZeroArray<uint>(m_.cellFaceLength);
    
    uint cellFaceIndex = 0;
    uint maxCell = m_.cellCount;
    for (uint i=1; i <= maxCell; i++) {
      m_.cellFaceIndex[i] = cellFaceIndex;
      cellFaceIndex += m_.cellFaceCount[i];
    }
    
    uint maxFace = m_.faceCount;
    for (uint i=1; i <= maxFace; i++) {
      AddFaceToCell(m_.faceCell[i*2], i);
      AddFaceToCell(m_.faceCell[i*2 + 1], i);
    }
  }
  
  inline void MeshPreparator::AddFaceToCell(uint& cell, uint& face) {
    if (cell == 0) {
      return;
    }
    uint index = m_.cellFaceIndex[cell];
    while (m_.cellFace[index] != 0) {
      index++;
    }
    m_.cellFace[index] = face;
  }
  
  
  void MeshPreparator::MarkGhosts() {
    std::cout << "Marking Ghost Cells and Faces" << std::endl;
    m_.cellGhost.clear();
    m_.cellGhost.resize(m_.cellCount - m_.cellGhostCount + 1, false);
    m_.cellGhost.resize(m_.cellCount + 1, true);
    m_.faceGhost.clear();
    m_.faceGhost.resize(m_.faceCount - m_.faceGhostCount + 1, false);
    m_.faceGhost.resize(m_.faceCount + 1, true);
  }
  
  void MeshPreparator::CreateGhostFaces() {
    std::cout << "Creating Ghost Faces" << std::endl;
    std::vector<bool> usedEdges;
    usedEdges.resize(m_.faceVertexLength, false);
    uint counter = 0;
    uint ghostFaceNumber = m_.faceCount - m_.faceGhostCount + 1;
    for (uint iCell = 1; iCell <= m_.cellCount; iCell++) {
      if (m_.cellGhost[iCell]) {
        uint startFace = m_.cellFace[m_.cellFaceIndex[iCell]];
        uint startFaceVertexIndex = m_.faceVertexIndex[startFace];
        uint startFaceVertexCount = m_.faceVertexCount[startFace];
        for (uint iEdge=0; iEdge < startFaceVertexCount; iEdge++) {
          if (!usedEdges[startFaceVertexIndex + iEdge]) {
            usedEdges[startFaceVertexIndex + iEdge] = true;
            uint vertexA = m_.faceVertex[startFaceVertexIndex + iEdge];
            uint vertexB = iEdge + 1 == startFaceVertexCount 
                    ? m_.faceVertex[startFaceVertexIndex] : m_.faceVertex[startFaceVertexIndex + iEdge + 1];
            uint nextCell = m_.faceCell[startFace * 2];
            uint nextFace = startFace;
            uint nextnextCell;
            uint nextnextFace;
            uint nextnextVertexIndex;
            bool nextnextIsPrimary;
            while (!m_.cellGhost[nextCell]) {
              TraverseEdge(vertexA, vertexB, 
                  nextCell, nextFace, &m_, nextnextCell, nextnextFace, 
                      nextnextVertexIndex, nextnextIsPrimary);
              nextCell = nextnextCell;
              nextFace = nextnextFace;
              usedEdges[nextnextVertexIndex] = true;
            }
            m_.faceCell[ghostFaceNumber * 2] = iCell;
            m_.faceCell[ghostFaceNumber * 2+1] = nextCell;
            m_.faceVertex[m_.faceVertexIndex[ghostFaceNumber]] = vertexA;
            m_.faceVertex[m_.faceVertexIndex[ghostFaceNumber] + 1] = vertexB;
            m_.faceVertex[m_.faceVertexIndex[ghostFaceNumber] + 2] = 0;
            AddFaceToCell(iCell, ghostFaceNumber);
            AddFaceToCell(nextCell, ghostFaceNumber);
            counter++;
            ghostFaceNumber++;
          }
        }
      }
    }
  }

  void MeshPreparator::ReadCellPositions() {
    std::cout << "Reading Cell Positions" << std::endl;
    PreAllocatedDataAcceptor prep(cellConMaps_);
    
    std::cout << "  Reading X Positions" << std::endl;
    double* posX = new double[m_.cellCount + 1];
    prep.SetFieldName("Centroid_0");
    prep.data = &posX[1];
    meshFile_->ReadFields(NULL, &prep);
    
    std::cout << "  Reading Y Positions" << std::endl;
    double* posY = new double[m_.cellCount + 1];
    prep.SetFieldName("Centroid_1");
    prep.data = &posY[1];
    meshFile_->ReadFields(NULL, &prep);
    
    std::cout << "  Reading Z Positions" << std::endl;
    double* posZ = new double[m_.cellCount + 1];
    prep.SetFieldName("Centroid_2");
    prep.data = &posZ[1];
    meshFile_->ReadFields(NULL, &prep);
    
    std::cout << "  Reordering Data" << std::endl;
    m_.cellPosition = new double[3*(m_.cellCount + 1)];
#pragma omp parallel for
    for (uint i=1; i <= m_.cellCount; i++) {
      m_.cellPosition[i*3] = posX[i];
      m_.cellPosition[i*3+1] = posY[i];
      m_.cellPosition[i*3+2] = posZ[i];
    }
    
    delete[] posX;
    delete[] posY;
    delete[] posZ;
  }
  
}

