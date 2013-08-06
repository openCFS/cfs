/* 
 * File:   DualMesher.cpp
 * Author: tz
 * 
 * Created on 3. Juni 2013, 18:29
 */

#include "DualMesher.hh"


#include <string>
#include <vector>
#include <iostream>

#include "Meshes.hh"
#include "divUtil.h"

#ifdef _OPENMP
#include "omp.h"
#endif

namespace CCM {
  
  DualMesher::DualMesher() {
  }
  
  DualMesher::DualMesher(const DualMesher& orig) {
  }
  
  DualMesher::~DualMesher() {
  }
  
  void DualMesher::DualizeMesh(DualizableMesh& mesh, DualizedMesh& dual, bool copyCoordinates) {
    std::cout << "Dualizing Mesh" << std::endl;
    mesh_ = mesh;
    dual_ = dual;
    
    dual_.cellCount = mesh_.vertexCount;
    dual_.vertexCount = mesh_.cellCount;
    dual_.cellGhostCount = 0;
    dual_.faceGhostCount = 0;
    dual_.cellGhost.clear();
    dual_.cellGhost.resize(dual_.cellCount + 1, false);
    
    CreateFaces();
    CreateCellFaces();
    CheckCollapsedFaces();
    if (copyCoordinates) {
      CopyCoordinates();
    } else {
      dual_.cellPosition = mesh_.vertexPosition;
      dual_.vertexPosition = mesh_.cellPosition;
    }
    
    mesh = mesh_;
    dual = dual_;
    std::cout << "Dualizing Mesh Finished" << std::endl;
  }
  
  void DualMesher::CreateFaces() {
    std::cout << "Creating Dual Faces" << std::endl;
    uint* faceVertexCount = new uint[(mesh_.faceVertexLength + 6)/3];
    uint* faceVertexIndex = new uint[(mesh_.faceVertexLength + 6)/3];
    dual_.faceVertexLength = mesh_.faceVertexLength;
    dual_.faceVertex = new uint[dual_.faceVertexLength];
    uint* faceCell = new uint[2*(mesh_.faceVertexLength + 6)/3];
    std::vector<bool> usedEdges(mesh_.faceVertexLength, false);
    
    uint dualFaceCount = 0;
    uint posFaceVertexIndex = 0;
    
    for (uint iCell=1; iCell <= mesh_.cellCount; iCell++) {
      uint iCellFaceIndex = mesh_.cellFaceIndex[iCell];
      uint maxiCellFaceIndex = iCellFaceIndex + mesh_.cellFaceCount[iCell];
      for (;iCellFaceIndex < maxiCellFaceIndex; iCellFaceIndex++) {
        uint iFace = mesh_.cellFace[iCellFaceIndex];
        if (mesh_.faceCell[iFace*2] == iCell) {
          uint miniFaceVertexIndex = mesh_.faceVertexIndex[iFace];
          uint maxiFaceVertexIndex = miniFaceVertexIndex + mesh_.faceVertexCount[iFace];
          for (uint iFaceVertexIndex = miniFaceVertexIndex; iFaceVertexIndex < maxiFaceVertexIndex; iFaceVertexIndex++) {
            if (!usedEdges[iFaceVertexIndex]) {
              dualFaceCount++;
              uint vertexA = mesh_.faceVertex[iFaceVertexIndex];
              uint vertexB = iFaceVertexIndex + 1 == maxiFaceVertexIndex 
                      ? mesh_.faceVertex[miniFaceVertexIndex] : mesh_.faceVertex[iFaceVertexIndex+1];
              faceVertexCount[dualFaceCount] = 1;
              faceVertexIndex[dualFaceCount] = posFaceVertexIndex;
              dual_.faceVertex[posFaceVertexIndex] = iCell;
              posFaceVertexIndex++;
              faceCell[dualFaceCount*2] = vertexA;
              faceCell[dualFaceCount*2+1] = vertexB;
              usedEdges[iFaceVertexIndex] = true;
              uint cell = mesh_.faceCell[iFace*2+1];
              uint face = iFace;
              uint nextCell;
              uint nextFace;
              uint nextVertexIndex;
              bool isPrimary;
              bool found = true;
              while (cell != iCell && found) {
                found = TraverseEdge(vertexA, vertexB, cell, face, &mesh_,
                    nextCell, nextFace, nextVertexIndex, isPrimary);
                usedEdges[nextVertexIndex] = true;
                dual_.faceVertex[posFaceVertexIndex] = cell;
                posFaceVertexIndex++;
                faceVertexCount[dualFaceCount]++;
                cell = nextCell;
                face = nextFace;
              }
            }
          }
        }
      }
    }
    
    dual_.faceCount = dualFaceCount;
    dual_.faceVertexCount = new uint[dualFaceCount + 1];
#pragma omp parallel for
    for (uint i=1; i <= dualFaceCount + 1; i++) {
      dual_.faceVertexCount[i] = faceVertexCount[i];
    }
    delete[] faceVertexCount;
    dual_.faceVertexIndex = new uint[dualFaceCount + 1];
#pragma omp parallel for
    for (uint i=1; i <= dualFaceCount + 1; i++) {
      dual_.faceVertexIndex[i] = faceVertexIndex[i];
    }
    delete[] faceVertexIndex;
    uint dualFaceCellLength = dualFaceCount*2 + 2;
    dual_.faceCell = new uint[dualFaceCellLength];
#pragma omp parallel for
    for (uint i=0; i < dualFaceCellLength; i++) {
      dual_.faceCell[i] = faceCell[i];
    }
    delete[] faceCell;
    
    dual_.faceGhost.clear();
    dual_.faceGhost.resize(dualFaceCount + 1, false);
  }
  
  void DualMesher::CreateCellFaces() {
    std::cout << "Creating Cell -> Faces" << std::endl;
    dual_.cellFaceCount = CreateZeroArray<uint>(dual_.cellCount + 1);
    for (uint iFace = 1; iFace <= dual_.faceCount; iFace++) {
      dual_.cellFaceCount[dual_.faceCell[iFace*2]]++;
      dual_.cellFaceCount[dual_.faceCell[iFace*2+1]]++;
    }
    
    dual_.cellFaceIndex = CreateZeroArray<uint>(dual_.cellCount + 1);
    uint numFaceCells = 0;
    for (uint iCell = 1; iCell <= dual_.cellCount; iCell++) {
      dual_.cellFaceIndex[iCell] = numFaceCells;
      numFaceCells += dual_.cellFaceCount[iCell];
    }
    
    dual_.cellFace = CreateZeroArray<uint>(numFaceCells);
    for (uint iFace = 1; iFace <= dual_.faceCount; iFace++) {
      AddFaceToCell(dual_.faceCell[iFace*2], iFace);
      AddFaceToCell(dual_.faceCell[iFace*2+1], iFace);
    }
    std::cout << "Finished Creating Cell -> Faces" << std::endl;
  }
  
  void DualMesher::AddFaceToCell(uint& cell, uint& face) {
    if (cell == 0) {
      return;
    }
    uint index = dual_.cellFaceIndex[cell];
    while (dual_.cellFace[index] != 0) {
      index++;
    }
    dual_.cellFace[index] = face;
  }

  void DualMesher::CopyCoordinates() {
    std::cout << "Copying cell positions to dual vertex positions"<<std::endl;
    uint length = (dual_.vertexCount + 1)*3;
    dual_.vertexPosition = new double[length];
#pragma omp parallel for
    for (uint i=0; i < length; i++) {
      dual_.vertexPosition[i] = mesh_.cellPosition[i];
    }
    std::cout << "Copying vertex positions to dual cell positions"<<std::endl;
    length = (dual_.cellCount + 1)*3;
    dual_.cellPosition = new double[length];
#pragma omp parallel for
    for (uint i=0; i < length; i++) {
      dual_.cellPosition[i] = mesh_.vertexPosition[i];
    }
  }
  
  void DualMesher::CheckCollapsedFaces() {
    for (uint iFace = 1; iFace <= dual_.faceCount; iFace++) {
      uint* vertex = &dual_.faceVertex[dual_.faceVertexIndex[iFace]];
      uint max = dual_.faceVertexCount[iFace];
      if (vertex[0] == vertex[max - 1]) {
        std::cout << "Collapsed " << std::endl;
      }
      for (uint i=0; i < max - 1; i++) {
        if (vertex[i] == vertex[i+1]) {
          std::cout << "Collapsed " << std::endl;
        }
      }
    }
  }
  
}

