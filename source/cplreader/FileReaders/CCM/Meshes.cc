

#include "Meshes.hh"

namespace CCM {
  
  bool TraverseEdge(uint primaryVertex, uint secondaryVertex, 
          uint cell, uint face, DualizedMesh* mesh,
          uint& nextCell, uint& nextFace, uint& nextVertexIndex, bool& isPrimary) {
    uint cellFaceCount = mesh->cellFaceCount[cell];
    uint faceIndex = mesh->cellFaceIndex[cell];
    uint maxFaceIndex = faceIndex + cellFaceCount;
    for (; faceIndex < maxFaceIndex; faceIndex++) {
      nextFace = mesh->cellFace[faceIndex];
      if (nextFace != face) {
        uint nextFaceVertexCount = mesh->faceVertexCount[nextFace];
        nextVertexIndex = mesh->faceVertexIndex[nextFace];
        uint maxNextFaceVertexIndex = nextVertexIndex + nextFaceVertexCount;
        for (; nextVertexIndex < maxNextFaceVertexIndex; nextVertexIndex++) {
          uint pVertex = mesh->faceVertex[nextVertexIndex];
          uint sVertex = 0;
          bool found = false;
          if (pVertex == primaryVertex) {
            sVertex = nextVertexIndex + 1 == maxNextFaceVertexIndex ? 
              mesh->faceVertex[mesh->faceVertexIndex[nextFace]] : mesh->faceVertex[nextVertexIndex + 1];
            if (sVertex == secondaryVertex) {
              found = true;
              isPrimary = true;
            }
          } else if (pVertex == secondaryVertex) {
            sVertex = nextVertexIndex + 1 == maxNextFaceVertexIndex ? 
              mesh->faceVertex[mesh->faceVertexIndex[nextFace]] : mesh->faceVertex[nextVertexIndex + 1];
            if (sVertex == primaryVertex) {
              found = true;
              isPrimary = false;
            }
          }
          if (found) {
            //nextCell = isPrimary ? mesh->faceCell[nextFace*2 + 1] : mesh->faceCell[nextFace*2];
            nextCell = mesh->faceCell[nextFace*2 + 1] == cell ?  mesh->faceCell[nextFace*2] : mesh->faceCell[nextFace*2 + 1];
            return true;
          }
        }
      }
    }
    return false;
  }
  
  void SplittableMesh::Clear() {
    delete[] cellFaceCount;
    cellFaceCount = NULL;
    delete[] cellFaceIndex;
    cellFaceIndex = NULL;
    delete[] cellFace;
    cellFace = NULL;
    std::vector<bool> emptyVectorCell;
    emptyVectorCell.swap(cellGhost);
    cellGhost.clear();
    
    delete[] faceVertexCount;
    faceVertexCount = NULL;
    delete[] faceVertexIndex;
    faceVertexIndex = NULL;
    delete[] faceVertex;
    faceVertex = NULL;
    std::vector<bool> emptyVectorFace;
    emptyVectorFace.swap(faceGhost);
    faceGhost.clear();
    
    delete[] vertexPosition;
    vertexPosition = NULL;
  }

  void DualizedMesh::Clear() {
    SplittableMesh::Clear();
    delete[] faceCell;
    faceCell = NULL;
    delete[] cellPosition;
    cellPosition = NULL;
  }
}
