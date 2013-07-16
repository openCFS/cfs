/* 
 * File:   MeshSplitter.cpp
 * Author: tz
 * 
 * Created on 3. Juni 2013, 10:57
 */

#include "MeshSplitter.hh"
#include "Meshes.hh"
#include "Geometrics.hh"

#include <vector>

#include "Domain/elem.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"



namespace CCM {
  
  MeshSplitter::MeshSplitter() {
    cplMesh_.TOPOLOGYDATA.clear();
  }
  
  MeshSplitter::MeshSplitter(const MeshSplitter& orig) {
  }
  
  MeshSplitter::~MeshSplitter() {
  }
  
  void MeshSplitter::SplitMesh(DualizedMesh& mesh, CPLMesh& cplMesh) {
    std::cout << "Splitting Mesh" << std::endl;
    mesh_ = mesh;

    InitTestSplit();
    SplitCells();
    
    InitSplit();
    SplitCells();
    VerboseSplitResults();
    CopyVertices();
    
    CopyMesh(cplMesh);
  }
  
  void MeshSplitter::InitTestSplit() {
    std::cout << "  Initializing Test Split" << std::endl;
    isTest = true;
    
    tetras = 0;
    pyramids = 0;
    wedges = 0;
    hexas = 0;
    sum = 0;

    invalidTetras = 0;
    invalidPyramids = 0;
    invalidWedges = 0;
    invalidHexas = 0;
    
    rawTetras = 0;
    rawPyramids = 0;
    rawWedges = 0;
    rawHexas = 0;
    rawPrisms = 0;
    rawPolyhedrons = 0;
  }
  
  void MeshSplitter::InitSplit() {
    std::cout << "  Initializing Mesh Split" << std::endl;
    isTest = false;
    cplMesh_.TOPOLOGYDATA.clear();
    cplMesh_.TOPOLOGYDATA.resize(8*sum);
    cplMesh_.elemTypes.clear();
    cplMesh_.elemTypes.resize(sum);
    cplMesh_.numRegions = 1;
    cplMesh_.regionElements.clear();
    cplMesh_.regionElements.resize(sum);
    cplMesh_.maxNumElemNodes = 8;
    
    tetras = 0;
    pyramids = 0;
    wedges = 0;
    hexas = 0;
    sum = 0;

    invalidTetras = 0;
    invalidPyramids = 0;
    invalidWedges = 0;
    invalidHexas = 0;
    
    rawTetras = 0;
    rawPyramids = 0;
    rawWedges = 0;
    rawHexas = 0;
    rawPrisms = 0;
    rawPolyhedrons = 0;
  }
  
  void MeshSplitter::SplitCells() {
    std::cout << "  Analyzing and splitting cells"<< std::endl;
    std::vector<uint> faces3;
    std::vector<uint> faces4;
    std::vector<uint> facesMore;
    
    
    for (uint iCell = 1; iCell <= mesh_.cellCount; iCell++) {
      if (!mesh_.cellGhost[iCell]) {
        faces3.clear();
        faces4.clear();
        facesMore.clear();
        bool isPolyhedron = true;
        
        uint miniCellFaceIndex = mesh_.cellFaceIndex[iCell];
        uint maxiCellFaceIndex = mesh_.cellFaceIndex[iCell] + mesh_.cellFaceCount[iCell];
        for (uint iCellFaceIndex = miniCellFaceIndex; iCellFaceIndex < maxiCellFaceIndex; iCellFaceIndex++) {
          uint face = mesh_.cellFace[iCellFaceIndex];
          if (mesh_.faceVertexCount[face] == 3) {
            faces3.push_back(face);
          } else if (mesh_.faceVertexCount[face] == 4) {
            faces4.push_back(face);
          } else {
            facesMore.push_back(face);
          }
        }
        
        if (faces4.size() == 6 && faces3.size() == 0 && facesMore.size() == 0) {
          if(IsHexa(iCell, faces4)) {
            rawHexas++;
            isPolyhedron = false;
          }
        } else if (faces4.size() == 3 && faces3.size() == 2 && facesMore.size() == 0) {
          if(IsWedge(iCell, faces4, faces3)) {
            rawWedges++;
            isPolyhedron = false;
          }
        } else if (faces4.size() == 1 && faces3.size() == 4 && facesMore.size() == 0) {
          CreatePyramid(iCell, faces4, faces3);
          rawPyramids++;
          isPolyhedron = false;
        } else if (faces4.size() == 0 && faces3.size() == 4 && facesMore.size() == 0) {
          CreateTetra(iCell, faces3);
          rawTetras++;
          isPolyhedron = false;
        } else if (faces3.size() == 0 && facesMore.size() == 2) {
          if (faces4.size() == mesh_.faceVertexCount[facesMore[0]]
                && faces4.size() == mesh_.faceVertexCount[facesMore[1]]) {
            if(IsPrism(iCell, faces4, facesMore)) {
              rawPrisms++;
              isPolyhedron = false;
            }
          }
        }
        if (isPolyhedron) {
          rawPolyhedrons++;
          SplitPolyhedron(iCell);
        }
      }
    }
  }
  
  inline bool MeshSplitter::IsHexa(uint& cell, std::vector<uint>& faces4) {
    return IsHexa(cell, faces4, (uint)0);
  }
 
  inline bool MeshSplitter::IsHexa(uint& cell, std::vector<uint>& faces4, uint depth) {
    std::vector<uint> topBottom;
    topBottom.push_back(faces4.back());
    faces4.pop_back();
    bool found = false;
    for (uint i=0; i < 5 && !found; i++) {
      if (!ShareFacesVertices(topBottom[0], faces4[i])) {
        topBottom.push_back(faces4[i]);
        faces4[i] = faces4[4];
        faces4.pop_back();
        found = true;
      }
    }
    std::vector<uint> verticesOrdered;
    bool isHexa = CheckPrism(cell, faces4, topBottom, verticesOrdered);
    if (isHexa) {
      if (!AddHexa(&verticesOrdered[0], (depth == 3)) && (depth < 3)) {
        double vectorTopBottom[3];
        double vectorTop[3];
        CreateVector(&mesh_.vertexPosition[verticesOrdered[0]*3], &mesh_.vertexPosition[verticesOrdered[1]*3], vectorTop);
        CreateVector(&mesh_.vertexPosition[verticesOrdered[0]*3], &mesh_.vertexPosition[verticesOrdered[4]*3], vectorTopBottom);
        if (Length(vectorTopBottom) < Length(vectorTop)) {
          IsPrism(cell, faces4, topBottom);
        } else {
          faces4.push_back(topBottom[0]);
          faces4.push_back(topBottom[1]);
          if (IsHexa(cell, faces4, depth + 1)) {
            return true;
          }
        }
      } else {
        SplitPolyhedron(cell);
      }
      return true;
    }
    return false;
  }
  
  inline bool MeshSplitter::IsWedge(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3) {
    std::vector<uint> verticesOrdered;
    bool isWedge = CheckPrism(cell, faces4, faces3, verticesOrdered);
    if (isWedge) {
      if (!AddWedge(&verticesOrdered[0])) {
        SplitPolyhedron(cell);
      }
      return true;
    }
    return false;
  }
  
  inline bool MeshSplitter::IsPrism(uint& cell, std::vector<uint>& faces4, std::vector<uint>& facesMore) {
    std::vector<uint> verticesOrdered;
    bool isPrism = CheckPrism(cell, faces4, facesMore, verticesOrdered);
    if (isPrism) {
      uint edges = verticesOrdered.size() / 2;
      uint nodes[6];
      nodes[0] = verticesOrdered[0];
      nodes[3] = verticesOrdered[edges];
      uint max = edges - 1;
      for (uint i=1; i < max; i++) {
        nodes[1] = verticesOrdered[i];
        nodes[2] = verticesOrdered[i+1];
        nodes[4] = verticesOrdered[i+edges];
        nodes[5] = verticesOrdered[i+edges+1];
        if (!AddWedge(nodes)) {
          uint tetranodes[4];
          tetranodes[0] = nodes[0];
          tetranodes[1] = nodes[1];
          tetranodes[2] = nodes[2];
          tetranodes[3] = nodes[3];
          AddTetra(tetranodes);
          tetranodes[0] = nodes[3];
          tetranodes[1] = nodes[4];
          tetranodes[2] = nodes[0];
          tetranodes[3] = nodes[5];
          AddTetra(tetranodes);
          tetranodes[0] = nodes[5];
          tetranodes[1] = nodes[4];
          tetranodes[2] = nodes[2];
          tetranodes[3] = nodes[0];
          AddTetra(tetranodes);
        }
      }
      return true;
    }
    return false;
  }

  inline void MeshSplitter::CreatePyramid(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3) {
    uint nodes[5];
    uint primFace = faces4[0];
    uint vertexIndex = mesh_.faceVertexIndex[primFace];
    uint maxIndex = vertexIndex + 4;
    uint i=0;
    for (; vertexIndex < maxIndex; vertexIndex++) {
      nodes[i] = mesh_.faceVertex[vertexIndex];
      i++;
    }
    if (mesh_.faceCell[primFace*2] == cell) {
      std::swap(nodes[0], nodes[3]);
      std::swap(nodes[1], nodes[2]);
    }
    
    i=0;
    vertexIndex = mesh_.faceVertexIndex[faces3[0]];
    maxIndex = vertexIndex + 3;
    for (;vertexIndex < maxIndex; vertexIndex++) {
      if (!ContainsFaceVertex(primFace, mesh_.faceVertex[vertexIndex])) {
        nodes[4] = mesh_.faceVertex[vertexIndex];
      }
    }
    if (!AddPyramid(nodes)) {
      SplitPolyhedron(cell);
    }
  }
  
  inline void MeshSplitter::CreateTetra(uint& cell, std::vector<uint>& faces3) {
    uint nodes[4];
    uint primFace = faces3[0];
    uint vertexIndex = mesh_.faceVertexIndex[primFace];
    uint maxIndex = vertexIndex + 3;
    uint i=0;
    for (; vertexIndex < maxIndex; vertexIndex++) {
      nodes[i] = mesh_.faceVertex[vertexIndex];
      i++;
    }
    if (mesh_.faceCell[primFace*2] == cell) {
      std::swap(nodes[0], nodes[2]);
    }
    
    i=0;
    vertexIndex = mesh_.faceVertexIndex[faces3[1]];
    maxIndex = vertexIndex + 3;
    for (;vertexIndex < maxIndex; vertexIndex++) {
      if (!ContainsFaceVertex(primFace, mesh_.faceVertex[vertexIndex])) {
        nodes[3] = mesh_.faceVertex[vertexIndex];
      }
    }
    AddTetra(nodes);
  }
  
  inline bool MeshSplitter::CheckPrism(uint& cell, std::vector<uint>& faces4, 
        std::vector<uint>& facesX, std::vector<uint>& verticesOrdered) {
    if (ShareFacesVertices(facesX[0],facesX[1])) {
      return false;
    }
    uint edges = faces4.size();
    verticesOrdered.resize(edges*2,(uint)0);
    
    
    // caring about top face;
    uint faceTop = facesX[0];
    std::vector<uint> topVertices;
    topVertices.resize(edges+1,(uint)0);
    std::vector<uint> topFaces;
    topFaces.resize(edges,(uint)0);
    uint faceTopVertexIndex = mesh_.faceVertexIndex[faceTop];
    for (uint i=0; i < edges; i++) {
      topVertices[i] = mesh_.faceVertex[faceTopVertexIndex+i];
    }
    bool topFacePrim = mesh_.faceCell[faceTop*2] == cell;
    if (topFacePrim) {
      for (uint i=0;i<edges/2; i++) {
        std::swap(topVertices[i],topVertices[edges-i-1]);
      }
    }
    topVertices[edges] = topVertices[0];
    for (uint i=0; i< edges; i++) {
      uint vertexA = topVertices[i];
      uint vertexB = topVertices[i+1];
      uint nextFace;
      uint dummy;
      bool bdummy;
      TraverseEdge(vertexA, vertexB, cell, faceTop, &mesh_, dummy, nextFace, dummy, bdummy); 
      topFaces[i] = nextFace;
    }
    for (uint i=0; i < edges; i++) {
      for (uint j=i+1; j < edges; j++) {
        if (topFaces[i] == topFaces[j]) {
          return false;
        }
      }
    }
    for (uint i=0; i < edges; i++) {
      verticesOrdered[i] = topVertices[i];
    }
    
    // caring about bottom face;
    uint faceBottom = facesX[1];
    std::vector<uint> bottomVertices;
    bottomVertices.resize(edges+1,(uint)0);
    std::vector<uint> bottomFaces;
    bottomFaces.resize(edges,(uint)0);
    uint faceBottomVertexIndex = mesh_.faceVertexIndex[faceBottom];
    for (uint i=0; i < edges; i++) {
      bottomVertices[i] = mesh_.faceVertex[faceBottomVertexIndex+i];
    }
    bool bottomFacePrim = mesh_.faceCell[faceBottom*2] == cell;
    if (!bottomFacePrim) {
      for (uint i=0;i<edges/2; i++) {
        std::swap(bottomVertices[i],bottomVertices[edges-i-1]);
      }
    }
    bottomVertices[edges] = bottomVertices[0];
    for (uint i=0; i< edges; i++) {
      uint vertexA = bottomVertices[i];
      uint vertexB = bottomVertices[i+1];
      uint nextFace;
      uint dummy;
      bool bdummy;
      TraverseEdge(vertexA, vertexB, cell, faceBottom, &mesh_, dummy, nextFace, dummy, bdummy); 
      bottomFaces[i] = nextFace;
    }
    for (uint i=0; i < edges; i++) {
      for (uint j=i+1; j < edges; j++) {
        if (bottomFaces[i] == bottomFaces[j]) {
          return false;
        }
      }
    }
    bool found = false;
    for (uint i=0; i< edges && !found; i++) {
      if (bottomFaces[i] == topFaces[0]) {
        found = true;
        for (uint j=0; j < edges; j++) {
          uint index = j + i >= edges ? j + i - edges : j + i;
          verticesOrdered[j + edges] = bottomVertices[index];
        }
      }
    }
    return true;
  }
  
  inline void MeshSplitter::SplitPolyhedron(uint& cell) {
    uint faceIndex = mesh_.cellFaceIndex[cell];
    uint maxFaceIndex = faceIndex + mesh_.cellFaceCount[cell];
    uint masterVertex = mesh_.faceVertex[mesh_.faceVertexIndex[mesh_.cellFace[faceIndex]]];
    faceIndex++;
    uint nodes[4];
    nodes[3] = masterVertex;
    
    for (;faceIndex < maxFaceIndex; faceIndex++) {
      uint face = mesh_.cellFace[faceIndex];
      uint sum = 0;
      if (!ContainsFaceVertex(face, masterVertex)) {
        bool master = mesh_.faceCell[face * 2] == cell;
        uint vertexIndex = mesh_.faceVertexIndex[face];
        uint maxVertexIndex = vertexIndex + mesh_.faceVertexCount[face] - 1;
        nodes[1] = mesh_.faceVertex[vertexIndex];
        vertexIndex++;
        for (;vertexIndex < maxVertexIndex; vertexIndex++) {
          if (master) {
            nodes[0] = mesh_.faceVertex[vertexIndex];
            nodes[2] = mesh_.faceVertex[vertexIndex+1];
          } else {
            nodes[2] = mesh_.faceVertex[vertexIndex];
            nodes[0] = mesh_.faceVertex[vertexIndex+1];
          }
          sum++;
          AddTetra(nodes);
        }
      }
    }
  }
  
  inline bool MeshSplitter::ContainsFaceVertex(uint& face, uint& vertex) {
    uint vertIndex = mesh_.faceVertexIndex[face];
    uint maxIndex = vertIndex + mesh_.faceVertexCount[face]; 
    for (;vertIndex < maxIndex; vertIndex++) {
      if (vertex == mesh_.faceVertex[vertIndex]) {
        return true;
      }
    }
    return false;
  }
  
  inline bool MeshSplitter::ShareFacesVertices(uint& faceA, uint& faceB) {
    uint vertIndexA = mesh_.faceVertexIndex[faceA];
    uint maxIndexA = vertIndexA + mesh_.faceVertexCount[faceA]; 
    
    uint startVertIndexB = mesh_.faceVertexIndex[faceB];
    uint maxIndexB = startVertIndexB + mesh_.faceVertexCount[faceB];
    uint vertIndexB = startVertIndexB;
    for (;vertIndexA < maxIndexA; vertIndexA++) {
      for (vertIndexB = startVertIndexB; vertIndexB < maxIndexB; vertIndexB++) {
        if (mesh_.faceVertex[vertIndexA] == mesh_.faceVertex[vertIndexB]) {
          return true;
        }
      }
    }
    return false;
  }
  
  bool MeshSplitter::AddTetra(uint* nodes) {
    double vectorA[3];
    double vectorB[3];
    double normVectorAB[3];
    double vectorC[3];
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[2]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[0]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[3]*3], vectorC);
    if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidTetras++;
        return false;
    }

    if (isTest) {
      tetras++;
      sum++;
      return true;
    }
    uint i;
    for (i=0; i < 4; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = nodes[i];
    }
    for (; i < 8; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = (uint)0;
    }
    cplMesh_.elemTypes[sum] = Elem::TET4;
    tetras++;
    AddedElement();
    return true;
  }
  
  bool MeshSplitter::AddPyramid(uint* nodes) {
    double vectorA[3];
    double vectorB[3];
    double normVectorAB[3];
    double vectorC[3];
    CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[1]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[3]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[4]*3], vectorC);
    if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidPyramids++;
        return false;
    }
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[2]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[0]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[4]*3], vectorC);
    if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidPyramids++;
        return false;
    }
    CreateVector(&mesh_.vertexPosition[nodes[2]*3], &mesh_.vertexPosition[nodes[3]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[2]*3], &mesh_.vertexPosition[nodes[1]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    CreateVector(&mesh_.vertexPosition[nodes[2]*3], &mesh_.vertexPosition[nodes[4]*3], vectorC);
    if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidPyramids++;
        return false;
    }
    CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[0]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[2]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[4]*3], vectorC);
    if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidPyramids++;
        return false;
    }

    if (isTest) {
      pyramids++;
      sum++;
      return true;
    }
    uint i;
    for (i=0; i < 5; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = nodes[i];
    }
    for (; i < 8; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = (uint)0;
    }
    cplMesh_.elemTypes[sum] = Elem::PYRA5;
    pyramids++;
    AddedElement();
    return true;
  }
  
  bool MeshSplitter::AddWedge(uint* nodes) {
    double vectorA[3];
    double vectorB[3];
    double normVectorAB[3];
    double vectorC[3];

    // face Upper
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[2]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[0]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=3; i < 6; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidWedges++;
        return false;
      }
    }
    // face Lower
    CreateVector(&mesh_.vertexPosition[nodes[4]*3], &mesh_.vertexPosition[nodes[3]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[4]*3], &mesh_.vertexPosition[nodes[5]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=0; i < 3; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[4]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        invalidWedges++;
        return false;
      }
    }


    if (isTest) {
      wedges++;
      sum++;
      return true;
    }
    uint i;
    for (i=0; i < 6; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = nodes[i];
    }
    for (; i < 8; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = (uint)0;
    }
    cplMesh_.elemTypes[sum] = Elem::WEDGE6;
    wedges++;
    AddedElement();
    return true;
  }
  
  bool MeshSplitter::AddHexa(uint* nodes, bool countInvalid) {
    double vectorA[3];
    double vectorB[3];
    double normVectorAB[3];
    double vectorC[3];

    // face Upper
    CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[1]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[3]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=4; i < 8; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[0]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[2]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[0]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=4; i < 8; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[2]*3], &mesh_.vertexPosition[nodes[3]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[2]*3], &mesh_.vertexPosition[nodes[1]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=4; i < 8; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[1]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[0]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[2]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=4; i < 8; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[3]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    // face Lower
    CreateVector(&mesh_.vertexPosition[nodes[4]*3], &mesh_.vertexPosition[nodes[7]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[4]*3], &mesh_.vertexPosition[nodes[5]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=0; i < 4; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[5]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[5]*3], &mesh_.vertexPosition[nodes[4]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[5]*3], &mesh_.vertexPosition[nodes[6]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=0; i < 4; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[5]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[6]*3], &mesh_.vertexPosition[nodes[5]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[6]*3], &mesh_.vertexPosition[nodes[7]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=0; i < 4; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[5]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }
    CreateVector(&mesh_.vertexPosition[nodes[7]*3], &mesh_.vertexPosition[nodes[6]*3], vectorA);
    CreateVector(&mesh_.vertexPosition[nodes[7]*3], &mesh_.vertexPosition[nodes[4]*3], vectorB);
    CrossProduct(vectorA, vectorB, normVectorAB);
    for (uint i=0; i < 4; i++) {
      CreateVector(&mesh_.vertexPosition[nodes[7]*3], &mesh_.vertexPosition[nodes[i]*3], vectorC);
      if (ScalarProduct(normVectorAB, vectorC) < 1E-15) {
        if (countInvalid) {
          invalidHexas++;
        }
        return false;
      }
    }


    if (isTest) {
      hexas++;
      sum++;
      return true;
    }
    for (uint i=0; i < 8; i++) {
      cplMesh_.TOPOLOGYDATA[sum*8+i] = nodes[i];
    }
    cplMesh_.elemTypes[sum] = Elem::HEXA8;
    hexas++;
    AddedElement();
    return true;
  }
  
  void MeshSplitter::AddedElement() {
    cplMesh_.regionElements[sum] = sum + 1;
    sum++;
  }
  
  void MeshSplitter::VerboseSplitResults() {
    std::cout << "    Initial Mesh Statistics" << std::endl;
    std::cout << "        Tetras:      " << rawTetras << std::endl;
    std::cout << "        Pyradmids:   " << rawPyramids << std::endl;
    std::cout << "        Wedges:      " << rawWedges << std::endl;
    std::cout << "        Hexas:       " << rawHexas << std::endl;
    std::cout << "        Prisms:      " << rawPrisms << std::endl;
    std::cout << "        Polyhedrons: " << rawPolyhedrons << std::endl;
    std::cout << "                     ---------" << std::endl;
    std::cout << "        SUM:         " << mesh_.cellCount << std::endl;
    std::cout << std::endl;
    std::cout << "    Splitted FEM Mesh Statistics" << std::endl;
    std::cout << "        Tetras:      " << tetras << std::endl;
    std::cout << "        Pyradmids:   " << pyramids << std::endl;
    std::cout << "        Wedges:      " << wedges << std::endl;
    std::cout << "        Hexas:       " << hexas << std::endl;
    std::cout << "                     ---------" << std::endl;
    std::cout << "        SUM:         " << sum << std::endl;
    std::cout << "        Invalid Tetras:   " << invalidTetras << std::endl;
    std::cout << "        Invalid Pyramids: " << invalidPyramids << std::endl;
    std::cout << "        Invalid Wedges:   " << invalidWedges << std::endl;
    std::cout << "        Invalid Hexas:    " << invalidHexas << std::endl;
    std::cout << std::endl;
  }
  
  void MeshSplitter::CopyVertices() {
    std::cout << "  Preparing Vertices " << mesh_.vertexCount << std::endl;
    cplMesh_.NODECOORD.clear();
    cplMesh_.NODECOORD.resize(mesh_.vertexCount*3);
    uint max = (mesh_.vertexCount+1)*3;
    for (uint i=3; i <= max; i++) {
      cplMesh_.NODECOORD[i-3] = mesh_.vertexPosition[i];
    }
  }
  
  void MeshSplitter::CopyMesh(CPLMesh& cplMesh) {
    cplMesh.maxNumElemNodes = cplMesh_.maxNumElemNodes;
    cplMesh.numRegions = cplMesh_.numRegions;
    cplMesh.NODECOORD.swap(cplMesh_.NODECOORD);
    cplMesh.regionElements.swap(cplMesh_.regionElements);
    cplMesh.TOPOLOGYDATA.swap(cplMesh_.TOPOLOGYDATA);
    cplMesh.elemTypes.swap(cplMesh_.elemTypes);
    cplMesh.numElemsPerRegion.clear();
    cplMesh.numElemsPerRegion.push_back(sum);
    cplMesh.numNodesPerRegion.clear();
    cplMesh.numNodesPerRegion.push_back(mesh_.vertexCount);
  }
  
}

