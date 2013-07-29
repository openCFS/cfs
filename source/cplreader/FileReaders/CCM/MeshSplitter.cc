/* 
 * File:   MeshSplitter.cpp
 * Author: Matthias Tautz
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

    isTest = true;
    ClearCounters();
    SplitCells();
    
    isTest = false;
    InitSplittedMesh();
    ClearCounters();
    SplitCells();
    
    VerboseSplitResults();
    
    CopyVertices();
    CopyMesh(cplMesh);
  }
  
  void MeshSplitter::SplitCells() {
    if (isTest) {
      std::cout << "  Analyzing and splitting cells (testing)"<< std::endl;
    } else {
      std::cout << "  Analyzing and splitting cells "<< std::endl;
    }
    std::cout << std::endl;
    
    std::vector<uint> faces3;
    std::vector<uint> faces4;
    std::vector<uint> facesMore;
    for (uint iCell = 1; iCell <= mesh_.cellCount; iCell++) {
      if (!mesh_.cellGhost[iCell]) {
        faces3.clear();
        faces4.clear();
        facesMore.clear();
        
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
          RunHexaTreatment(iCell, faces4);
        } else if (faces4.size() == 3 && faces3.size() == 2 && facesMore.size() == 0) {
          RunWedgeTreatment(iCell, faces4, faces3);
        } else if (faces4.size() == 1 && faces3.size() == 4 && facesMore.size() == 0) {
          RunPyramidTreatment(iCell, faces4, faces3);
        } else if (faces4.size() == 0 && faces3.size() == 4 && facesMore.size() == 0) {
          RunTetraTreatment(iCell, faces3);
        } else if (faces3.size() == 0 && facesMore.size() == 2) {
          if (faces4.size() == mesh_.faceVertexCount[facesMore[0]]
                && faces4.size() == mesh_.faceVertexCount[facesMore[1]]) {
            RunPrismTreatment(iCell, faces4, facesMore);
          } else {
            RunPolyTreatment(iCell);
          }
        } else {
          RunPolyTreatment(iCell);
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
  
  inline void PrintValidInvalidSplit(uint& valid, uint& invalid) {
    std::cout << "        Valid:       " << valid << std::endl;
    std::cout << "        Invalid:     " << invalid << std::endl << std::endl;
  }
  
  void MeshSplitter::VerboseSplitResults() {
    std::cout << "    Topological Mesh Statistics" << std::endl;
    std::cout << "        Polyhedrons: " << topologicalPolyhedrons << std::endl;
    std::cout << "        Prisms:      " << topologicalPrisms << std::endl;
    std::cout << "        Hexas:       " << topologicalHexas << std::endl;
    std::cout << "        Wedges:      " << topologicalWedges << std::endl;
    std::cout << "        Pyradmids:   " << topologicalPyramids << std::endl;
    std::cout << "        Tetras:      " << topologicalTetras << std::endl;
    std::cout << "                     ---------" << std::endl;
    std::cout << "        SUM:         " << mesh_.cellCount << std::endl;
    std::cout << std::endl;
    
    std::cout << "    Invalid Topological FEM Cells" << std::endl;
    std::cout << "        Hexas:       " << invalidHexas << std::endl;
    std::cout << "        Wedges:      " << invalidWedges << std::endl;
    std::cout << "        Pyradmids:   " << invalidPyramids << std::endl;
    std::cout << "        Tetras:      " << invalidTetras << std::endl;
    std::cout << std::endl;
    
    if (topologicalPolyhedrons > 0 || topologicalPrisms > 0 || invalidWedgesFromPrisms > 0 || invalidHexas > 0 || invalidWedgesFromHexas > 0 ||
      invalidWedges > 0 || invalidPyramids > 0) {
      std::cout << "    Split Statistics:" << std::endl;
      if (topologicalPolyhedrons > 0) {
        std::cout << "    Polyhedrons -> Tetras" << std::endl;
	PrintValidInvalidSplit(tetrasFromPolyhedrons, invalidTetrasFromPolyhedrons);
      }
      if (topologicalPrisms > 0) {
        std::cout << "    Prisms -> Wedges" << std::endl;
	PrintValidInvalidSplit(wedgesFromPrisms, invalidWedgesFromPrisms);
      }
      if (invalidWedgesFromPrisms > 0) {
        std::cout << "    Prisms -> Invalid Wedges -> Tetras" << std::endl;
	PrintValidInvalidSplit(tetrasFromPrisms, invalidTetrasFromPrisms);
      }
      if (invalidHexas > 0) {
        std::cout << "    Invalid Hexas -> Wedges" << std::endl;
	PrintValidInvalidSplit(wedgesFromHexas, invalidWedgesFromHexas);
      }
      if (invalidWedgesFromHexas > 0) {
        std::cout << "    Invalid Hexas -> Invalid Wedges -> Tetras" << std::endl;
	PrintValidInvalidSplit(tetrasFromHexas, invalidTetrasFromHexas);
      }
      if (invalidWedges > 0) {
        std::cout << "    Invalid Wedges -> Tetras" << std::endl;
	PrintValidInvalidSplit(tetrasFromWedges, invalidTetrasFromWedges);
      }
      if (invalidPyramids > 0) {
        std::cout << "    Invalid Pyramids -> Tetras" << std::endl;
	PrintValidInvalidSplit(tetrasFromPyramids, invalidTetrasFromPyramids);
      }
    }
    
    std::cout << "    Splitted FEM Mesh Statistics" << std::endl;
    std::cout << "        Hexas:       " << hexas << std::endl;
    std::cout << "        Wedges:      " << wedges << std::endl;
    std::cout << "        Pyradmids:   " << pyramids << std::endl;
    std::cout << "        Tetras:      " << tetras << std::endl;
    std::cout << "                     ---------" << std::endl;
    std::cout << "        SUM:         " << sum << std::endl << std::endl;
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
  
  inline void MeshSplitter::NodeNumsToNodeCoords(uint* nodes, double** coords, uint nodeCount) {
    for (uint i = 0; i < nodeCount; i++) {
      coords[i] = &mesh_.vertexPosition[nodes[i]*3];
    }
  }

  inline bool MeshSplitter::IsValidTetra(uint* nodes) {
    double* coords[12];
    NodeNumsToNodeCoords(nodes, coords, 4);
    return AreTriangleNormal<0, 1, 2, 3, -1, -1, -1>(coords);
  }
  
  inline bool MeshSplitter::IsValidPyramid(uint* nodes) {
    double* coords[15];
    NodeNumsToNodeCoords(nodes, coords, 5);
    return AreQuadNormal<0, 1, 2, 3, 4, -1, -1, -1>(coords);
  }
  
  inline bool MeshSplitter::IsValidWedge(uint* nodes) {
    double* coords[18];
    NodeNumsToNodeCoords(nodes, coords, 6);
    return AreTriangleNormal<0, 1, 2, 3, 4, 5, -1>(coords) && 
          AreTriangleNormal<5, 4, 3, 0, 1, 2, -1>(coords) &&
          AreQuadNormal<3, 4, 1, 0, 2, 5, -1, -1>(coords) &&
          AreQuadNormal<4, 5, 2, 1, 0, 3, -1, -1>(coords) &&
          AreQuadNormal<5, 3, 0, 2, 1, 4, -1, -1>(coords);
  }
  
  inline bool MeshSplitter::IsValidHexa(uint* nodes) {
    double* coords[24];
    NodeNumsToNodeCoords(nodes, coords, 8);
    return AreQuadNormal<0, 1, 2, 3, 4, 5, 6, 7>(coords) && 
          AreQuadNormal<7, 6, 5, 4, 0, 1, 2, 3>(coords) &&
          AreQuadNormal<0, 3, 7, 4, 1, 2, 5, 6>(coords) &&
          AreQuadNormal<5, 6, 2, 1, 0, 3, 7, 4>(coords) &&
          AreQuadNormal<1, 0, 4, 5, 2, 3, 6, 7>(coords) &&
          AreQuadNormal<6, 7, 3, 2, 0, 1, 4, 5>(coords);
  }
  
  inline void MeshSplitter::BuildTetra(uint* nodes) {
    tetras++;
    BuildCell(nodes, 4, CoupledField::Elem::TET4);
  }
  
  inline void MeshSplitter::BuildPyramid(uint* nodes) {
    pyramids++;
    BuildCell(nodes, 5, CoupledField::Elem::PYRA5);
  }
  
  inline void MeshSplitter::BuildWedge(uint* nodes) {
    wedges++;
    BuildCell(nodes, 6, CoupledField::Elem::WEDGE6);
  }
  
  inline void MeshSplitter::BuildHexa(uint* nodes) {
    hexas++;
    BuildCell(nodes, 8, CoupledField::Elem::HEXA8);
  }
  
  inline void MeshSplitter::BuildCell(uint* nodes, uint nodeCount, uint elemType) {
    if (isTest) {
      sum++;
      return;
    }
    uint i=0;
    for (i=0; i < nodeCount; i++) {
      cplMesh_.TOPOLOGYDATA[sum * maxNumElemNodes + i] = nodes[i];
    }
    for (; i < maxNumElemNodes; i++) {
      cplMesh_.TOPOLOGYDATA[sum * maxNumElemNodes + i] = (uint)0;
    }
    cplMesh_.regionElements[sum] = sum + 1;
    cplMesh_.elemTypes[sum] = elemType;
    sum++;
  }
  
  inline bool MeshSplitter::IsPrismaticCell(uint& cell, std::vector<uint>& faces4, 
        std::vector<uint>& facesTopBottom, std::vector<uint>& verticesOrdered) {
    if (ShareFacesVertices(facesTopBottom[0],facesTopBottom[1])) {
      return false;
    }
    uint edges = faces4.size();
    verticesOrdered.resize(edges*2,(uint)0);
    
    
    // caring about top face;
    uint faceTop = facesTopBottom[0];
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
    uint faceBottom = facesTopBottom[1];
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
  
  inline bool MeshSplitter::SortHexaFaces(std::vector<uint>& faces4, 
        std::vector<uint>& facesTopBottom, std::vector<uint>& facesSides) {
    facesTopBottom.clear();
    facesSides.clear();
    facesTopBottom.push_back(faces4[0]);
    bool found = false;
    for (uint i = 1; i < 6; i++) {
      if (!ShareFacesVertices(facesTopBottom[0], faces4[i])) {
        facesTopBottom.push_back(faces4[i]);
        found = true;
      } else {
        facesSides.push_back(faces4[i]);
      }
    }
    return found;
  }
  
  inline void MeshSplitter::PrismaticNodesToWedges(std::vector<uint>& verticesOrdered, uint& wedgeCount, uint& invalidWedgeCount,
        uint& tetraCount, uint& invalidTetraCount) {
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
      if (IsValidWedge(nodes)) {
        BuildWedge(nodes);
        wedgeCount++;
      } else {
        invalidWedgeCount++;
        uint tetranodes[4];
        
        tetranodes[0] = nodes[0];
        tetranodes[1] = nodes[1];
        tetranodes[2] = nodes[2];
        tetranodes[3] = nodes[3];
        TryBuildTetra(tetranodes, tetraCount, invalidTetraCount);
        
        tetranodes[0] = nodes[3];
        tetranodes[1] = nodes[4];
        tetranodes[2] = nodes[0];
        tetranodes[3] = nodes[5];
        TryBuildTetra(tetranodes, tetraCount, invalidTetraCount);
        
        tetranodes[0] = nodes[5];
        tetranodes[1] = nodes[4];
        tetranodes[2] = nodes[2];
        tetranodes[3] = nodes[0];
        TryBuildTetra(tetranodes, tetraCount, invalidTetraCount);
      }
    }
  }
  
  inline void MeshSplitter::TryBuildTetra(uint* tetraNodes, uint& tetraCount, uint& invalidTetraCount) {
    if (IsValidTetra(tetraNodes)) {
      BuildTetra(tetraNodes);
      tetraCount++;
    } else {
      invalidTetraCount++;
    }
  }
  
  inline void MeshSplitter::SortHexaNodesForAspectRatio(uint* nodes) {
    double LengthA = CalcNodeDistance<0, 4>(nodes) + CalcNodeDistance<1, 5>(nodes) + CalcNodeDistance<2, 6>(nodes) + CalcNodeDistance<3, 7>(nodes);
    double LengthB = CalcNodeDistance<1, 2>(nodes) + CalcNodeDistance<5, 6>(nodes) + CalcNodeDistance<4, 7>(nodes) + CalcNodeDistance<0, 3>(nodes);
    double LengthC = CalcNodeDistance<0, 1>(nodes) + CalcNodeDistance<2, 3>(nodes) + CalcNodeDistance<4, 5>(nodes) + CalcNodeDistance<6, 7>(nodes);
    
    if (LengthC < LengthA && LengthC < LengthB) {
      std::swap(nodes[0], nodes[3]);
      std::swap(nodes[1], nodes[7]);
      std::swap(nodes[2], nodes[4]);
      std::swap(nodes[5], nodes[6]);
    } else if (LengthB < LengthA && LengthB < LengthC) {
      std::swap(nodes[0], nodes[1]);
      std::swap(nodes[2], nodes[4]);
      std::swap(nodes[3], nodes[5]);
      std::swap(nodes[6], nodes[7]);
    }
  }

  inline void MeshSplitter::RunPolyTreatment(uint& cell) {
    topologicalPolyhedrons++;
    
    uint faceIndex = mesh_.cellFaceIndex[cell];
    uint maxFaceIndex = faceIndex + mesh_.cellFaceCount[cell];
    uint masterVertex = mesh_.faceVertex[mesh_.faceVertexIndex[mesh_.cellFace[faceIndex]]];
    faceIndex++;
    uint nodes[4];
    nodes[3] = masterVertex;
    
    for (;faceIndex < maxFaceIndex; faceIndex++) {
      uint face = mesh_.cellFace[faceIndex];
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
          TryBuildTetra(&nodes[0], tetrasFromPolyhedrons, invalidTetrasFromPolyhedrons);
        }
      }
    }
  }

  inline void MeshSplitter::RunPrismTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& facesMore) {
    std::vector<uint> verticesOrdered;
    bool success = IsPrismaticCell(cell, faces4, facesMore, verticesOrdered);
    if (!success) {
      RunPolyTreatment(cell);
      return;
    }
    topologicalPrisms++;
    PrismaticNodesToWedges(verticesOrdered, wedgesFromPrisms, invalidWedgesFromPrisms,
              tetrasFromPrisms, invalidTetrasFromPrisms);
  }
  
  inline void MeshSplitter::RunHexaTreatment(uint& cell, std::vector<uint>& faces4) {
    std::vector<uint> facesTopBottom;
    std::vector<uint> facesSides;
    bool success = true;
    success = SortHexaFaces(faces4, facesTopBottom, facesSides);
    if (!success) {
      RunPolyTreatment(cell);
      return;
    }
    std::vector<uint> verticesOrdered;
    success = IsPrismaticCell(cell, facesSides, facesTopBottom, verticesOrdered);
    if (!success) {
      RunPolyTreatment(cell);
      return;
    }
    topologicalHexas++;
    success = IsValidHexa(&verticesOrdered[0]);
    if (!success) {
      invalidHexas++;
      SortHexaNodesForAspectRatio(&verticesOrdered[0]);
      PrismaticNodesToWedges(verticesOrdered, wedgesFromHexas, invalidWedgesFromHexas,
            tetrasFromHexas, invalidTetrasFromHexas);
      return;
    }
    BuildHexa(&verticesOrdered[0]);
  }
  
  inline void MeshSplitter::RunWedgeTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3) {
    std::vector<uint> verticesOrdered;
    bool success = IsPrismaticCell(cell, faces4, faces3, verticesOrdered);
    if (!success) {
      RunPolyTreatment(cell);
      return;
    }
    topologicalWedges++;
    uint dummy = (uint)0;
    PrismaticNodesToWedges(verticesOrdered, dummy, invalidWedges, tetrasFromWedges, invalidTetrasFromWedges);
  }
  
  inline void MeshSplitter::RunPyramidTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3) {
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
    
    topologicalPyramids++;
    if (IsValidPyramid(nodes)) {
      BuildPyramid(nodes);
    } else {
      std::swap(nodes[3], nodes[4]);
      TryBuildTetra(nodes, tetrasFromPyramids, invalidTetrasFromPyramids);
      std::swap(nodes[1], nodes[4]);
      std::swap(nodes[1], nodes[2]);
      TryBuildTetra(nodes, tetrasFromPyramids, invalidTetrasFromPyramids);
    }
  }

  inline void MeshSplitter::RunTetraTreatment(uint& cell, std::vector<uint>& faces3) {
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
    
    topologicalTetras++;
    uint dummy = (uint)0;
    TryBuildTetra(nodes, dummy, invalidTetras);
  }
  
  inline void MeshSplitter::ClearCounters() {
    topologicalTetras = 0;
    topologicalPyramids = 0;
    topologicalWedges = 0;
    topologicalHexas = 0;
    topologicalPrisms = 0;
    topologicalPolyhedrons = 0;
    
    invalidTetras = 0;
    invalidPyramids = 0;
    invalidWedges = 0;
    invalidHexas = 0;
    
    tetrasFromPolyhedrons = 0;
    invalidTetrasFromPolyhedrons = 0;
    
    wedgesFromPrisms = 0;
    invalidWedgesFromPrisms = 0;
    tetrasFromPrisms = 0;
    invalidTetrasFromPrisms = 0;
   
    wedgesFromHexas = 0;
    invalidWedgesFromHexas = 0;
    tetrasFromHexas = 0;
    invalidTetrasFromHexas = 0;
    
    tetrasFromWedges = 0;
    invalidTetrasFromWedges = 0;
    
    tetrasFromPyramids = 0;
    invalidTetrasFromPyramids = 0;
    
    tetras = 0;
    pyramids = 0;
    wedges = 0;
    hexas = 0;
    sum = 0;
  }
 
  inline void MeshSplitter::InitSplittedMesh() {
    if (hexas > 0) {
      maxNumElemNodes = 8;
    } else if (wedges > 0) {
      maxNumElemNodes = 6;
    } else if (pyramids > 0) {
      maxNumElemNodes = 5;
    } else if (tetras > 0) {
      maxNumElemNodes = 4;
    }
    
    cplMesh_.TOPOLOGYDATA.clear();
    cplMesh_.TOPOLOGYDATA.resize(maxNumElemNodes*sum);
    cplMesh_.elemTypes.clear();
    cplMesh_.elemTypes.resize(sum);
    cplMesh_.numRegions = 1;
    cplMesh_.regionElements.clear();
    cplMesh_.regionElements.resize(sum);
    cplMesh_.maxNumElemNodes = maxNumElemNodes;
  }
 
}

