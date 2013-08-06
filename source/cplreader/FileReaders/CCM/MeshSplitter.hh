/* 
 * File:   MeshSplitter.h
 * Author: tz
 *
 * Created on 3. Juni 2013, 10:57
 */

#ifndef MESHSPLITTER_H
#define	MESHSPLITTER_H

#include <vector>
#include "Meshes.hh"
#include "Domain/elem.hh"
#include "Geometrics.hh"


namespace CCM {
  
  class MeshSplitter {
  public:
    MeshSplitter();
    MeshSplitter(const MeshSplitter& orig);
    virtual ~MeshSplitter();
    
    void SplitMesh(DualizedMesh& mesh, CPLMesh& cplMesh);
  private:
    void SplitCells();
    
    inline bool ContainsFaceVertex(uint& face, uint& vertex);
    
    inline bool ShareFacesVertices(uint& faceA, uint& faceB);
    
    void VerboseSplitResults();
    
    void CopyVertices();
    
    void CopyMesh(CPLMesh& cplMesh);
    
    inline void NodeNumsToNodeCoords(uint* nodes, double** coords, uint nodeCount);
    
    template<int faceIdx, int nodeIdx>
    inline bool IsNodeNormal(double** coords, double* normABC) {
      double vectorFaceNode[3];
      CreateVector<double>(coords[faceIdx], coords[nodeIdx], vectorFaceNode);
      return !(ScalarProduct(vectorFaceNode, normABC) < 1E-15);
    }
    
    template<int faceIdxA, int faceIdxB, int faceIdxC, int nodeIdxA, int nodeIdxB, int nodeIdxC, int nodeIdxD>
    inline bool AreTriangleNormal(double** coords) {
      double vectorAB[3];
      double vectorBC[3];
      double normABC[3];
      CreateVector(coords[faceIdxA], coords[faceIdxB], vectorAB);
      CreateVector(coords[faceIdxB], coords[faceIdxC], vectorBC);
      CrossProduct(vectorAB, vectorBC, normABC);
      return (nodeIdxA < 0 || IsNodeNormal<faceIdxA, nodeIdxA>(coords, &normABC[0])) &&
            (nodeIdxB < 0 || IsNodeNormal<faceIdxA, nodeIdxB>(coords, &normABC[0])) &&
            (nodeIdxC < 0 || IsNodeNormal<faceIdxA, nodeIdxC>(coords, &normABC[0])) &&
            (nodeIdxD < 0 || IsNodeNormal<faceIdxA, nodeIdxD>(coords, &normABC[0]));
    }
    
    template<int faceIdxA, int faceIdxB, int faceIdxC, int faceIdxD, int nodeIdxA, int nodeIdxB, int nodeIdxC, int nodeIdxD>
    inline bool AreQuadNormal(double** coords) {
      return AreTriangleNormal<faceIdxA, faceIdxB, faceIdxC, nodeIdxA, nodeIdxB, nodeIdxC, nodeIdxD>(coords) &&
            AreTriangleNormal<faceIdxB, faceIdxC, faceIdxD, nodeIdxA, nodeIdxB, nodeIdxC, nodeIdxD>(coords) &&
            AreTriangleNormal<faceIdxC, faceIdxD, faceIdxA, nodeIdxA, nodeIdxB, nodeIdxC, nodeIdxD>(coords) &&
            AreTriangleNormal<faceIdxD, faceIdxA, faceIdxB, nodeIdxA, nodeIdxB, nodeIdxC, nodeIdxD>(coords);
    }
    
    template<uint nodeIdxA, uint nodeIdxB>
    inline double CalcNodeDistance(uint* nodes) {
      double vector[3];
      CreateVector(&mesh_.vertexPosition[nodes[nodeIdxA]*3], &mesh_.vertexPosition[nodes[nodeIdxB]*3], vector);
      return Length(vector);
    }
    
    inline bool IsValidTetra(uint* nodes);
    inline bool IsValidPyramid(uint* nodes);
    inline bool IsValidWedge(uint* nodes);
    inline bool IsValidHexa(uint* nodes);
    
    inline void BuildTetra(uint* nodes);
    inline void BuildPyramid(uint* nodes);
    inline void BuildWedge(uint* nodes);
    inline void BuildHexa(uint* nodes);
    
    inline void BuildCell(uint* nodes, uint nodeCount, uint elemType);
    
    inline bool IsPrismaticCell(uint& cell, std::vector<uint>& faces4, 
        std::vector<uint>& facesTopBottom, std::vector<uint>& verticesOrdered);
    
    inline bool SortHexaFaces(std::vector<uint>& faces4, std::vector<uint>& facesTopBottom, std::vector<uint>& facesSides);
    
    inline void PrismaticNodesToWedges(std::vector<uint>& verticesOrdered, uint& wedgeCount, uint& invalidWedgeCount,
        uint& tetraCount, uint& invalidTetraCount);

    inline void TryBuildTetra(uint* tetraNodes, uint& tetraCount, uint& invalidTetraCount);
    
    inline void SortHexaNodesForAspectRatio(uint* nodes);
    
    inline void RunPolyTreatment(uint& cell);
    inline void RunPrismTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& facesMore);
    inline void RunHexaTreatment(uint& cell, std::vector<uint>& faces4);
    inline void RunWedgeTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3);
    inline void RunPyramidTreatment(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3);
    inline void RunTetraTreatment(uint& cell, std::vector<uint>& faces3);
    
    
    inline void ClearCounters();
    
    inline void InitSplittedMesh();
    
    DualizedMesh mesh_;
    CPLMesh cplMesh_;
    
    uint topologicalTetras;
    uint topologicalPyramids;
    uint topologicalWedges;
    uint topologicalHexas;
    uint topologicalPrisms;
    uint topologicalPolyhedrons;
    
    uint invalidTetras;
    uint invalidPyramids;
    uint invalidWedges;
    uint invalidHexas;
    
    uint tetrasFromPolyhedrons;
    uint invalidTetrasFromPolyhedrons;
    
    uint wedgesFromPrisms;
    uint invalidWedgesFromPrisms;
    uint tetrasFromPrisms;
    uint invalidTetrasFromPrisms;
   
    uint wedgesFromHexas;
    uint invalidWedgesFromHexas;
    uint tetrasFromHexas;
    uint invalidTetrasFromHexas;
    
    uint tetrasFromWedges;
    uint invalidTetrasFromWedges;
    
    uint tetrasFromPyramids;
    uint invalidTetrasFromPyramids;
    
    
    uint tetras;
    uint pyramids;
    uint wedges;
    uint hexas;
    uint sum;
    
    uint maxNumElemNodes;

    bool isTest;
  };
  
}

#endif	/* MESHSPLITTER_H */

