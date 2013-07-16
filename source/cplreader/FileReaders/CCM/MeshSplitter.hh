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

namespace CCM {
  
  class MeshSplitter {
  public:
    MeshSplitter();
    MeshSplitter(const MeshSplitter& orig);
    virtual ~MeshSplitter();
    
    void SplitMesh(DualizedMesh& mesh, CPLMesh& cplMesh);
  private:
    void SplitCells();
    
    void InitTestSplit();
    void InitSplit();
    
    inline bool IsHexa(uint& cell, std::vector<uint>& faces4);
    inline bool IsHexa(uint& cell, std::vector<uint>& faces4, uint depth);
    inline bool IsWedge(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3);
    inline bool IsPrism(uint& cell, std::vector<uint>& faces4, std::vector<uint>& facesMore);
    
    inline bool CheckPrism(uint& cell, std::vector<uint>& faces4, 
        std::vector<uint>& facesX, std::vector<uint>& verticesOrdered);
    
    inline void CreatePyramid(uint& cell, std::vector<uint>& faces4, std::vector<uint>& faces3);
    inline void CreateTetra(uint& cell, std::vector<uint>& faces3);
    
    inline void SplitPolyhedron(uint& cell);
    
    inline bool ContainsFaceVertex(uint& face, uint& vertex);
    
    inline bool ShareFacesVertices(uint& faceA, uint& faceB);
    
    bool AddTetra(uint* nodes);
    bool AddPyramid(uint* nodes);
    bool AddWedge(uint* nodes);
    bool AddHexa(uint* nodes, bool countInvalid);
    void AddedElement();
    
    void VerboseSplitResults();
    
    void CopyVertices();
    
    void CopyMesh(CPLMesh& cplMesh);
    

    DualizedMesh mesh_;
    CPLMesh cplMesh_;
    
    uint rawTetras;
    uint rawPyramids;
    uint rawWedges;
    uint rawHexas;
    uint rawPrisms;
    uint rawPolyhedrons;
    
    uint invalidTetras;
    uint invalidPyramids;
    uint invalidWedges;
    uint invalidHexas;
   
    uint tetras;
    uint pyramids;
    uint wedges;
    uint hexas;
    uint sum;

    bool isTest;
  };
  
}

#endif	/* MESHSPLITTER_H */

