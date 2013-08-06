/* 
 * File:   Meshes.h
 * Author: tz
 *
 * Created on 31. Mai 2013, 10:17
 */

#ifndef MESHES_H
#define	MESHES_H

#include <vector>
#include <string>
#include "General/environment.hh"

namespace CCM {
  
  //typedef unsigned int uint;
  
  struct SplittableMesh {
    std::string name;
    
    uint cellCount;
    uint* cellFaceCount;
    uint* cellFaceIndex;
    uint* cellFace;
    uint cellGhostCount;
    std::vector<bool> cellGhost;
            
    uint faceCount;
    uint* faceVertexCount;
    uint* faceVertexIndex;
    uint* faceVertex;
    uint faceGhostCount;
    std::vector<bool> faceGhost;
    
    uint vertexCount;
    double* vertexPosition;
    
    virtual void Clear();
    
  };
  
  struct DualizedMesh : SplittableMesh {
    uint cellFaceLength;
    double* cellPosition;
    
    uint* faceCell;
    uint faceVertexLength;
    virtual void Clear();
  };
  
  struct DualizableMesh : DualizedMesh {
    
  };
  
  struct ConsecutiveMap {
    uint id;
    std::string label;
    uint size;
    uint startValue;
  };
  
  bool TraverseEdge(uint primaryVertex, uint secondaryVertex, 
          uint cell, uint face, DualizedMesh* mesh,
          uint& nextCell, uint& nextFace, uint& nextVertexIndex, bool& isPrimary);
  
  struct CPLMesh {
    std::vector<double> NODECOORD;
    std::vector<uint> TOPOLOGYDATA;
    uint maxNumElemNodes;
    std::vector<uint> elemTypes;
    std::vector<uint> regionElements;
    std::vector<uint> numNodesPerRegion;
    std::vector<uint> numElemsPerRegion;
    uint numRegions;
  };
  
}



#endif	/* MESHES_H */

