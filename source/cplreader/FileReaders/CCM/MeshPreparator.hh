/* 
 * File:   MeshPreparator.hh
 * Author: Matthias Tautz
 *
 * Created on 2. Juni 2013, 15:15
 */

#ifndef MESHPREPARATOR_H
#define	MESHPREPARATOR_H

#include <string>
#include <vector>

#include "MeshFile.hh"
#include "Utilities.hh"
#include "Meshes.hh"



namespace CCM {
  
  class MeshPreparator {
  public:
    MeshPreparator();
    MeshPreparator(const MeshPreparator& orig);
    virtual ~MeshPreparator();
    
    void PrepareMesh(MeshFile* meshFile, DualizableMesh& mesh, std::vector<ConsecutiveMap>& valueConMaps);
  private:
    void ReadMesh();
    
    void RemapRawData();
    Map GetMap(MappedData* mappedData);
    void ClearMaps();
    
    void CreateVertices();
    void ClearVertices();
    
    void DetermineFaceCounts();
    void CreateFaceToVertex();
    void CreateFaceToCell();
    void ClearFaces();
    
    void CreateCellToFace();
    inline void AddFaceToCell(uint& cell, uint& face);
    
    void MarkGhosts();
    void CreateGhostFaces();
    
    void ReadCellPositions();
    
    void PrintExportRecommends();
    
    DualizableMesh m_;
    DualizableMesh* mesh_;
    MeshFile* meshFile_;
    
    // raw mesh data;
    std::vector<Map> maps_;
    std::vector<FaceData> faces_;
    std::vector<VertexData> vertices_;
    
    std::vector<ConsecutiveMap> cellConMaps_;
    std::vector<ConsecutiveMap> cellConMapsForFaces_;
    std::vector<ConsecutiveMap> vertexConMaps_;
  };
  
}

#endif	/* MESHPREPARATOR_H */

