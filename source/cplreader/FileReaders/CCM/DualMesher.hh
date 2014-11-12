/* 
 * File:   DualMesher.h
 * Author: tz
 *
 * Created on 3. Juni 2013, 18:29
 */

//#ifndef DUALMESHER_H
//#define	DUALMESHER_H

#include <string>
#include <vector>

#include "Meshes.hh"

namespace CCM {
  
  class DualMesher {
  public:
    DualMesher();
    DualMesher(const DualMesher& orig);
    virtual ~DualMesher();
    
    void DualizeMesh(DualizableMesh& mesh, DualizedMesh& dual, bool copyCoordinates);
  private:
    void CreateFaces();
    void CreateCellFaces();
    void AddFaceToCell(uint& cell, uint& face);
    void CopyCoordinates();
    void CheckCollapsedFaces();

    DualizableMesh mesh_;
    DualizedMesh dual_;
  };
  
}

//#endif	/* DUALMESHER_H */

