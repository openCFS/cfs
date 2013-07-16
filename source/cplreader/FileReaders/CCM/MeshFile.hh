/* 
 * File:   CCMMeshFile.h
 * Author: Mace
 *
 * Created on 31. Mai 2013, 23:56
 */

#ifndef CCMMESHFILE_H
#define	CCMMESHFILE_H

#include <vector>
#include <string>

#include "LibCCMIO/ccmio.h"
#include "LibCCMIO/ccmiocore.h"
#include "LibCCMIO/ccmiotypes.h"

#include "DataFile.hh"
#include "Utilities.hh"

namespace CCM {
  
  class MeshFile : public DataFile {
  public:
    MeshFile(std::string fileName);
    MeshFile(const MeshFile& orig);
    virtual ~MeshFile();
    
    void ReadFaces(std::vector<FaceData>* faceDatas);
    void ReadVertices(std::vector<VertexData>* vertexDatas);
    void ReadCells(std::vector<CellData>* cellDatas);
  protected:
    std::vector<CCMIOID> topologies_;
  
  };
  
}
#endif	/* CCMMESHFILE_H */

