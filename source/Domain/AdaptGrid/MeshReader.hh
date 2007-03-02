// - C++ -
/***************************************************************************
    File        : MeshReader.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Jan 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef MESHREADER_H
#define MESHREADER_H


// System
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <cstring>
// Grid
#include "Stack.hh"
#include "Vertex.hh"
#include "Element.hh"
#include "Triangle.hh"
#include "Quadrangle.hh"
#include "Tetrahedron.hh"
#include "Hexahedron.hh"


namespace grd {


class MeshReader {
public:
  MeshReader() {}
  ~MeshReader() {}

  // Unit Square
  void createTriangleMesh();
  void createUnitSquareWithTris();
  void createUnitSquareWithQuad();
  // Unit cube
  void createUnitCubeWithTets();
  void createUnitCubeWithHexs();

  // Read Vertex and Element for
  // CAPA meshes
  void readCAPA(std::string& filename);

  // Read Vertex and Element for
  // GRUMMP meshes
  void readGRUMMP(const std::string& filename);
  
  // Read Vertex and Element for LSE-Mesh
  void readLSEMesh(const std::string& filename);

  // Access vertex and element lists
  std::vector<Vertex*>&  getVertexList()   { return vertex; }
  std::vector<Element*>& getElementList() { return element; }

private:
  std::vector<Vertex*>  vertex;
  std::vector<Element*> element;

  // CAPA meshes
  void readCAPAVertex(std::ifstream& theFile,std::vector<int>& vertexId,std::vector<double*>& vertexValues);
  void readCAPAElement(std::ifstream& theFile,std::vector<int>& elemId,std::vector<int*>& elemTopology);
  void setCAPAVertexAndElement(std::vector<int>& vertexId,std::vector<double*>& vertexValues,
                               std::vector<int>& elemId,std::vector<int*>& elemTopology);

  // Private read functions
  void setVertexAndElement(std::vector<int*>& elem,std::vector<double*>& vtr);
};

} // namespace grd

#endif // MESHREADER_H
