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
#include "Stack.h"
#include "Vertex.h"
#include "Element.h"
#include "Triangle.h"
#include "Quadrangle.h"
#include "Tetrahedron.h"


namespace grd {


class MeshReader {
public:
  MeshReader() {}
  ~MeshReader() {}

  // Read Vertex and Element for
  // SRG (Hierarchical Tri-Mesh
  void readSRG(string& filename);

  // Read Vertex and Element for
  // CAPA meshes
  void readCAPA(string& filename);

  // Read Vertex and Element for
  // GRUMMP meshes
  void readGRUMMP(string& filename);

  // Access vertex and element lists
  vector<Vertex*>&  getVertexList()   { return vertex; }
  vector<Element*>& getElementList() { return element; }
  Stack<double**>&  getVertexStack()  { return mlVertex; }
  Stack<int**>&     getElementStack() { return mlElement; }

private:
  vector<Vertex*>  vertex;
  vector<Element*> element;
  Stack<double**>  mlVertex;
  Stack<int**>     mlElement;

  // SRG meshes
  void readCoarseMeshSRG(ifstream& theFile);
  void readMeshLevelSRG(ifstream& theFile,int level);
  void buildMultilevelBoundary();

  // CAPA meshes
  void readCAPAVertex(ifstream& theFile,vector<int>& vertexId,vector<double*>& vertexValues);
  void readCAPAElement(ifstream& theFile,vector<int>& elemId,vector<int*>& elemTopology);
  void setCAPAVertexAndElement(vector<int>& vertexId,vector<double*>& vertexValues,
                               vector<int>& elemId,vector<int*>& elemTopology);

  // Private read functions
  void setVertexAndElement(vector<int*>& elem,vector<double*>& vtr);
};

} // namespace grd

#endif // MESHREADER_H
