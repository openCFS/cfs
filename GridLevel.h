// - C++ -
/***************************************************************************
    File        : GridLevel.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef GRIDLEVEL_H
#define GRIDLEVEL_H


// System
#include <iostream>
#include <vector>
#include <list>

// Grid
//#include "Stack.h"
#include "Vertex.h"
#include "Edge.h"

namespace grd {

class Element;


class GridLevel {
public:
  // Constructors
  GridLevel();

  // Destructor
  virtual ~GridLevel();

  // Member Functions
  // Mesh refinement
  void refineBoundary(GridLevel& nextLevel);
  void coarsenBoundary();
  void closeBoundary();
  void refine(GridLevel& nextLevel);
  void coarsen();
  void close(GridLevel& nextLevel);

  void refVerticesAndEdges();
  void removeUnrefVerticesAndEdges();
  void unrefVerticesAndEdges();
  void resetDescendant();

  // Mesh sizes
  inline bool empty();
  inline int  noOfVertices();
  inline int  noOfEdges();
  inline int  noOfTriangles();
  inline int  noOfQuadrangles();
  inline int  noOfTetrahedra();
  inline int  noOfOctahedra();
  inline int  noOfHexahedra();
  inline int  noOfPrisms();
  inline int  noOfPyramids();
  inline int  noOfElements();
  inline int  noOfBoundaryElements();
  inline int  noOfVolumeElements();

  // set global vertex numbers
  void setVertexNumbers(int& vtId);

  // boundary
  void extractBoundary();

  // list operations
  inline void appendVertex(Vertex* vt);
  inline void appendEdge(Edge* edg);
  inline void appendTriangle(Element* tri);
  inline void appendQuadrangle(Element* quad);
  inline void appendTetrahedron(Element* tetra);
  inline void appendOctahedron(Element* octa);
  inline void appendHexahedron(Element* hexa);
  inline void appendPrism(Element* pri);
  inline void appendPyramid(Element* pyz);

  inline void removeVertex(Vertex* vt);
  inline void removeEdge(Edge* edg);
  inline void removeTriangle(Element* tri);
  inline void removeQuadrangle(Element* qua);
  inline void removeTetrahedron(Element* tetra);
  inline void removeOctahedron(Element* octa);
  inline void removeHexahedron(Element* hexa);
  inline void removePrism(Element* pris);
  inline void removePyramid(Element* pyra);

  inline std::list<Vertex*>*  getVertexList() { return &vertex; }
  inline std::list<Edge*>*    getEdgeList()   { return &edge;   }

  inline std::list<Element*>** getElementList()         { return allElement; }
  inline std::list<Element*>** getBoundaryElementList() { return boundary;   }
  inline std::list<Element*>** getVolumeElementList()   { return volume;     }

  std::list<Element*>* getTriangleList()    { return &triangle;    }
  std::list<Element*>* getQuadrangleList()  { return &quadrangle;  }
  std::list<Element*>* getTetrahedronList() { return &tetrahedron; }
  std::list<Element*>* getOctahedronList()  { return &octahedron;  }
  std::list<Element*>* getHexahedronList()  { return &hexahedron;  }
  std::list<Element*>* getPrismList()       { return &prism;       }
  std::list<Element*>* getPyramidList()     { return &pyramid;     }

  // iteration
  template<class T> inline void forEachElement(T& f);

private:
  std::list<Vertex*>  vertex;
  std::list<Edge*>    edge;
  std::list<Element*> triangle;
  std::list<Element*> quadrangle;
  std::list<Element*> tetrahedron;
  std::list<Element*> octahedron;
  std::list<Element*> hexahedron;
  std::list<Element*> prism;
  std::list<Element*> pyramid;

  // Lists for element iteratio
  std::list<Element*>* boundary[3];
  std::list<Element*>* volume[6];
  std::list<Element*>* allElement[8];

};


// ____________________________________________________________
//
//	
//	         inline functions
//
// ____________________________________________________________
//


// Description
//
inline bool
GridLevel::empty()
{
  return (vertex.empty());
}


// Desctiption
//
inline int
GridLevel::noOfVertices()
{
  return (vertex.size());
}


// Description
//
inline int
GridLevel::noOfEdges()
{
  return (edge.size());
}


// Description
//
inline int
GridLevel::noOfTriangles()
{
  return (triangle.size());
}


// Description
//
inline int
GridLevel::noOfQuadrangles()
{
  return (quadrangle.size());
}


// Description
//
inline int
GridLevel::noOfTetrahedra()
{
  return (tetrahedron.size());
}


// Description
//
inline int
GridLevel::noOfOctahedra()
{
  return (octahedron.size());
}

// Description
//
inline int
GridLevel::noOfHexahedra()
{
  return (hexahedron.size());
}

// Description
//
inline int
GridLevel::noOfPrisms()
{
  return (prism.size());
}

// Description
//
inline int
GridLevel::noOfPyramids()
{
  return (pyramid.size());
}


// Description
//
inline int
GridLevel::noOfElements()
{
  int res = 0;
  res += triangle.size();
  res += quadrangle.size();
  res += tetrahedron.size();
  res += octahedron.size();
  res += hexahedron.size();
  res += prism.size();
  res += pyramid.size();

  return res;
}

// Description
//
inline int
GridLevel::noOfBoundaryElements()
{
  int res = 0;
  res += triangle.size();
  res += quadrangle.size();

  return res;
}


// Description
//
inline int
GridLevel::noOfVolumeElements()
{
  int res = 0;
  res += tetrahedron.size();
  res += octahedron.size();
  res += hexahedron.size();
  res += prism.size();
  res += pyramid.size();
  return res;
}

// Description
//
inline void
GridLevel::appendVertex(Vertex* vt)
{
  vertex.push_back(vt);
}


// Description
//
inline void
GridLevel::appendEdge(Edge* edg)
{
  edge.push_back(edg);
}

// Description
//
inline void
GridLevel::appendTriangle(Element* tri)
{
  triangle.push_back(tri);
}


// Description
//
inline void
GridLevel::appendQuadrangle(Element* quad)
{
  quadrangle.push_back(quad);
}


// Description
//
inline void
GridLevel::appendTetrahedron(Element* tetra)
{
  tetrahedron.push_back(tetra);
}

// Description
//
inline void
GridLevel::appendOctahedron(Element* octa)
{
  octahedron.push_back(octa);
}

// Description
//
inline void
GridLevel::appendHexahedron(Element* hexa)
{
  hexahedron.push_back(hexa);
}


// Description
//
inline void
GridLevel::appendPrism(Element* pri)
{
  prism.push_back(pri);
}

// Description
//
inline void
GridLevel::appendPyramid(Element* pyr)
{
  pyramid.push_back(pyr);
}

// Description
//
inline void
GridLevel::removeVertex(Vertex* vt)
{
  vertex.remove(vt);
}

// Description
//
inline void
GridLevel::removeEdge(Edge* edg)
{
  edge.remove(edg);
}

// Description
//
inline void
GridLevel::removeTriangle(Element* tri)
{
  triangle.remove(tri);
}

// Description
//
inline void
GridLevel::removeQuadrangle(Element* qua)
{
  quadrangle.remove(qua);
}

// Description
//
inline void
GridLevel::removeTetrahedron(Element* tetra)
{
  tetrahedron.remove(tetra);
}

// Description
//
inline void
GridLevel::removeOctahedron(Element* octa)
{
  octahedron.remove(octa);
}

//
//
// Description
//
inline void
GridLevel::removeHexahedron(Element* hexa)
{
  octahedron.remove(hexa);
}

//
//
// Description
//
inline void
GridLevel::removePrism(Element* pris)
{
  octahedron.remove(pris);
}

//
//
// Description
//
inline void
GridLevel::removePyramid(Element* pyra)
{
  octahedron.remove(pyra);
}


//
//
//
template<class T> inline void
GridLevel::forEachElement(T& f)
{
  typedef std::list<Element*>::iterator ElmI;

  std::list<Element*>** lt = getElementList();
  int type = 0;
  while (lt[type]) {
    for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
    {
      // call functor
      f(*p);
    }
    // next element type
    type++;
  }

}

} // namespace grd


#endif // GRIDLEVEL_H

