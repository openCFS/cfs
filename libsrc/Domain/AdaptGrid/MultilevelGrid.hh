// - C++ -
/***************************************************************************
    File        : MultileveGrid.h
    Description :

 ---------------------------------------------------------------------------
    Begin       : Wed Dec 12 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/
#ifndef MULTILEVELGRID_H
#define MULTILEVELGRID_H


// System
#include <vector>


// Grid
#include "Stack.hh"
#include "Vertex.hh"
#include "Edge.hh"
#include "Element.hh"
#include "Triangle.hh"
#include "Quadrangle.hh"
#include "Tetrahedron.hh"
#include "Octahedron.hh"
#include "Hexahedron.hh"
//#include "Prism.hh"
#include "Pyramid.hh"
#include "GridLevel.hh"


namespace grd {


struct __VolumeIterator {
  typedef std::list<Element*>::iterator iterator;
  inline Element* operator* () { return *p; }
  inline __VolumeIterator operator= (const __VolumeIterator& vi)
  {
    level = vi.level;
    type  = vi.type;
    p     = vi.p;
    for (int nn = 0; nn < (int)vi.grid.size(); nn++)
      grid.push(vi.grid[nn]);

    return *this;
  }
  inline bool operator== (const __VolumeIterator& vi) { return (p == vi.p); }
  inline bool operator!= (const __VolumeIterator& vi) { return (p != vi.p); }
  
  inline void operator++ ()
  {
    if (p != lt[type]->end())
      ++p;
    for ( ; p != lt[type]->end(); ++p)
    {
      if ((*p)->isRegular() || (*p)->isClosure())
      {
        return;
      }
    }
    type++;
    while (lt[type])
    {
      for (p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular() || (*p)->isClosure())
        {
          return;
        }
      }
      type++;
    }
    level++;
    for ( ; level < (unsigned int)grid.size(); level++)
    {
      type = 2;
      lt = grid[level]->getElementList();
      while (lt[type])
      {
        for (p = lt[type]->begin(); p != lt[type]->end(); ++p)
        {
          if ((*p)->isRegular() || (*p)->isClosure())
          {
            return;
          }
        }
        type++;
      }
    }
    p = lt[type-1]->end();
    return;
  }
  

  // data
  unsigned int level;
  unsigned int type;
  iterator p;
  std::list<Element*>** lt;
  Stack<GridLevel*> grid;
};

class MultilevelGrid {
public:
  typedef __VolumeIterator VolumeIterator;
public:
  // Constructors
  MultilevelGrid();

  // Destructor
  virtual ~MultilevelGrid();

  // Mesh refinement
  void refine();
  void uniformRefine();

  // No closure
  void setGridWithClosure()    { closure = true;  }
  void setGridWithoutClosure() { closure = false; }
  bool isGridWithClosure()     { return closure;  }
  // Boundary elements
  void setGridWithBoundary()    { boundary = true;  }
  void setGridWithoutBoundary() { boundary = false; }
  bool isGridWithBoundary()     { return boundary;  }
  
  // Mesh sizes
  inline int getNoOfLevels();
  inline int getTopLevel();

  inline int getNoOfVertices();
  inline int getNoOfVertices(int lvl);

  inline int getNoOfEdges();
  inline int getNoOfEdges(int lvl);

  inline int getNoOfTriangles();
  inline int getNoOfTriangles(int lvl);
  inline int getNoOfQuadrangles();
  inline int getNoOfQuadrangles(int lvl);
  inline int getNoOfTetrahedra();
  inline int getNoOfTetrahedra(int lvl);
  inline int getNoOfOctahedra();
  inline int getNoOfOctahedra(int lvl);
  inline int getNoOfHexahedra();
  inline int getNoOfHexahedra(int lvl);
  inline int getNoOfPrisms();
  inline int getNoOfPrisms(int lvl);
  inline int getNoOfPyramids();
  inline int getNoOfPyramids(int lvl);

  inline int getNoOfElements();

  // set global vertex numbers
  inline void setVertexNumbers();

  // Generate the corase mesh from a list of
  // Vertices and a list of elements
  void buildCoarseGrid(std::vector<Vertex*>& vertex,std::vector<Element*>& element);

  // Topology
  inline bool isVolume();

  // Access operations
  GridLevel* getCoarseGrid()       { return grid[0];   }
  GridLevel* getGridLevel(int lvl) { return grid[lvl]; }

  std::list<Element*>*  getTriangleList(int lvl)    { return (grid[lvl]->getTriangleList());    }
  std::list<Element*>*  getQuadrangleList(int lvl)  { return (grid[lvl]->getQuadrangleList());  }
  std::list<Element*>*  getTetrahedronList(int lvl) { return (grid[lvl]->getTetrahedronList()); }
  std::list<Element*>*  getOctahedronList(int lvl)  { return (grid[lvl]->getOctahedronList());  }
  std::list<Element*>*  getHexahedronList(int lvl)  { return (grid[lvl]->getHexahedronList());  }
  std::list<Element*>*  getPrismList(int lvl)       { return (grid[lvl]->getPrismList());       }
  std::list<Element*>*  getPyramidList(int lvl)     { return (grid[lvl]->getPyramidList());     }
  std::list<Element*>** getElementList(int lvl)     { return (grid[lvl]->getElementList());     }
  void getVolumeGrid(std::vector<Element*>& grd);
  
  // Iteration
//   template<class T> inline void forEachElement(T& f);
//   template<class T> inline void forEachSurfaceElement(T& f);
//   template<class T> inline void forEachVolumeElement(T& f);
//   template<class T> inline void forEachVertex(T& f);
//   template<class T> inline void forEachEdge(T& f);

  // Test methods
  bool checkGridConsistencyAndConnectivity();

public:
  inline VolumeIterator beginVolume()
  {
    __VolumeIterator p;
    for (int nn = 0; nn <= topLevel; nn++) p.grid.push(grid[nn]);
    p.lt = grid[0]->getElementList();
    p.level = 0;
    p.type = 2; // tetrahedron
    p.p = p.lt[p.type]->begin();
    if (p.p == p.lt[p.type]->end()) ++p;

    return p;
  }
  inline VolumeIterator endVolume()
  {
    __VolumeIterator p;
    p.p = grid[topLevel]->getPyramidList()->end();
    return p;
  }
private:
  // Marks and tags
  bool closure;
  bool boundary;
  // no. of levels
  int topLevel;

  // no. of elements
  Stack<int> noOfVertices;
  Stack<int> noOfEdges;
  Stack<int> noOfTriangles;
  Stack<int> noOfQuadrangles;
  Stack<int> noOfTetrahedra;
  Stack<int> noOfOctahedra;
  Stack<int> noOfHexahedra;
  Stack<int> noOfPrisms;
  Stack<int> noOfPyramids;

  // multilevel grid
  Stack<GridLevel*> grid;

  // Private Methods
  // recompute the number of elements in the mesh
  void updateNoOfElements();
  // Extract the boundary mesh
  void extractBoundary();
  // Control Consistency of Element Refinement Marks
  void evaluateRefinementMarks();
  // Remove conforming closure
  void prepareGridForRefinement();
  // Compute the conforming clousre
  void completeGridWithClosureElements();
  // control function
  void checkEdgesOfAllElements(char* text);

};



// ____________________________________________________________
//
//
//           inline functions
//
// ____________________________________________________________



// Description
//
inline int
MultilevelGrid::getNoOfLevels()
{
  return (topLevel+1);
}


// Description
//
inline int
MultilevelGrid::getTopLevel()
{
  return topLevel;
}


// Description
//
inline int
MultilevelGrid::getNoOfVertices()
{
  return noOfVertices[topLevel];
}


// Description
//
inline int
MultilevelGrid::getNoOfVertices(int lvl)
{
  return noOfVertices[lvl];
}


// Description
//
inline int
MultilevelGrid::getNoOfEdges()
{
  return noOfEdges[topLevel];
}


// Description
//
inline int
MultilevelGrid::getNoOfEdges(int lvl)
{
  return noOfEdges[lvl];
}


// Description
//
inline int
MultilevelGrid::getNoOfTriangles()
{
  return noOfTriangles[topLevel];
}


// Description
//
inline int
MultilevelGrid::getNoOfTriangles(int lvl)
{
  return noOfTriangles[lvl];
}

//
//
inline int
MultilevelGrid::getNoOfQuadrangles()
{
  return noOfQuadrangles[topLevel];
}

//
//
inline int
MultilevelGrid::getNoOfQuadrangles(int lvl)
{
  return noOfQuadrangles[lvl];
}

// Description
//
inline int
MultilevelGrid::getNoOfTetrahedra()
{
    return noOfTetrahedra[topLevel];
}


// Description
//
inline int
MultilevelGrid::getNoOfTetrahedra(int lvl)
{
    return noOfTetrahedra[lvl];
}


// Description
//
inline int
MultilevelGrid::getNoOfOctahedra()
{
  return noOfOctahedra[topLevel];
}


// Description
//
inline int
MultilevelGrid::getNoOfOctahedra(int lvl)
{
  return noOfOctahedra[lvl];
}


// Description
//
inline int
MultilevelGrid::getNoOfHexahedra()
{
  return noOfHexahedra[topLevel];
}

// Description
//
inline int
MultilevelGrid::getNoOfHexahedra(int lvl)
{
  return noOfHexahedra[topLevel];
}

//
//
inline int
MultilevelGrid::getNoOfElements()
{
  int res = 0;
  res += noOfTriangles[topLevel];
  res += noOfQuadrangles[topLevel];
  res += noOfTetrahedra[topLevel];
  res += noOfOctahedra[topLevel];
  res += noOfHexahedra[topLevel];

  return res;
}

// set vertex numbers beginning with T_0 and continuing
// with T_1 \ T_0, T_2 \ T_1, ...
inline void
MultilevelGrid::setVertexNumbers()
{
  int vertexId = 0;

  for (int i = 0; i <= topLevel; i++)
    grid[i]->setVertexNumbers(vertexId);
}

// Description
// Topology
inline bool
MultilevelGrid::isVolume()
{
  if (grid[0]->noOfTetrahedra())
    return true;
  else if (grid[0]->noOfOctahedra())
    return true;
  else if (grid[0]->noOfHexahedra())
    return true;
  else if (grid[0]->noOfPrisms())
    return true;
  else if (grid[0]->noOfPyramids())
    return true;
  else
    return false;
}



// ____________________________________________________________
//
//  Iteration
//  Use function objects
// ____________________________________________________________

// // Iteration over Elements
// template<class T> inline void
// MultilevelGrid::forEachElement(T& f)
// {
//   typedef std::list<Element*>::iterator ElmI;
//   for (int j = 0; j <= topLevel; j++)
//   {
//     std::list<Element*>** lt = grid[j]->getElementList();
//     int type = 0;
//     while (lt[type])
//     {
//       if (type != 3)
//       {  // excludes octahedra from iteration
//         for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
//         {
//           if (!(*p)->isRefined())
//           {
//             f(*p);
//           }
//         }
//       } // if (type !=3)
//       // next element type
//       type++;
//     }
//   }
// }
// 
// 
// // Iteration over Vertices
// template<class T> inline void
// MultilevelGrid::forEachVertex(T& f)
// {
//   typedef std::list<Vertex*>::iterator VrtI;
//   for (int j = 0; j <= topLevel; j++) {
//     std::list<Vertex*>* lt = grid[j]->getVertexList();
//     for (VrtI p = lt->begin(); p != lt->end(); ++p)
//     {
//       f(*p);
//     }
//   }
// }
// 
// // Iteration over Edges
// template<class T> inline void
// MultilevelGrid::forEachEdge(T& f)
// {
//   typedef std::list<Edge*>::iterator EdgI;
//   for (int j = 0; j <= topLevel; j++) {
//     std::list<Edge*>* lt = grid[j]->getEdgeList();
//     for (EdgI p = lt->begin(); p != lt->end(); ++p)
//     {
//       if (!(*p)->isRefined())
//         f(*p);
//     }
//   }
// }
// 
// 


} // namespace grd

#endif // MULTILEVELGRID_H
