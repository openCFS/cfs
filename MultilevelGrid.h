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
#include "Stack.h"
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "GridLevel.h"
#include "ConformingClosure.h"

// Surface elements
#include "Triangle.h"
#include "TriSkeleton.h"


namespace grd {


class MultilevelGrid {
public:
  // Constructors
  MultilevelGrid();

  // Destructor
  virtual ~MultilevelGrid();

  // Mesh refinement
  void refineBoundary();
  void refine();

  // Mesh sizes
  inline int getNoOfLevels();
  inline int getTopLevel();

  inline int getNoOfVertices();
  inline int getNoOfVertices(int lvl);

  inline int getNoOfEdges();
  inline int getNoOfEdges(int lvl);

  inline int getNoOfTriangles();
  inline int getNoOfTriangles(int lvl);
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
  void buildCoarseMesh(std::vector<Vertex*>& vertex,vector<Element*>& element);
  void setBoundaryMesh(std::vector<Element*>& element);
  void buildMultilevelBoundary(Stack<double**>& mlVertex,Stack<int**>& mlElement);

  // Initialize Grid
  GridLevel* getCoarseGrid() { return grid[0]; }
  // recompute the number of elements in the mesh
  void updateNoOfElements();

  // Grid levels access
  GridLevel* getGridLevel(int lvl) { return grid[lvl]; }

  // Iteration
  template<class T> inline void forEachElement(T& f);
  template<class T> inline void forEachBoundaryElement(T& f);
  template<class T> inline void forEachVolumeElement(T& f);

  template<class T> inline void forEachUnrefinedBoundaryElement(T& f);
  template<class T> inline void forEachUnrefinedVolumeElement( T& f);
  template<class T> inline void forEachRegAndIrregElement(T& f);
  template<class T> inline void forEachUnrefinedElement(T& f);

  template<class T> inline void forEachVertex(T& f);
  template<class T> inline void forEachEdge(T& f);
  template<class T> inline void forEachEdge(int lv, T& f);

private:
  // no. of levels
  int topLevel;
  int bndTopLevel;

  // element sets
  int noOfFaceSets;

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

  // Conforming closure
  ConformingClosure closure;
};



// ____________________________________________________________
//
//	
//	         inline functions
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


// ____________________________________________________________
//
//  Iteration
//  Use function objects
// ____________________________________________________________

// Iteration over Elements
template<class T> inline void
MultilevelGrid::forEachElement(T& f)
{
  typedef list<Element*>::iterator ElmI;
  typedef ConformingClosure::triangleIterator    TriI;
  typedef ConformingClosure::tetrahedronIterator TetI;
  for (int j = 0; j <= topLevel; j++) {
    list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type]) {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular()) {
          f(*p);
        }
        else if ((*p)->isIrregular()) {
          (*p)->close(closure);
          for (TriI tri = closure.beginTriangle(); tri != closure.endTriangle(); ++tri)
            f(*tri);
          for (TetI tet = closure.beginTetrahedron(); tet != closure.endTetrahedron(); ++tet)
            f(*tet);
        }
      }
      // next element type
      type++;
    }
  }
}


//
//
template<class T> inline void
MultilevelGrid::forEachRegAndIrregElement(T& f)
{
  typedef list<Element*>::iterator ElmI;
  typedef ConformingClosure::tetrahedronIterator TetI;
  for (int j = 0; j <= topLevel; j++)
  {
    list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular()) {
          f(*p);
        }
        else if ((*p)->isIrregular()) {
          f(*p);
        }
      }
      // next element type
      type++;
    } // while element type
  } // for each grid level
}


//
//
template<class T> inline void
MultilevelGrid::forEachUnrefinedElement(T& f)
{
  typedef list<Element*>::iterator ElmI;
  typedef ConformingClosure::tetrahedronIterator TetI;
  for (int j = 0; j <= topLevel; j++)
  {
    list<Element*>** lt = grid[j]->getElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if (!(*p)->isRefined()) {
          f(*p);
        }
      }
      // next element type
      type++;
    } // while element type
  } // for each grid level
}

//
//
//
template<class T> inline void
MultilevelGrid::forEachUnrefinedBoundaryElement( T& f)
{
  typedef list<Element*>::iterator ElmI;
  for (int j = 0; j <= topLevel; j++)
  {
    list<Element*>* lt = grid[j]->getTriangleList();
    for (ElmI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined()) {
        f(*p);
      }
    }

    lt = grid[j]->getQuadrangleList();
    for (ElmI q = lt->begin(); q != lt->end(); ++q)
    {
      if (!(*q)->isRefined()) {
        f(*q);
      }
    }
  } // for each grid level
}


//
//
//
template<class T> inline void
MultilevelGrid::forEachUnrefinedVolumeElement( T& f)
{
  typedef list<Element*>::iterator ElmI;
  for (int j = 0; j <= bndTopLevel; j++) {
    list<Element*>** lt = grid[j]->getVolumeElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular())
        {
          f(*p);
        }
        else if ((*p)->isIrregular())
        {
          f(*p);
        }
      }
      // get next element type
      type++;
    } // while element type
  } // for topLevele
}


//
//
//
template<class T> inline void
MultilevelGrid::forEachBoundaryElement(T& f)
{
  typedef list<Element*>::iterator ElmI;
  typedef ConformingClosure::triangleIterator TriI;
  for (int j = 0; j <= bndTopLevel; j++) {
    list<Element*>** lt = grid[j]->getBoundaryElementList();
    for (int n = 0; n < 2; n++)
    {
      for (ElmI p = lt[n]->begin(); p != lt[n]->end(); ++p)
      {
        if ((*p)->isRegular())
        {
          f(*p);
        }
        else if ((*p)->isIrregular())
        {
          (*p)->close(closure);
          for (TriI tri = closure.beginTriangle(); tri != closure.endTriangle(); ++tri)
            f(*tri);
        }
      }
    } // for element type
  } // for topLevele
}


//
//
//
template<class T> inline void
MultilevelGrid::forEachVolumeElement(T& f)
{
  typedef list<Element*>::iterator ElmI;
  typedef ConformingClosure::tetrahedronIterator TetI;
  for (int j = 0; j <= bndTopLevel; j++) {
    list<Element*>** lt = grid[j]->getVolumeElementList();
    int type = 0;
    while (lt[type])
    {
      for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      {
        if ((*p)->isRegular())
        {
          f(*p);
        }
        else if ((*p)->isIrregular())
        {
          (*p)->close(closure);
          for (TetI tet = closure.beginTetrahedron(); tet != closure.endTetrahedron(); ++tet)
            f(*tet);
        }
      }
      // get next element type
      type++;
    } // while element type
  } // for topLevele
}


// Iteration over Vertices
template<class T> inline void
MultilevelGrid::forEachVertex(T& f)
{
  typedef list<Vertex*>::iterator VrtI;
  for (int j = 0; j <= topLevel; j++) {
    list<Vertex*>* lt = grid[j]->getVertexList();
    for (VrtI p = lt->begin(); p != lt->end(); ++p)
    {
      f(*p);
    }
  }
}

// Iteration over Edges
template<class T> inline void
MultilevelGrid::forEachEdge(T& f)
{
  typedef list<Edge*>::iterator EdgI;
  for (int j = 0; j <= topLevel; j++) {
    list<Edge*>* lt = grid[j]->getEdgeList();
    for (EdgI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        f(*p);
    }
  }
}


// Iteration over Edges
template<class T> inline void
MultilevelGrid::forEachEdge(int lvl,T& f)
{
  typedef list<Edge*>::iterator EdgI;
  if (lvl > topLevel)
  {
    cerr << "ERROR: can't iterate edges up to level: " << lvl << '\n';
    exit(0);
  }

  for (int j = 0; j <= lvl; j++) {
    list<Edge*>* lt = grid[j]->getEdgeList();
    for (EdgI p = lt->begin(); p != lt->end(); ++p)
    {
      if (!(*p)->isRefined())
        f(*p);
    }
  }
}


} // namespace grd

#endif // MULTILEVELGRID_H
