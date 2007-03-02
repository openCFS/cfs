/***************************************************************************
    File        : GridLevel.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/


#include "GridLevel.hh"
#include "Element.hh"
#include "Triangle.hh"
#include "Quadrangle.hh"
#include "Tetrahedron.hh"
#include "Octahedron.hh"
#include "Hexahedron.hh"

// Utilities
#include "Utility.hh"

namespace grd {


// Description
// Constructor
GridLevel::GridLevel()
{
  // All elements
  allElement[0] = &triangle;
  allElement[1] = &quadrangle;
  allElement[2] = &tetrahedron;
  allElement[3] = &octahedron;
  allElement[4] = &hexahedron;
  allElement[5] = &prism;
  allElement[6] = &pyramid;
  allElement[7] = 0;
}

// Description:
// delete all elements in the list
GridLevel::~GridLevel()
{
  // remove vertices
  typedef std::list<Vertex*>::iterator vtLI;
  for (vtLI vt = vertex.begin(); vt != vertex.end(); ++vt)
    delete *vt;
  vertex.clear();
  // remove edges
  typedef std::list<Edge*>::iterator edgLI;
  for (edgLI edg = edge.begin(); edg != edge.end(); ++edg)
    delete *edg;
  edge.clear();

  // delete all elements and clear lists
  typedef std::list<Element*>::iterator ElmLI;
  std::list<Element*>** lt = getElementList();
  int type = 0;
  while (lt[type])
  {
    for (ElmLI p = lt[type]->begin(); p != lt[type]->end(); ++p)
      delete *p;
    lt[type]->clear();
    // next element type
    type++;
  }
}



//===============================================================
//
//   Inline Functions for Mesh Refinement
//
//===============================================================
inline void refineLevelElement(std::list<Element*>& lt,GridLevel& gridLevel);
inline void coarseLevelElement(std::list<Element*>& lt);
inline void closeLevelElement(std::list<Element*>& lt,GridLevel& gridLevel);
inline void markLevelBoundary(std::list<Element*>& lt);
//===============================================================

//===============================================================
// Description
// refine elements
// Input:
//   nextLevel: pointer to next level grid to insert new elements
void
GridLevel::refine(GridLevel& gridLevel)
{
  // Refine volume element first
  refineLevelElement(tetrahedron,gridLevel);
  refineLevelElement(octahedron,gridLevel);
  refineLevelElement(hexahedron,gridLevel);
  refineLevelElement(prism,gridLevel);
  refineLevelElement(pyramid,gridLevel);
  // Refine boundary elements
  markLevelBoundary(triangle);
  markLevelBoundary(quadrangle);
  refineLevelElement(triangle,gridLevel);
  refineLevelElement(quadrangle,gridLevel);

} // refine()


// Description
// virtual closure
void
GridLevel::coarsen()
{
  // Coarse volume elements
  coarseLevelElement(tetrahedron);
  coarseLevelElement(octahedron);
  coarseLevelElement(hexahedron);
  coarseLevelElement(prism);
  coarseLevelElement(pyramid);
  // Coarse boundary elements
  coarseLevelElement(triangle);
  coarseLevelElement(quadrangle);

}



// Description
// virtual closure, care of new nodes
// at the barycenter
void
GridLevel::close(GridLevel& gridLevel)
{
  closeLevelElement(triangle,gridLevel);
  closeLevelElement(quadrangle,gridLevel);
  closeLevelElement(tetrahedron,gridLevel);
  closeLevelElement(octahedron,gridLevel);
  closeLevelElement(hexahedron,gridLevel);
  closeLevelElement(prism,gridLevel);
  closeLevelElement(pyramid,gridLevel);
}

// Description
//  split octahedra into four tetrahedra
void
GridLevel::splitOctahedra()
{
  std::list<Element*>::iterator p;
  // Split octas
  for (p = octahedron.begin(); p != octahedron.end(); ++p)
  {
    if ((*p)->isRegular())
      ((Octahedron*)(*p))->split(*this);
  }
}


// Description
//  unsplit octahedra
//  remove tetrahedra from mesh
void
GridLevel::unsplitOctahedra()
{
  std::list<Element*>::iterator p;
  
  // Unplit octahedra
  // mark unused tetrahedra and edges
  for (p = octahedron.begin(); p != octahedron.end(); ++p)
  {
    if ((*p)->isRegular())
      ((Octahedron*)(*p))->unsplit();   
  }
}

// Description
// remove the unused elements for the octahedron split
void 
GridLevel::removeOctahedronSplitElements()
{
  // remove unused tetrahedra and edges
  std::list<Element*>::iterator p = tetrahedron.begin();
  while (p != tetrahedron.end())
  {
    if ((*p)->isRegular() && (*p)->queryRef())
    {
      delete *p;  
      p = tetrahedron.erase(p);
    }
    else
      ++p;
  }

  std::list<Edge*>::iterator q = edge.begin();
  while (q != edge.end())
  {
    if ((*q)->queryRef())
    {
      delete *q;
      q = edge.erase(q);
    }
    else
      ++q;
  }
}


// Description
// reference vertices and edges
void
GridLevel::refVerticesAndEdges()
{
  int i,noVertices,noEdges;
  Vertex*  theVertex;
  Edge*    theEdge;
  Element* theElement;

  typedef std::list<Element*>::iterator IT;
  std::list<Element*>** lt = getElementList();

  // mark unreferenced edges and vertices
  int et = 0;
  while(lt[et])
  {
    for (IT p = lt[et]->begin(); p != lt[et]->end(); ++p)
    {
      theElement = *p;
      // process edges
      noEdges = theElement->getNoOfEdges();
      for (i = 0; i < noEdges; i++)
      {
        theEdge = theElement->getEdge(i);
        if (!theEdge)
          std::cout << "ERROR: element has NULL edge" << std::endl;
        theEdge->ref();
      }
      // process vertices
      noVertices = theElement->getNoOfVertices();
      for (i= 0; i < noVertices; i++)
      {
        theVertex = theElement->getVertex(i);
        theVertex->ref();
      }
    }
    // get next element type
    et++;
  }
}



// Description
// remove unused vertices and edges
void
GridLevel::removeUnrefVerticesAndEdges()
{
  // remove all unused edges
  std::list<Edge*>::iterator p = edge.begin();
  while (p != edge.end())
  {
    // if edge is not referenced
    // remove from list
    if (!(*p)->queryRef())
    {
      delete (*p);
      p = edge.erase(p);
    }
    else
    {
      // if children will be removed from list
      // reset the parent child hierarchy here
      if ((*p)->isRefined())
      {
        Edge* child0;
        Edge* child1;
        (*p)->getChildren(&child0,&child1);
        if (!child0->queryRef() && !child1->queryRef())
          (*p)->resetDescendant();
      }
      // update pointer
      ++p;
    }

  }

  // remove all unused vertices
  std::list<Vertex*>::iterator q = vertex.begin();;
  while (q != vertex.end())
  {
    // if vertex not referenced
    // remove from list
    if (!(*q)->queryRef())
    {
      std::cout << "Erase element from list" << std::endl;
      delete (*q);
      q = vertex.erase(q);
    }
    else
    {
      // update pointer
      ++q;
    }
  }
}



// Description
// unset the reference tag
void
GridLevel::unrefVerticesAndEdges()
{
  std::list<Edge*>::iterator p;
  for (p = edge.begin(); p != edge.end(); ++p)
  {
    (*p)->unref();
  }
  std::list<Vertex*>::iterator q;
  for (q = vertex.begin(); q != vertex.end(); ++q)
    (*q)->unref();
}



// Descritpion
// set global vertex numbers
void
GridLevel::setVertexNumbers(int& vertexId)
{
  std::list<Vertex*>::iterator p;
  for (p = vertex.begin(); p != vertex.end(); ++p)
  {
    (*p)->setId(vertexId);
    vertexId++;
  }
}



//===============================================================
//
//   Inline Functions for Mesh Refinement
//
//===============================================================


// Description
// Refine for grid level
//inline void
inline void
refineLevelElement(std::list<Element*>& elemList,GridLevel& gridLevel)
{
  std::list<Element*>::iterator p;

  // refine
  for (p = elemList.begin(); p != elemList.end(); ++p)
  {
    if ((*p)->isMarkedForRefinement()) {
      (*p)->refine(gridLevel);
    }
    else if ((*p)->isIrregular()) {
      if ((*p)->queryChildrenEdgeRefinement())
        (*p)->refine(gridLevel);
    }
    else if ((*p)->isRefined()) {
      if ((*p)->queryChildrenMarks())
        (*p)->setRegular();
    }
  }
}

// Description
//  Close element of grid level
inline void
coarseLevelElement(std::list<Element*>& lt)
{
  std::list<Element*>::iterator p = lt.begin();
  while (p != lt.end())
  {
    if ((*p)->isMarkedForCoarsening())
    {
      delete *p;
      p = lt.erase(p);
    }
    else
      ++p;
  }
}



// Description
//  Close element of grid level
//  Add children to next level grid
inline void
closeLevelElement(std::list<Element*>& lt,GridLevel& gridLevel)
{
  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p)
  {
    (*p)->unsetAllRefinementMarks();
    if ((*p)->isRegular())
    {
      if ((*p)->queryEdgeRefinement())
      {
        (*p)->close();
      }
    }
    else if ((*p)->isIrregular())
    {
      if (!((*p)->queryEdgeRefinement()))
        (*p)->setRegular();
    }
  }
}


inline void
markLevelBoundary(std::list<Element*>& lt)
{
  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p)
  {
    if ((*p)->isRefined())
      continue;
    Element* neighbor = (*p)->getNeighbor(0);
    if (neighbor && neighbor->isRefined())
      (*p)->markForRefinement();
  }
}



// Description
//  check neighbors
void checkNeighbor(GridLevel& gridLevel)
{
  int counter = 0;
  std::list<Element*>* lt;
  typedef std::list<Element*>::iterator ElmI;

  lt = gridLevel.getTetrahedronList();
  for (ElmI ei = lt->begin(); ei != lt->end(); ++ei) {
    Element* t = *ei;
    int noOfFaces = t->getNoOfFaces();
    grd::Element* p = t->getParent();
    for (int i = 0; i < noOfFaces; i++)
    {
      grd::Element* n = t->getNeighbor(i);
      //int type = n->type();
      //if (type == GRD_TRIANGLE)
      if (n == 0 && p != 0) counter++;
    }
  }

  lt = gridLevel.getOctahedronList();
  for (ElmI ei = lt->begin(); ei != lt->end(); ++ei) {
    Element* t = *ei;
    int noOfFaces = t->getNoOfFaces();
    grd::Element* p = t->getParent();
    for (int i = 0; i < noOfFaces; i++)
    {
      grd::Element* n = t->getNeighbor(i);
      //int type = n->type();
      //if (type == GRD_TRIANGLE)
      if (n == 0 && p != 0) counter++;
    }
  }
}



} // namespace grd
