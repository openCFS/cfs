/***************************************************************************
    File        : GridLevel.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/


#include "GridLevel.h"
#include "Element.h"
#include "Triangle.h"
#include "Quadrangle.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "Hexahedron.h"

namespace grd {


// Description
// Constructor
GridLevel::GridLevel()
{
  // Boundary 2D elements
  boundary[0] = &triangle;
  boundary[1] = &quadrangle;
  boundary[2] = 0;
  // Volume 3D elements
  volume[0] = &tetrahedron;
  volume[1] = &octahedron;
  volume[2] = &hexahedron;
  volume[3] = &prism;
  volume[4] = &pyramid;
  volume[5] = 0;

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


// Description
// Refine or coarse elements of grid level
inline void
refineLevelElement(std::list<Element*>& elemList,GridLevel& gridLevel)
{
  Element* element;
  std::list<Element*>::iterator p;

  // refine
  for (p = elemList.begin(); p != elemList.end(); ++p)
  {
    element = *p;
    if (element->isMarkedForRefinement()) {
      element->refine(gridLevel);
    }
    else if (element->isIrregular()) {
      if (element->queryChildrenEdgeRefinement())
        element->refine(gridLevel);
    }
    else if (element->isRefined()) {
      if (element->queryChildrenMarks())
        //element->coarsen(gridLevel);
        element->setRegular();
    }
  }
}

// Description
//  Close element of grid level
inline void
removeLevelElement(std::list<Element*>& lt)
{
  std::list<Element*>::iterator p = lt.begin();
  while (p != lt.end())
  {
    std::list<Element*>::iterator q = p;
    ++p;
    Element* element = *q;
    if (element->isMarkedForCoarsening())
    {
      delete *q;
      lt.erase(q);
    }
  }
}


// Description
//  Close element of grid level
inline void
closeLevelElement(std::list<Element*>& lt)
{
  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p)
  {
    Element* element = *p;
    if (element->isRegular())
    {
      if (element->queryEdgeRefinement())
      {
        element->close();
      }
    }
    else if (element->isIrregular())
    {
      if (!(element->queryEdgeRefinement()))
        element->setRegular();
    }
  }
}



// Description
//  Close element of grid level
inline void
closeLevelElement(std::list<Element*>& lt,GridLevel& gridLevel)
{
  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p)
  {
    Element* element = *p;
    if (element->isRegular())
    {
      if (element->queryEdgeRefinement())
      {
        ((Hexahedron*)element)->close(gridLevel);
      }
    }
    else if (element->isIrregular())
    {
      if (!(element->queryEdgeRefinement()))
        element->setRegular();
    }
  }
}


// Description
// refine boundary elements
// due to changes in the volume
inline void
refineLevelBoundary(std::list<Element*>& lt,GridLevel& gridLevel)
{
  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p)
  {
    Element* element  = *p;
    Element* neighbor = element->getNeighbor(0);
    if (element->isMarkedForRefinement())
    {
      if (neighbor) {
        if (!neighbor->isRefined())
          neighbor->markForRefinement();
      }
    }
    else if (element->isRegular()) {
      if (neighbor) {
        if (neighbor->isMarkedForRefinement()) {
          element->refine(gridLevel);
        }
      }
    }
    else if (element->isIrregular()) {
      if (element->queryChildrenEdgeRefinement())
        element->refine(gridLevel);
    }
  }
}



// Description
//  update boundary elements
//  due to changes in the volume
inline void
updateLevelBoundary(std::list<Element*>& lt,GridLevel& gridLevel)
{
  Element* element;
  Element* neighbor;

  typedef std::list<Element*>::iterator ElmI;
  for (ElmI p = lt.begin(); p != lt.end(); ++p) {
    element = *p;
    if (element->isRegular()) {
      neighbor = element->getNeighbor(0);
      if (neighbor) {
        if (neighbor->isRefined())
          element->refine(gridLevel);
      }
    }
    else if (element->isIrregular()) {
      neighbor = element->getNeighbor(0);
      if (neighbor) {
        if (neighbor->isRefined())
          element->refine(gridLevel);
      }
    }
    else if (element->isRefined()) {
      neighbor = element->getNeighbor(0);
      if (neighbor) {
        if (neighbor->isRegular())
        {
          std::cerr << "Oops!\n";
          element->coarsen(gridLevel);
        }
      }
    }
  }
}


//===============================================================

//===============================================================

// Description
// boundary refinement, elements already marked
void
GridLevel::refineBoundary(GridLevel& gridLevel)
{
  // Process each element list
  refineLevelElement(triangle,gridLevel);
  refineLevelElement(quadrangle,gridLevel);
}


// Description
// remove marked boundary elements
void
GridLevel::coarsenBoundary( )
{
  // Process each element list
  removeLevelElement(triangle);
  removeLevelElement(quadrangle);
}

// Description
// boundary refinement
void
GridLevel::closeBoundary()
{
  // Process each element list
  closeLevelElement(triangle);
  closeLevelElement(quadrangle);
}




void checkNeighbor(GridLevel& gridLevel)
{
  int counter = 0;
  std::list<Element*>* lt;
  typedef std::list<Element*>::iterator ElmI;

  lt = gridLevel.getTetrahedronList();
  for (ElmI p = lt->begin(); p != lt->end(); ++p) {
    Element* t = *p;
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
  for (ElmI p = lt->begin(); p != lt->end(); ++p) {
    Element* t = *p;
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

  cerr << "counter: " << counter << '\n';
}

// Description
// refine marked elements.
// if parent irregular mark parent for refinement
// insert new elements, vertices, etc. in next level
// Input:
//   nextLevel
//           pointer to next level grid to insert new elements
void
GridLevel::refine(GridLevel& gridLevel)
{
  // process first boundary elements
  refineLevelBoundary(triangle,gridLevel);
  refineLevelBoundary(quadrangle,gridLevel);

  // refine volume elements
  refineLevelElement(triangle,gridLevel);
  refineLevelElement(quadrangle,gridLevel);
  refineLevelElement(tetrahedron,gridLevel);
  refineLevelElement(octahedron,gridLevel);
  refineLevelElement(hexahedron,gridLevel);

  checkNeighbor(gridLevel);

  // update boundary
  // refine elements due to changes in the volume
  updateLevelBoundary(triangle,gridLevel);
  updateLevelBoundary(quadrangle,gridLevel);
}


// Description
// virtual closure
void
GridLevel::coarsen()
{
  removeLevelElement(triangle);
  removeLevelElement(quadrangle);
  removeLevelElement(tetrahedron);
  removeLevelElement(octahedron);
  removeLevelElement(hexahedron);
  removeLevelElement(prism);
  removeLevelElement(pyramid);
}



// Description
// virtual closure, care of new nodes
// at the barycenter
void
GridLevel::close(GridLevel& gridLevel)
{
  closeLevelElement(triangle);
  closeLevelElement(quadrangle);
  closeLevelElement(tetrahedron);
  closeLevelElement(octahedron);
  closeLevelElement(hexahedron,gridLevel);
  closeLevelElement(prism);
  closeLevelElement(pyramid);
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

  std::list<Element*>::iterator p;
  std::list<Element*>** lt = getElementList();
  // mark unreferenced edges and vertices
  for (int n = 0; n < 7; n++)
  {
    for (p = lt[n]->begin(); p != lt[n]->end(); ++p) {
      theElement = *p;
      noEdges = theElement->getNoOfEdges();
      for (i = 0; i < noEdges; i++) {
        theEdge = theElement->getEdge(i);
        theEdge->ref();
      }
      noVertices = theElement->getNoOfVertices();
      for (i= 0; i < noVertices; i++) {
        theVertex = theElement->getVertex(i);
        theVertex->ref();
      }
      int type = theElement->type();
      if (type == GRD_HEXAHEDRON)
      {
        theVertex = ((Hexahedron*)theElement)->getBarycenter();
        if (theVertex)
          theVertex->ref();
      }
    }
  }
}



// Description
// remove unused vertices and edges
void
GridLevel::removeUnrefVerticesAndEdges()
{
  Vertex* theVertex;
  Edge*   theEdge;
  Edge*   child0;
  Edge*   child1;

  std::list<Edge*>::iterator edg = edge.begin();
  while (edg != edge.end())
  {
    std::list<Edge*>::iterator p = edg;
    ++edg;
    theEdge = *p;
    if (!theEdge->queryRef()) {
      edge.erase(p);
      //delete theEdge;
    }
    else {
      //theEdge->ref();
      if (theEdge->isRefined()) {
        theEdge->getChildren(&child0,&child1);
        if (!child0->queryRef())
          theEdge->resetDescendant();
      }
    }

  }

  std::list<Vertex*>::iterator vt = vertex.begin();
  while (vt != vertex.end())
  {
    std::list<Vertex*>::iterator p = vt;
    ++vt;
    theVertex = *p;
    if (!theVertex->queryRef()) {
      vertex.erase(p);
      //delete theVertex;
    }
    else
      theVertex->ref();
  }
}




// Description
// unset the reference tag
void
GridLevel::unrefVerticesAndEdges()
{
  std::list<Edge*>::iterator edg;
  for (edg = edge.begin(); edg != edge.end(); ++edg)
    (*edg)->unref();
  std::list<Vertex*>::iterator vt;
  for (vt = vertex.begin(); vt != vertex.end(); ++vt)
    (*vt)->unref();
}



// Descritpion
// set global vertex numbers
void
GridLevel::setVertexNumbers(int& vertexId)
{
  std::list<Vertex*>::iterator p;
  for (p = vertex.begin(); p != vertex.end(); ++p) {
    (*p)->setId(vertexId);
    vertexId++;
  }
}



//
//
//
inline void
assembleTriangle(Triangle* tri,Vertex* vt[3],Edge* eg[3])
{
  int i;
  Vertex* vv1;
  Vertex* vv2;

  // set edges
  for (int edge = 0; edge < 3; edge++)
  {
    tri->setEdge(edge,eg[edge]);
    vv1 = eg[edge]->getVertex(0);
    vv2 = eg[edge]->getVertex(1);
    for (i = 0; i < 3; i++)
    {
      if (vt[i] != vv1 && vt[i] != vv2)
        tri->setVertex(edge,vt[i]);
    }
  }

  // edge 0
//  tri->setEdge(0,eg[0]);
//  vv1 = eg[0]->getVertex(0);
//  vv2 = eg[0]->getVertex(1);
//  for (i = 0; i < 3; i++)
//  {
//    if (vt[i] != vv1 && vt[i] != vv2)
//      tri->setVertex(0,vt[i]);
//  }

  // edge 1
//  tri->setEdge(1,eg[1]);
//  vv1 = eg[1]->getVertex(0);
//  vv2 = eg[1]->getVertex(1);
//  for (i = 0; i < 3; i++)
//  {
//    if (vt[i] != vv1 && vt[i] != vv2)
//      tri->setVertex(1,vt[i]);
//  }
//
//  // edge 2
//  tri->setEdge(2,eg[2]);
//  vv1 = eg[2]->getVertex(0);
//  vv2 = eg[2]->getVertex(1);
//  for (i = 0; i < 3; i++)
//  {
//    if (vt[i] != vv1 && vt[i] != vv2)
//      tri->setVertex(2,vt[i]);
//  }
}



// Description
//  compute the boundary from the volume
//  elements
void
GridLevel::extractBoundary()
{
  int i;
  int face;

  Vertex*  theVertex;
  Edge*    theEdge;
  Element* theElement;
  Element* theNeighbor;

  typedef std::list<Element*>::iterator ElmI;
  std::list<Element*>** lt = getVolumeElementList();
  int type = 0;
  while (lt[type])  // for each element type
  {
    for (ElmI p = lt[type]->begin(); p != lt[type]->end(); ++p)
    {
      theElement = *p;
      for (face = 0; face < theElement->getNoOfFaces(); face++)
      {
        theNeighbor = theElement->getNeighbor(face);
        if (!theNeighbor)
        {
          int noOfVertex = theElement->noOfVertexAtFace(face);
          // The neighbor is a triangle
          if ( noOfVertex == 3)
          {
            Triangle* tri = new Triangle;
            Vertex* vt[3];
            Edge*   eg[3];
            for (i = 0; i < 3; i++) {
              vt[i] = theElement->getVertexAtFace(face,i);
              eg[i] = theElement->getEdgeAtFace(face,i);
            }
            assembleTriangle(tri,vt,eg);
            // check triangle
            theEdge = tri->getEdge(0);
            Vertex* vv1 = tri->getVertex(1);
            Vertex* vv2 = tri->getVertex(2);
            theVertex = theEdge->getVertex(0);
            if (theVertex != vv1 && theVertex != vv2)
              std::cerr << "Panic!\n";
            triangle.push_back(tri);
            tri->setNeighbor(0,theElement);
            theElement->setNeighbor(face,tri);
          }
          // The neighbor is a quadrangle
          if ( noOfVertex == 4)
          {
            Quadrangle* qua = new Quadrangle;

            for (i = 0; i < 4; i++) {
              theVertex = theElement->getVertexAtFace(face,i);
              theEdge   = theElement->getEdgeAtFace(face,i);
              qua->setVertex(i,theVertex);
              qua->setEdge(i,theEdge);
            }

            quadrangle.push_back(qua);
            qua->setNeighbor(0,theElement);
            theElement->setNeighbor(face,qua);
          }
        } // the neighbor is null
      } // for each element's face
    } // for each element

    // get next element type
    type++;
  } // while element type

}


} // namespace grd
