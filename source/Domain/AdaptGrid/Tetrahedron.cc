/***************************************************************************
    File        : Tetrahedron.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Tetrahedron.hh"
#include "Octahedron.hh"
#include "GridLevel.hh"

namespace grd {

// (static obj) initialize memory management here
//Pool<Tetrahedron> Tetrahedron::pool;

// Description:
//
Tetrahedron::Tetrahedron( )
{
  int i;
  // Vertices and neighbors
  for (i = 0; i < 4; i++) {
    v[i] = 0;
    n[i] = 0;
  }
  // Edges
  for (i = 0; i < 6; i++)
    e[i] = 0;

  // set no. of children
  // chcount = 0; <- set in abstract class Element
}



//=============================================================================
//  REFINEMENT
//=============================================================================

// Description
//
// refine tetrahedron regular into 4 tetrahedra and one octahedra
// insert new elements in the lists.
// input: tesselation - contains pointer to the element lists
// output: bool true or false
// -* 1. process parent vertices and get new vertices
// -* 2. process parent edges and get new vertices + new edges
// -* 3. process parent faces and get new edges + neighbors
// *- 4. assemble children from new elements
void
Tetrahedron::refine(GridLevel& gridLevel)
{
  int ch;
  int vertCount,edgCount;
  int vertex,edge,face,neighbor;
  int vertexIndex,edgeIndex,neighIndex,childIndex;
  int nVertex,nFace;
  int neighFace[16],theChildFace[4];

  Vertex*  newVertex[11];
  Edge*    newEdge[24];
  Element* newNeighbor[21];
  Vertex*  theVertex[3];
  Edge*    theEdge[3];
  Element* theChild[4];

  // -* init counters
  vertCount = edgCount = 0;

  // -* 1. process parent vertices
  for (vertex = 0; vertex < 4; vertex++)
    newVertex[vertCount++] = v[vertex];

  // -* 2. process parent edges
  for (edge = 0; edge < 6; edge++) {
    if (e[edge]->isRefined()) {
      e[edge]->getChildren(newVertex[tetVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[tetVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      newVertex[vertCount] = e[edge]->getMidPoint();

      // pend new elements in list
      gridLevel.appendVertex(newVertex[vertCount]);
      gridLevel.appendEdge(newEdge[edgCount]);
      gridLevel.appendEdge(newEdge[edgCount+1]);

      edgCount += 2;
      vertCount++;
    }
  }

//  // -* get vertex at barycenter of parent element
//
//  Vertex* baryCenter = new Vertex;
//
//  (*baryCenter)[0] = 0.25*((*v[0])[0]+(*v[1])[0]+(*v[2])[0]+(*v[3])[0]);
//  (*baryCenter)[1] = 0.25*((*v[0])[1]+(*v[1])[1]+(*v[2])[1]+(*v[3])[1]);
//  (*baryCenter)[2] = 0.25*((*v[0])[2]+(*v[1])[2]+(*v[2])[2]+(*v[3])[2]);
//  gridLevel.appendVertex(baryCenter);
//
//  newVertex[10] = baryCenter;


  // -* 3. process parent faces

  for (face = 0; face < 4; face++) {
    // if not at boundary
    if (n[face]) {
      // if neighbor regular refined
      if (n[face]->isRefined()) {
//        int type = n[face]->type();
//        if (type == GRD_TRIANGLE)
//          tri = true;
        // get neighbor common face
        nFace = n[face]->getCommonFace(this);
        n[face]->getRefElements(nFace,theVertex,theEdge,theChild,theChildFace);

        // for all vertices of neighbor at neighbor's face
        int index = 2;
        vertex    = 2;
        for (nVertex = 0; nVertex < 3; nVertex++) {
          // find common vertex of parent and neighbor
          vertex      = ((index--)%3);
          vertexIndex = tetVertexOfFace[face][vertex];
          if (index < 0) index = 2;
          while(v[vertexIndex] != theVertex[nVertex]) {
            vertex      = ((index--)%3);
            vertexIndex = tetVertexOfFace[face][vertex];
            if (index < 0) index = 2;
          }

          // set child edge
          edgeIndex = tetEdgeFromParentFace[face][vertex];
          newEdge[edgeIndex] = theEdge[nVertex];
          // set child neighbor
          neighIndex = tetNeighborOfParentFace[face][vertex];
          newNeighbor[neighIndex] = theChild[nVertex];
          // remember which face of neighbor to connenct
          neighFace[neighIndex] = theChildFace[nVertex];
        }
        // set last neighbor from face middle
        neighIndex = tetNeighborOfParentFace[face][3];
        newNeighbor[neighIndex] = theChild[3];
        neighFace[neighIndex] = theChildFace[3];
      }
      else {
        // neighbor is not refined
        // create new edges on face
        // set children neighbors at face to NULL
        int v1,v2;
        for (edge = 0; edge < 3; edge++) {
          Edge* tmpEdge = new Edge;
          v1 = tetVertexOfEdgeOfParentFace[face][edge][0];
          v2 = tetVertexOfEdgeOfParentFace[face][edge][1];
          tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
          newEdge[tetEdgeFromParentFace[face][edge]] = tmpEdge;
          gridLevel.appendEdge(tmpEdge);
        }
        for (neighbor = 0; neighbor < 4; neighbor++) {
          neighIndex = tetNeighborOfParentFace[face][neighbor];
          newNeighbor[neighIndex] = 0;
        }
      }
    }
    else {
      // neighbor is not refined
      // create new edges on face
      // set children neighbors at face to NULL
      int v1,v2;
      for (edge = 0; edge < 3; edge++) {
        Edge* tmpEdge = new Edge;
        v1 = tetVertexOfEdgeOfParentFace[face][edge][0];
        v2 = tetVertexOfEdgeOfParentFace[face][edge][1];
        tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
        newEdge[tetEdgeFromParentFace[face][edge]] = tmpEdge;
        gridLevel.appendEdge(tmpEdge);
      }
      for (neighbor = 0; neighbor < 4; neighbor++) {
        neighIndex = tetNeighborOfParentFace[face][neighbor];
        newNeighbor[neighIndex] = 0;
      }
    }
  }


  // create children
  child = new Element*[5];

  child[0] = (Element*) new Tetrahedron;
  child[1] = (Element*) new Tetrahedron;
  child[2] = (Element*) new Tetrahedron;
  child[3] = (Element*) new Tetrahedron;
  child[4] = (Element*) new Octahedron;

  gridLevel.appendTetrahedron(child[0]);
  gridLevel.appendTetrahedron(child[1]);
  gridLevel.appendTetrahedron(child[2]);
  gridLevel.appendTetrahedron(child[3]);
  gridLevel.appendOctahedron(child[4]);

  // set newNeighbor for interior connections
  newNeighbor[16] = child[0];
  newNeighbor[17] = child[1];
  newNeighbor[18] = child[2];
  newNeighbor[19] = child[3];
  newNeighbor[20] = child[4];

  // set children parent
  child[0]->setParent(this);
  child[1]->setParent(this);
  child[2]->setParent(this);
  child[3]->setParent(this);
  child[4]->setParent(this);


  // assemble tetra children from new elements
  for (ch = 0; ch < 4; ch++) {
    // set vertices
    for (vertex = 0; vertex < 4; vertex++) {
      vertexIndex = tetVertexOfTetraChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 6; edge++) {
      edgeIndex = tetEdgeOfTetraChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
    // connect child
    int neighIndex;
    for (face = 0; face < 4; face++) {
      neighIndex = tetNeighOfFaceOfTetraChild[ch][face];
      child[ch]->setNeighbor(face,newNeighbor[neighIndex]);

    }
  }

  // assemble octa; don't forget barycenter
  for (vertex = 0; vertex < 6; vertex++) {
    vertexIndex = tetVertexOfOctaChild[vertex];
    child[4]->setVertex(vertex,newVertex[vertexIndex]);
  }
  for (edge = 0; edge < 12; edge++) {
    edgeIndex = tetEdgeOfOctaChild[edge];
    child[4]->setEdge(edge,newEdge[edgeIndex]);
  }
  for (face = 0; face < 8; face++) {
    neighIndex = tetNeighOfFaceOfOctaChild[face];
    child[4]->setNeighbor(face,newNeighbor[neighIndex]);
  }

  // connect external neighbor with children
  for (neighbor = 0; neighbor < 16; neighbor++) {

    if (newNeighbor[neighbor]) {
      childIndex = tetChildNeighOfNeighbor[neighbor];
      newNeighbor[neighbor]->setNeighbor(neighFace[neighbor],child[childIndex]);
    }
  }

  // set tag
  setRefined();
  // set no. of children
  chcount = 5;

  // set value
  for (ch = 0; ch < 5; ch++) child[ch]->setValue(getValue());
}



// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Tetrahedron::coarsen(GridLevel& gridLevel)
{
  int ch;
  int face;
  int nFace;
  Element* theNeighbor;

  // *- tets children
  for (ch = 0; ch < 4; ch++)
  {
    // *- set neighbor's pointer
    for (face = 0; face < 4; face++)
    {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor) {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeTetrahedron(child[ch]);
  }

  // *- octa child
  for (face = 0; face < 8; face++)
  {
    theNeighbor = child[4]->getNeighbor(face);
    if (theNeighbor) {
      nFace = theNeighbor->getCommonFace(child[4]);
      theNeighbor->setNeighbor(nFace,0);
    }
  }
  gridLevel.removeOctahedron(child[4]);

  // Delete the childre from memory
  for (ch = 0; ch < 5; ch++)
    delete child[ch];
  // Delete the pointers
  delete [] child;
  // Initialize pointer
  child = 0;

  // set tag
  setRegular();
  // set no. of children
  chcount = 0;

} // coarsen()


// Description
//
void
Tetrahedron::close(GridLevel& gridLevel)
{
  int refPattern = 0;
  int refMask    = 1;
  int noOfElem;
  int i0,i1,i2,i3;

  Vertex* aa;
  Vertex* bb;
  Vertex* cc;
  Vertex* dd;
  Vertex* ee;

  Vertex* v4;
  Vertex* v5;
  Vertex* v6;
  Vertex* v7;
  Vertex* v8;
  Vertex* v9;

  // Init element buffer
  Element* theElem[8];

  // get refinement pattern
  for (size_t i = 0; i < 6; i++)
  {
    if (e[i]->isRefined())
    {
      refPattern = (refPattern | refMask);
    }
    refMask = (refMask<<1);
  }

  switch(refPattern) {
    case 0:
      noOfElem = 0;
      break;
// case: one edge refined
    case 1:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v[1]);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,v[3]);
      noOfElem = 2;
      break;
    case 2:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v[1]);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,v[0]);
      theElem[1]->setVertex(1,aa);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,v[3]);
      noOfElem = 2;
      break;
    case 4:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[2]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,aa);
      theElem[0]->setVertex(2,v[2]);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,v[3]);
      noOfElem = 2;
      break;
    case 8:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[3]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v[1]);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,v[0]);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,aa);
      noOfElem = 2;
      break;
    case 16:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[4]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,aa);
      theElem[0]->setVertex(2,v[2]);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,v[0]);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,aa);
      noOfElem = 2;
      break;
    case 32:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;

      aa = e[5]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v[1]);
      theElem[0]->setVertex(2,v[2]);
      theElem[0]->setVertex(3,aa);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,v[3]);
      noOfElem = 2;
      break;
// case: two edges refined
    case 10:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      bb = e[3]->getMidPoint();
      if ((v[1]->getId()) > (v[3]->getId()))
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,bb);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,bb);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,bb);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      noOfElem = 3;
      break;
    case 18:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      bb = e[4]->getMidPoint();
      if ((v[2]->getId()) > v[3]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,bb);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 24:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[3]->getMidPoint();
      bb = e[4]->getMidPoint();
      if ((v[1]->getId()) > v[2]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,aa);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,aa);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 33:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[5]->getMidPoint();
      if ((v[2]->getId()) > v[3]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 9:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[3]->getMidPoint();
      if ((v[0]->getId()) > v[3]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,bb);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,bb);

        theElem[2]->setVertex(0,v[1]);
        theElem[2]->setVertex(1,v[2]);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,bb);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[1]);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      noOfElem = 3;
      break;
    case 40:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[3]->getMidPoint();
      bb = e[5]->getMidPoint();
      if ((v[0]->getId()) > v[2]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,aa);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,aa);
      }
      noOfElem = 3;
      break;
    case 20:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[2]->getMidPoint();
      bb = e[4]->getMidPoint();
      if ((v[0]->getId()) > v[3]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,bb);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      else {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,aa);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,bb);
      }
      noOfElem = 3;
      break;
    case 36:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[2]->getMidPoint();
      bb = e[5]->getMidPoint();
      if ((v[1]->getId()) > v[3]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,aa);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,bb);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,aa);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 48:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[4]->getMidPoint();
      bb = e[5]->getMidPoint();
      if ((v[0]->getId()) > v[1]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,aa);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,aa);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,bb);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,aa);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 3:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      if ((v[0]->getId()) > v[1]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,bb);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,bb);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 5:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      if ((v[1]->getId()) > v[2]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,v[1]);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
    case 6:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      if ((v[0]->getId()) > v[2]->getId())
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,v[1]);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,bb);
        theElem[2]->setVertex(1,aa);
        theElem[2]->setVertex(2,v[2]);
        theElem[2]->setVertex(3,v[3]);
      }
      noOfElem = 3;
      break;
// case: two oposite edges
    case 17:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[4]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,bb);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,v[0]);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,aa);
      theElem[1]->setVertex(3,bb);

      theElem[2]->setVertex(0,aa);
      theElem[2]->setVertex(1,bb);
      theElem[2]->setVertex(2,v[2]);
      theElem[2]->setVertex(3,v[3]);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,v[1]);
      theElem[3]->setVertex(2,v[2]);
      theElem[3]->setVertex(3,bb);
      noOfElem = 4;
      break;
    case 34:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      bb = e[5]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v[1]);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,bb);

      theElem[1]->setVertex(0,bb);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,aa);
      theElem[1]->setVertex(3,v[3]);

      theElem[2]->setVertex(0,v[0]);
      theElem[2]->setVertex(1,aa);
      theElem[2]->setVertex(2,v[2]);
      theElem[2]->setVertex(3,bb);

      theElem[3]->setVertex(0,bb);
      theElem[3]->setVertex(1,aa);
      theElem[3]->setVertex(2,v[2]);
      theElem[3]->setVertex(3,v[3]);
      noOfElem = 4;
      break;
    case 12:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[2]->getMidPoint();
      bb = e[3]->getMidPoint();
      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,aa);
      theElem[0]->setVertex(2,v[2]);
      theElem[0]->setVertex(3,bb);

      theElem[1]->setVertex(0,v[0]);
      theElem[1]->setVertex(1,aa);
      theElem[1]->setVertex(2,bb);
      theElem[1]->setVertex(3,v[3]);

      theElem[2]->setVertex(0,aa);
      theElem[2]->setVertex(1,v[1]);
      theElem[2]->setVertex(2,v[2]);
      theElem[2]->setVertex(3,bb);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,v[1]);
      theElem[3]->setVertex(2,bb);
      theElem[3]->setVertex(3,v[3]);
      noOfElem = 4;
      break;
// case: three edges on a face
    case 26:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[1]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();

      theElem[0]->setVertex(0,v[1]);
      theElem[0]->setVertex(1,cc);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[0]);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,bb);
      theElem[1]->setVertex(2,v[2]);
      theElem[1]->setVertex(3,v[0]);

      theElem[2]->setVertex(0,cc);
      theElem[2]->setVertex(1,v[3]);
      theElem[2]->setVertex(2,bb);
      theElem[2]->setVertex(3,v[0]);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,cc);
      theElem[3]->setVertex(2,bb);
      theElem[3]->setVertex(3,v[0]);
      noOfElem = 4;
      break;
    case 41:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[5]->getMidPoint();

      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,aa);
      theElem[0]->setVertex(2,cc);
      theElem[0]->setVertex(3,v[1]);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,v[2]);
      theElem[1]->setVertex(2,bb);
      theElem[1]->setVertex(3,v[1]);

      theElem[2]->setVertex(0,cc);
      theElem[2]->setVertex(1,bb);
      theElem[2]->setVertex(2,v[3]);
      theElem[2]->setVertex(3,v[1]);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,bb);
      theElem[3]->setVertex(2,cc);
      theElem[3]->setVertex(3,v[1]);
      noOfElem = 4;
      break;
    case 52:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[2]->getMidPoint();
      bb = e[4]->getMidPoint();
      cc = e[5]->getMidPoint();

      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,cc);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[2]);

      theElem[1]->setVertex(0,aa);
      theElem[1]->setVertex(1,bb);
      theElem[1]->setVertex(2,v[1]);
      theElem[1]->setVertex(3,v[2]);

      theElem[2]->setVertex(0,cc);
      theElem[2]->setVertex(1,v[3]);
      theElem[2]->setVertex(2,bb);
      theElem[2]->setVertex(3,v[2]);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,cc);
      theElem[3]->setVertex(2,bb);
      theElem[3]->setVertex(3,v[2]);
      noOfElem = 4;
      break;
    case 7:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();

      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,cc);
      theElem[0]->setVertex(2,aa);
      theElem[0]->setVertex(3,v[3]);

      theElem[1]->setVertex(0,cc);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,bb);
      theElem[1]->setVertex(3,v[3]);

      theElem[2]->setVertex(0,aa);
      theElem[2]->setVertex(1,bb);
      theElem[2]->setVertex(2,v[2]);
      theElem[2]->setVertex(3,v[3]);

      theElem[3]->setVertex(0,aa);
      theElem[3]->setVertex(1,cc);
      theElem[3]->setVertex(2,bb);
      theElem[3]->setVertex(3,v[3]);
      noOfElem = 4;
      break;
// case: three edges not on the same face
    case 11:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[3]->getMidPoint();

      if (i0 > i1 && i0 > i3)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
      }
      else if (i1 > i0 && i1 > i3)
      {
        if (i0 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
      }
      else if (i3 > i0 && i3 > i1)
      {
        if (i0 > i1)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);
        }
      }
      noOfElem = 4;
      break;
    case 13:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      else //if (i3 > i0)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      noOfElem = 5;
      break;
    case 14:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();

      if (i0 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      else //if (i2 > i0)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      noOfElem = 5;
      break;
    case 19:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[4]->getMidPoint();

      if (i0 > i1)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      else //if (i1 > i0)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      noOfElem = 5;
      break;
    case 21:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      else //if (i3 > i0)
      {
        if (i1 > i2) {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      noOfElem = 5;
      break;
    case 22:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      i0 = v[0]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[4]->getMidPoint();

      if (i0 > i2 && i0 > i3)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);
        }
      }
      else if (i2 > i0 && i2 > i3)
      {
        if (i0 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);
        }
      }
      else if (i3 > i0 && i3 > i2)
      {
        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      noOfElem = 4;
      break;
    case 25:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      else //if (i3 > i0)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      noOfElem = 5;
      break;
    case 28:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[2]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[2]);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,cc);
        }
      }
      else //if (i3 > i0)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      noOfElem = 5;
      break;
    case 35:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      else //if (i1 > i0)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      noOfElem = 5;
      break;
    case 37:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i1 > i2 && i1 > i3)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      else if (i2 > i1 && i2 > i3)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      else if (i3 > i1 && i3 > i2)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[1]);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      noOfElem = 4;
      break;
    case 38:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      else //if (i2 > i0)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      noOfElem = 5;
      break;
    case 42:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      else //if (i2 > i0)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      noOfElem = 5;
      break;
    case 44:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[2]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      else //if (i2 > i0)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      noOfElem = 5;
      break;
    case 49:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[4]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,cc);
          theElem[0]->setVertex(1,v[0]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,cc);
          theElem[0]->setVertex(1,v[0]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
      }
      else //if (i1 > i0)
      {
        if (i2 > i3) {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[2]);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);
        }
      }
      noOfElem = 5;
      break;
    case 50:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[4]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[0]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[0]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      else //if (i1 > i0)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);
        }
      }
      noOfElem = 5;
      break;
    case 56:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();

      aa = e[3]->getMidPoint();
      bb = e[4]->getMidPoint();
      cc = e[5]->getMidPoint();

      if (i0 > i1 && i0 > i2)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[0]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,aa);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,aa);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[0]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      else if (i1 > i0 && i1 > i2)
      {
        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,aa);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,aa);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      else if (i2 > i0 && i2 > i1)
      {
        if (i0 > i1)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,bb);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[0]);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,aa);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,cc);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,aa);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);
        }
      }
      noOfElem = 4;
      break;
// case: four edges refined
    case 15:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[3]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,dd);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,dd);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
      }
      else
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,dd);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,dd);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,cc);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
      }
      noOfElem = 6;
      break;
    case 23:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,cc);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,dd);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,cc);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,dd);
        }
      }
      else
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,cc);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,dd);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,cc);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,dd);
        }
      }
      noOfElem = 6;
      break;
    case 27:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();

      if (i0 > i1)
      {
        if (i0 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
      }
      else
      {
        if (i0 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,v[0]);
          theElem[5]->setVertex(1,dd);
          theElem[5]->setVertex(2,cc);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,dd);
          theElem[5]->setVertex(2,cc);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      noOfElem = 6;
      break;
    case 29:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[1]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
      }
      else
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,v[1]);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,dd);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,dd);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,dd);
        }
      }
      noOfElem = 6;
      break;
    case 30:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();

      if (i0 > i3)
      {
        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,dd);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,v[3]);

          noOfElem = 5;
        }
        else
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,aa);
          theElem[5]->setVertex(3,dd);

          noOfElem = 6;
        }
      }
      else
      {
        theElem[5] = new Tetrahedron;

        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,aa);
          theElem[5]->setVertex(3,dd);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,v[3]);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,dd);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,aa);
          theElem[5]->setVertex(3,dd);
        }
        noOfElem = 6;
      }
      break;
    case 39:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i1 > i3)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,dd);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      else
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,cc);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      noOfElem = 6;
      break;
    case 43:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i1)
      {
        theElem[5] = new Tetrahedron;

        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[0]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[1]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,cc);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,dd);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,dd);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,cc);
          theElem[5]->setVertex(3,v[3]);
        }
        noOfElem = 6;
      }
      else
      {
        if (i1 > i3)
        {
          theElem[0] = new Tetrahedron;
          theElem[1] = new Tetrahedron;
          theElem[2] = new Tetrahedron;
          theElem[3] = new Tetrahedron;
          theElem[4] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,v[3]);

          noOfElem = 5;
        }
        else
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,dd);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);

          noOfElem = 6;
        }
      }
      break;
    case 45:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i1 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,v[1]);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          noOfElem = 5;
        }
        else
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[1]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);

          noOfElem = 6;
        }
      }
      else
      {
        theElem[5] = new Tetrahedron;

        if (i1 > i3) {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,cc);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,dd);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,cc);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        noOfElem = 6;
      }
      break;
    case 46:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i2)
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,dd);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,v[0]);
          theElem[5]->setVertex(1,aa);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,cc);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,bb);
          theElem[5]->setVertex(1,v[1]);
          theElem[5]->setVertex(2,aa);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      else
      {
        if (i1 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,cc);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,dd);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,v[3]);

          theElem[2]->setVertex(0,cc);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,v[3]);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,cc);
        }
      }
      noOfElem = 6;
      break;
    case 51:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      else
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,bb);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      noOfElem = 6;
      break;
    case 53:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i1 = v[1]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i1 > i2)
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[1]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      else
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,bb);
          theElem[1]->setVertex(2,cc);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,bb);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,bb);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,aa);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      noOfElem = 6;
      break;
    case 54:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i2)
      {
        theElem[5] = new Tetrahedron;

        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,v[3]);

          theElem[4]->setVertex(0,v[0]);
          theElem[4]->setVertex(1,aa);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,aa);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);
        }
        noOfElem = 6;
      }
      else
      {
        if (i2 > i3)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,v[3]);

          noOfElem = 5;
        }
        else
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,bb);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,bb);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,bb);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,bb);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,v[3]);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,aa);
          theElem[5]->setVertex(2,v[2]);
          theElem[5]->setVertex(3,v[3]);

          noOfElem = 6;
        }
      }
      break;
    case 57:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();

      aa = e[0]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,cc);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,v[1]);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,cc);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      else
      {
        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,cc);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
      }
      noOfElem = 6;
      break;
    case 58:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();

      aa = e[1]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i1)
      {
        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          noOfElem = 5;
        }
        else
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,cc);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,aa);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,v[0]);
          theElem[3]->setVertex(1,aa);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);

          noOfElem = 6;
        }
      }
      else
      {
        theElem[5] = new Tetrahedron;

        if (i0 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[0]);
          theElem[2]->setVertex(1,aa);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,bb);

          theElem[3]->setVertex(0,v[1]);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,aa);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,aa);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,v[1]);
          theElem[0]->setVertex(2,aa);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,v[1]);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,aa);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,dd);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,bb);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,dd);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        noOfElem = 6;
      }
      break;
    case 60:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();
      i2 = v[2]->getId();

      aa = e[2]->getMidPoint();
      bb = e[3]->getMidPoint();
      cc = e[4]->getMidPoint();
      dd = e[5]->getMidPoint();

      if (i0 > i2)
      {
        theElem[5] = new Tetrahedron;

        if (i1 > i2)
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,bb);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,v[0]);
          theElem[1]->setVertex(1,aa);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,bb);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,cc);
          theElem[2]->setVertex(2,bb);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[2]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,cc);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);
        }
        noOfElem = 6;
      }
      else
      {
        if (i1 > i2)
        {
          theElem[5] = new Tetrahedron;

          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,bb);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,bb);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,dd);

          theElem[3]->setVertex(0,aa);
          theElem[3]->setVertex(1,v[1]);
          theElem[3]->setVertex(2,bb);
          theElem[3]->setVertex(3,cc);

          theElem[4]->setVertex(0,aa);
          theElem[4]->setVertex(1,v[1]);
          theElem[4]->setVertex(2,v[2]);
          theElem[4]->setVertex(3,bb);

          theElem[5]->setVertex(0,dd);
          theElem[5]->setVertex(1,cc);
          theElem[5]->setVertex(2,bb);
          theElem[5]->setVertex(3,v[3]);

          noOfElem = 6;
        }
        else
        {
          theElem[0]->setVertex(0,v[0]);
          theElem[0]->setVertex(1,aa);
          theElem[0]->setVertex(2,v[2]);
          theElem[0]->setVertex(3,dd);

          theElem[1]->setVertex(0,aa);
          theElem[1]->setVertex(1,cc);
          theElem[1]->setVertex(2,v[2]);
          theElem[1]->setVertex(3,dd);

          theElem[2]->setVertex(0,aa);
          theElem[2]->setVertex(1,v[1]);
          theElem[2]->setVertex(2,v[2]);
          theElem[2]->setVertex(3,cc);

          theElem[3]->setVertex(0,dd);
          theElem[3]->setVertex(1,cc);
          theElem[3]->setVertex(2,v[2]);
          theElem[3]->setVertex(3,bb);

          theElem[4]->setVertex(0,dd);
          theElem[4]->setVertex(1,cc);
          theElem[4]->setVertex(2,bb);
          theElem[4]->setVertex(3,v[3]);

          noOfElem = 5;
        }
      }
      break;
    case 31:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i0 = v[0]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[3]->getMidPoint();
      ee = e[4]->getMidPoint();

      if (i0 > i3)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,v[0]);
        theElem[1]->setVertex(1,ee);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,dd);

        theElem[2]->setVertex(0,v[0]);
        theElem[2]->setVertex(1,ee);
        theElem[2]->setVertex(2,dd);
        theElem[2]->setVertex(3,v[3]);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,cc);
        theElem[3]->setVertex(2,bb);
        theElem[3]->setVertex(3,ee);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,ee);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,cc);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,bb);
        theElem[5]->setVertex(3,ee);

        theElem[6]->setVertex(0,aa);
        theElem[6]->setVertex(1,bb);
        theElem[6]->setVertex(2,v[2]);
        theElem[6]->setVertex(3,dd);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,v[3]);

        theElem[1]->setVertex(0,cc);
        theElem[1]->setVertex(1,ee);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,v[3]);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,ee);
        theElem[2]->setVertex(2,dd);
        theElem[2]->setVertex(3,v[3]);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,cc);
        theElem[3]->setVertex(2,bb);
        theElem[3]->setVertex(3,ee);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,ee);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,cc);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,bb);
        theElem[5]->setVertex(3,ee);

        theElem[6]->setVertex(0,aa);
        theElem[6]->setVertex(1,bb);
        theElem[6]->setVertex(2,v[2]);
        theElem[6]->setVertex(3,dd);
      }
      noOfElem = 7;
      break;
    case 47:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i1 = v[1]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[3]->getMidPoint();
      ee = e[5]->getMidPoint();

      if (i1 > i3)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,cc);
        theElem[1]->setVertex(2,bb);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,ee);
        theElem[2]->setVertex(2,bb);
        theElem[2]->setVertex(3,dd);

        theElem[3]->setVertex(0,cc);
        theElem[3]->setVertex(1,v[1]);
        theElem[3]->setVertex(2,bb);
        theElem[3]->setVertex(3,ee);

        theElem[4]->setVertex(0,ee);
        theElem[4]->setVertex(1,v[1]);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,aa);
        theElem[5]->setVertex(1,bb);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,dd);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,v[1]);
        theElem[6]->setVertex(2,dd);
        theElem[6]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,cc);
        theElem[1]->setVertex(2,bb);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,ee);
        theElem[2]->setVertex(2,bb);
        theElem[2]->setVertex(3,dd);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,bb);
        theElem[3]->setVertex(2,v[2]);
        theElem[3]->setVertex(3,dd);

        theElem[4]->setVertex(0,ee);
        theElem[4]->setVertex(1,cc);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,v[3]);

        theElem[5]->setVertex(0,ee);
        theElem[5]->setVertex(1,bb);
        theElem[5]->setVertex(2,dd);
        theElem[5]->setVertex(3,v[3]);

        theElem[6]->setVertex(0,cc);
        theElem[6]->setVertex(1,v[1]);
        theElem[6]->setVertex(2,bb);
        theElem[6]->setVertex(3,v[3]);
      }
      noOfElem = 7;
      break;
    case 55:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i2 = v[2]->getId();
      i3 = v[3]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();
      dd = e[4]->getMidPoint();
      ee = e[5]->getMidPoint();

      if (i2 > i3)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,cc);
        theElem[1]->setVertex(2,dd);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,cc);
        theElem[2]->setVertex(2,bb);
        theElem[2]->setVertex(3,dd);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,dd);
        theElem[3]->setVertex(2,v[2]);
        theElem[3]->setVertex(3,ee);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,bb);
        theElem[4]->setVertex(2,v[2]);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,cc);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,bb);
        theElem[5]->setVertex(3,dd);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,v[2]);
        theElem[6]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,cc);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,cc);
        theElem[1]->setVertex(2,dd);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,cc);
        theElem[2]->setVertex(2,bb);
        theElem[2]->setVertex(3,dd);

        theElem[3]->setVertex(0,cc);
        theElem[3]->setVertex(1,v[1]);
        theElem[3]->setVertex(2,bb);
        theElem[3]->setVertex(3,dd);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,bb);
        theElem[4]->setVertex(2,v[2]);
        theElem[4]->setVertex(3,v[3]);

        theElem[5]->setVertex(0,ee);
        theElem[5]->setVertex(1,dd);
        theElem[5]->setVertex(2,aa);
        theElem[5]->setVertex(3,v[3]);

        theElem[6]->setVertex(0,aa);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,bb);
        theElem[6]->setVertex(3,v[3]);
      }
      noOfElem = 7;
      break;
    case 59:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i0 = v[0]->getId();
      i1 = v[1]->getId();

      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();
      ee = e[5]->getMidPoint();

      if (i0 > i1)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,dd);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,dd);
        theElem[1]->setVertex(2,cc);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,dd);
        theElem[2]->setVertex(2,bb);
        theElem[2]->setVertex(3,cc);

        theElem[3]->setVertex(0,v[0]);
        theElem[3]->setVertex(1,bb);
        theElem[3]->setVertex(2,aa);
        theElem[3]->setVertex(3,dd);

        theElem[4]->setVertex(0,v[0]);
        theElem[4]->setVertex(1,v[1]);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,aa);
        theElem[5]->setVertex(1,bb);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,cc);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,v[1]);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,v[1]);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,dd);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,dd);
        theElem[2]->setVertex(2,cc);
        theElem[2]->setVertex(3,ee);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,dd);
        theElem[3]->setVertex(2,bb);
        theElem[3]->setVertex(3,cc);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,bb);
        theElem[4]->setVertex(2,v[2]);
        theElem[4]->setVertex(3,cc);

        theElem[5]->setVertex(0,aa);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,dd);
        theElem[5]->setVertex(3,ee);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      noOfElem = 7;
      break;
    case 61:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i1 = v[1]->getId();
      i2 = v[2]->getId();

      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();
      ee = e[5]->getMidPoint();

      if (i1 > i2)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,dd);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,dd);
        theElem[2]->setVertex(2,cc);
        theElem[2]->setVertex(3,ee);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,bb);
        theElem[3]->setVertex(2,v[1]);
        theElem[3]->setVertex(3,dd);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,v[1]);
        theElem[4]->setVertex(2,cc);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,aa);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,cc);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,aa);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,dd);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,aa);
        theElem[2]->setVertex(1,dd);
        theElem[2]->setVertex(2,cc);
        theElem[2]->setVertex(3,ee);

        theElem[3]->setVertex(0,aa);
        theElem[3]->setVertex(1,cc);
        theElem[3]->setVertex(2,dd);
        theElem[3]->setVertex(3,v[2]);

        theElem[4]->setVertex(0,aa);
        theElem[4]->setVertex(1,dd);
        theElem[4]->setVertex(2,bb);
        theElem[4]->setVertex(3,v[2]);

        theElem[5]->setVertex(0,bb);
        theElem[5]->setVertex(1,v[1]);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,dd);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      noOfElem = 7;
      break;
    case 62:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;

      i0 = v[0]->getId();
      i2 = v[2]->getId();

      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      cc = e[3]->getMidPoint();
      dd = e[4]->getMidPoint();
      ee = e[5]->getMidPoint();

      if (i0 > i2)
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,aa);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,ee);
        theElem[1]->setVertex(1,bb);
        theElem[1]->setVertex(2,aa);
        theElem[1]->setVertex(3,dd);

        theElem[2]->setVertex(0,ee);
        theElem[2]->setVertex(1,dd);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,cc);

        theElem[3]->setVertex(0,v[0]);
        theElem[3]->setVertex(1,aa);
        theElem[3]->setVertex(2,cc);
        theElem[3]->setVertex(3,ee);

        theElem[4]->setVertex(0,bb);
        theElem[4]->setVertex(1,v[1]);
        theElem[4]->setVertex(2,aa);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,v[0]);
        theElem[5]->setVertex(1,aa);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,cc);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      else
      {
        theElem[0]->setVertex(0,v[0]);
        theElem[0]->setVertex(1,bb);
        theElem[0]->setVertex(2,v[2]);
        theElem[0]->setVertex(3,ee);

        theElem[1]->setVertex(0,bb);
        theElem[1]->setVertex(1,aa);
        theElem[1]->setVertex(2,v[2]);
        theElem[1]->setVertex(3,ee);

        theElem[2]->setVertex(0,ee);
        theElem[2]->setVertex(1,bb);
        theElem[2]->setVertex(2,aa);
        theElem[2]->setVertex(3,dd);

        theElem[3]->setVertex(0,ee);
        theElem[3]->setVertex(1,dd);
        theElem[3]->setVertex(2,aa);
        theElem[3]->setVertex(3,cc);

        theElem[4]->setVertex(0,bb);
        theElem[4]->setVertex(1,v[1]);
        theElem[4]->setVertex(2,aa);
        theElem[4]->setVertex(3,dd);

        theElem[5]->setVertex(0,ee);
        theElem[5]->setVertex(1,aa);
        theElem[5]->setVertex(2,v[2]);
        theElem[5]->setVertex(3,cc);

        theElem[6]->setVertex(0,ee);
        theElem[6]->setVertex(1,dd);
        theElem[6]->setVertex(2,cc);
        theElem[6]->setVertex(3,v[3]);
      }
      noOfElem = 7;
      break;
// case: full pattern
    case 63:
      theElem[0] = new Tetrahedron;
      theElem[1] = new Tetrahedron;
      theElem[2] = new Tetrahedron;
      theElem[3] = new Tetrahedron;
      theElem[4] = new Tetrahedron;
      theElem[5] = new Tetrahedron;
      theElem[6] = new Tetrahedron;
      theElem[7] = new Tetrahedron;

      v4 = e[0]->getMidPoint();
      v5 = e[1]->getMidPoint();
      v6 = e[2]->getMidPoint();
      v7 = e[3]->getMidPoint();
      v8 = e[4]->getMidPoint();
      v9 = e[5]->getMidPoint();

      theElem[0]->setVertex(0,v[0]);
      theElem[0]->setVertex(1,v6);
      theElem[0]->setVertex(2,v4);
      theElem[0]->setVertex(3,v9);

      theElem[1]->setVertex(0,v6);
      theElem[1]->setVertex(1,v[1]);
      theElem[1]->setVertex(2,v5);
      theElem[1]->setVertex(3,v8);

      theElem[2]->setVertex(0,v4);
      theElem[2]->setVertex(1,v5);
      theElem[2]->setVertex(2,v[2]);
      theElem[2]->setVertex(3,v7);

      theElem[3]->setVertex(0,v9);
      theElem[3]->setVertex(1,v8);
      theElem[3]->setVertex(2,v7);
      theElem[3]->setVertex(3,v[3]);

      theElem[4]->setVertex(0,v9);
      theElem[4]->setVertex(1,v6);
      theElem[4]->setVertex(2,v4);
      theElem[4]->setVertex(3,v8);

      theElem[5]->setVertex(0,v9);
      theElem[5]->setVertex(1,v4);
      theElem[5]->setVertex(2,v7);
      theElem[5]->setVertex(3,v8);

      theElem[6]->setVertex(0,v4);
      theElem[6]->setVertex(1,v6);
      theElem[6]->setVertex(2,v5);
      theElem[6]->setVertex(3,v8);

      theElem[7]->setVertex(0,v4);
      theElem[7]->setVertex(1,v5);
      theElem[7]->setVertex(2,v7);
      theElem[7]->setVertex(3,v8);
      noOfElem = 8;
      break;
    default:
      noOfElem = 0;
      break;
  }

  // Connect hierarchy
  child = new Element*[noOfElem];
  for (int nn = 0; nn < noOfElem; nn++)
  {
    theElem[nn]->setClosure();
    theElem[nn]->setValue(getValue());
    theElem[nn]->setParent(this);
    child[nn] = theElem[nn];
    gridLevel.appendTetrahedron(theElem[nn]);
  }

  // set no. of children
  chcount = noOfElem;
  // set irregular mark
  setIrregular();
  // set neighbors and create new edges
  connect(gridLevel);
}



// Description
//  computes the edge refinement pattern
int
Tetrahedron::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;
  for (int i = 0; i < 6; i++) {
    if (e[i]->isRefined()) {
      refPattern = (refPattern | refMask);
    }
    refMask = (refMask<<1);
  }
  return refPattern;
}


//*************************************************************
// ____________________________________________________________
//
//
//         virtual functions
//
// ____________________________________________________________
//*************************************************************


//=============================================================================
//  TOPOLOGY
//=============================================================================

// Description
//
bool
Tetrahedron::queryEdgeRefinement() const
{
  for (int i = 0; i < 6; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Tetrahedron::queryChildrenEdgeRefinement() const
{
  int i;
  for (i = 0; i < 6; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  // check if neighbor's child is refined
  for (i = 0; i < 4; i++) {
    if (n[i]) {
      if (n[i]->queryChildrenRefinement())
        return true;
    }
  }

  return false;
}


// Description
int
Tetrahedron::getCommonFace(Element* theNeighbor) const
{
  for (int face = 0; face < 4; face++)
    if (n[face] == theNeighbor)
      return face;

  return -1;
}



// Description
//
void
Tetrahedron::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
{
  int vertex,index, childIndex;

  for (vertex = 0; vertex < 3; vertex++) {
    index = tetVertexOfFace[face][vertex];
    theVertex[vertex] = v[index];

    childIndex = tetChildOfParentFace[face][vertex];
    theChild[vertex] = child[childIndex];
    theChildFace[vertex] = tetChildFaceOfParentFace[face][vertex];

    index = tetChildEdgeOfParentFace[face][vertex];
    theEdge[vertex] = child[childIndex]->getEdge(index);
  }

  childIndex = tetChildOfParentFace[face][3];
  theChild[3] = child[childIndex];
  theChildFace[3] = tetChildFaceOfParentFace[face][3];
}



// Description
//
void
Tetrahedron::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
void
Tetrahedron::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 4; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Tetrahedron::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Tetrahedron::getVertices(Vertex** vt) const
{
  for (int i = 0; i < 4; i++)
    vt[i] = v[i];
}


// Description
//
void
Tetrahedron::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
void
Tetrahedron::setEdges(Edge** edg)
{
  for (int i = 0; i < 6; i++)
    e[i] = edg[i];
}

// Description
//
Edge*
Tetrahedron::getEdge(int i) const
{
  return e[i];
}


// Description
//
void
Tetrahedron::getEdges(Edge** edg) const
{
  for (int i = 0; i < 6; i++)
    edg[i] = e[i];
}



// Description
//
void
Tetrahedron::setNeighbor(int face,Element* theNeighbor)
{
  n[face] = theNeighbor;
}


// Description
//
Element*
Tetrahedron::getNeighbor(int face) const
{
  return n[face];
}


//=============================================================================
//  TOPOLOGY
//=============================================================================

// Description
//
int
Tetrahedron::type() const
{
  return GRD_TETRAHEDRON;
}

// Description
//
int
Tetrahedron::noOfVerticesAtFace(int fc) const
{
  return 3;
}


// Description
//
int
Tetrahedron::topologicalVertexAtFace(int fc,int vt) const
{
  return tetVertexOfFace[fc][vt];
}


// Description
//
Vertex*
Tetrahedron::getVertexAtFace(int fc,int vt) const
{
  return (v[tetVertexOfFace[fc][vt]]);
}

// Description
//
int
Tetrahedron::noOfEdgesAtFace(int fc) const
{
  return 3;
}
// Description
//
Vertex*
Tetrahedron::getVertexAtEdge(int ed,int vt) const
{
  return (v[tetVertexOfEdge[ed][vt]]);
}

//
//
Edge*
Tetrahedron::getEdgeAtFace(int fc,int ed) const
{
  return e[tetEdgeOfFace[fc][ed]];
}


// Description
//
void
Tetrahedron::connect(GridLevel& gridLevel)
{
  // Create List of objects to be connected
  std::list<Element*> nelt;
  std::list<Element*> nelc;

  // Collect all neighbors
  for (int nn = 0; nn < 4; nn++)
  {
    if (n[nn]->getNoOfChildren() > 0)
    {
      for (int ch = 0; ch < n[nn]->getNoOfChildren(); ch++)
        nelt.push_back(n[nn]->getChild(ch));
    }
    else
    {
      if (!n[nn]->isIrregular())
        nelt.push_back(n[nn]);
    }
  }

  // Connect elements in the list
  typedef std::list<Element*>::iterator EI;
  for (int ch = 0; ch < static_cast<int>(chcount); ch++)
  {
    for (int nn = 0; nn < child[ch]->getNoOfFaces(); ++nn)
    {
      Element* neighbor = child[ch]->getNeighbor(nn);
      if (neighbor)
        continue;
      else
      {
        Vertex* vt[3];
        vt[0] = child[ch]->getVertexAtFace(nn,0);
        vt[1] = child[ch]->getVertexAtFace(nn,1);
        vt[2] = child[ch]->getVertexAtFace(nn,2);

        // connect children
        for (int chq = ch+1; chq < static_cast<int>(chcount); chq++)
        {
          for (int fc = 0; fc < child[chq]->getNoOfFaces(); fc++)
          {
            Vertex* vv[3];
            vv[0] = child[chq]->getVertexAtFace(fc,0);
            vv[1] = child[chq]->getVertexAtFace(fc,1);
            vv[2] = child[chq]->getVertexAtFace(fc,2);

            // compare vertices
            int count = 0;
            for (int i = 0; i < 3; i++)
            {
              for (int j = 0; j < 3; j++)
              {
                if (vt[i] == vv[j])
                {
                  count++;
                  break;
                }
              }
            }
            if (count == 3)
            {
              // neighbour found
              child[ch]->setNeighbor(nn,child[chq]);
              child[chq]->setNeighbor(fc,child[ch]);
              break;
            }
          }
        } // for q child list
        // connect children
        EI q = nelt.begin();
        while (q != nelt.end())
        {
          bool flag = false;
          for (int fc = 0; fc < (*q)->getNoOfFaces(); fc++)
          {
            Vertex* vv[3];
            vv[0] = (*q)->getVertexAtFace(fc,0);
            vv[1] = (*q)->getVertexAtFace(fc,1);
            vv[2] = (*q)->getVertexAtFace(fc,2);

            // compare vertices
            int count = 0;
            for (int i = 0; i < 3; i++)
            {
              for (int j = 0; j < 3; j++)
              {
                if (vt[i] == vv[j])
                {
                  count++;
                  break;
                }
              }
            }
            if (count == 3)
            {
              // neighbour found
              child[ch]->setNeighbor(nn,*q);
              (*q)->setNeighbor(fc,child[ch]);
              flag = true;
              break;
            }  
          }
          if (flag)
          {
            nelc.push_back(*q);
            q = nelt.erase(q);
          }
          else ++q;
        } // for q neighbor list
      } // else neighbor not connected
    } // for nn faces
  } // for p

  // Check edges
  // Create List of objects to be connected
  typedef std::list<Edge*>::iterator EGI;
  std::list<Edge*> edglt;
  for (int nn = 0; nn < 6; nn++)
  {
    if (!e[nn]->isRefined())
      edglt.push_back(e[nn]);
    else
    {
      edglt.push_back(e[nn]->getChild(0));
      edglt.push_back(e[nn]->getChild(1));
    }
  }
  // get edge from neighbors
  for (EI q = nelc.begin(); q != nelc.end(); ++q)
  {
    for (int nn = 0; nn < (*q)->getNoOfEdges(); nn++)
    {
      Edge* tmp = (*q)->getEdge(nn);
      if (tmp) edglt.push_back(tmp);
    }
  }

  // create new edges
  for (int nn = 0; nn < static_cast<int>(chcount); nn++)
  {
    for (int i = 0; i < child[nn]->getNoOfEdges(); i++)
    {
      Edge* edge = child[nn]->getEdge(i);
      if (edge == 0)
      {
        Vertex* vt0 = child[nn]->getVertexAtEdge(i,0);
        Vertex* vt1 = child[nn]->getVertexAtEdge(i,1);
        bool flag = false;
        for (EGI ei = edglt.begin(); ei != edglt.end(); ++ei)
        {
          Vertex* ev0 = (*ei)->getVertex(0);
          Vertex* ev1 = (*ei)->getVertex(1);
          if ((vt0 == ev0 && vt1 == ev1) || (vt1 == ev0 && vt0 == ev1))
          {
            child[nn]->setEdge(i,*ei);
            flag = true;
            break;
          }
        } // for EI p

        // check if edge found
        // if not, create a new edge
        if (!flag)
        {
          Edge* tmp = new Edge(vt0,vt1);
          child[nn]->setEdge(i,tmp);
          edglt.push_back(tmp);
          gridLevel.appendEdge(tmp);
        }
      } // if (edge == 0)
    } // for i
  } // for nn chcount


  // clear list
  nelt.clear();
  nelc.clear();
  edglt.clear();
} // connect()

//// Description
////
//void
//Tetrahedron::getGlobalCoords(const double* lc,double* coords)
//{
//  double* pos[4];
//
//  double l1 = lc[0];
//  double l2 = lc[1];
//  double l3 = lc[2];
//  double l0 = 1.0 - l1 - l2 - l3;
//
//
//  pos[0] = v[0]->getPosition();
//  pos[1] = v[1]->getPosition();
//  pos[2] = v[2]->getPosition();
//  pos[3] = v[3]->getPosition();
//
//  coords[0] = l0*pos[0][0] + l1*pos[1][0] + l2*pos[2][0] + l3*pos[3][0];
//  coords[1] = l0*pos[0][1] + l1*pos[1][1] + l2*pos[2][1] + l3*pos[3][1];
//  coords[2] = l0*pos[0][2] + l1*pos[1][2] + l2*pos[2][2] + l3*pos[3][2];
//
//}
//
//
//// Description
////
//void
//Tetrahedron::getJacobian(Jacobian& Jac)
//{
//  Jac(0,0) = (*v[1])[0] - (*v[0])[0];
//  Jac(0,1) = (*v[2])[0] - (*v[0])[0];
//  Jac(0,2) = (*v[3])[0] - (*v[0])[0];
//
//  Jac(1,0) = (*v[1])[1] - (*v[0])[1];
//  Jac(1,1) = (*v[2])[1] - (*v[0])[1];
//  Jac(1,2) = (*v[3])[1] - (*v[0])[1];
//
//  Jac(2,0) = (*v[1])[2] - (*v[0])[2];
//  Jac(2,1) = (*v[2])[2] - (*v[0])[2];
//  Jac(2,2) = (*v[3])[2] - (*v[0])[2];
//
//}



} // namespace grd
