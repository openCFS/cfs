// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/***************************************************************************
    File        : Quadrangle.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 8 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Quadrangle.hh"
#include "Triangle.hh"
#include "GridLevel.hh"

#include "Utility.hh"

namespace grd {

// (static obj) initialize memory management here
//Pool<Quadrangle> Quadrangle::pool;

// Description:
//
Quadrangle::Quadrangle( )
{
  int i;
  // Vertices and edges
  for (i = 0; i < 4; i++) {
    v[i] = 0;
    e[i] = 0;
  }
  // Pointer to volume element
  n = 0;

  // set no. of children
  // chcount = 0; <- set in abstract class Element
}


//=============================================================================
//  REFINEMENT
//=============================================================================


// Description
// Consider refinement of isotropic and anisotropic
// elements
void
Quadrangle::refine(GridLevel& gridLevel)
{
  int i,ch;
  int vertCount,edgCount;
  int vertex,edge;
  int vertexIndex;

  int index;
  int offset;
  int nFace;

  int theChildFace[4];

  Vertex*  newVertex[9];
  Edge*    newEdge[12];
  Vertex*  theVertex[5];
  Edge*    theEdge[4];
  Element* theChild[4];

  // -* init indices
  offset = 0;
  index = -1;
  nFace = -1;

  // -* init counters
  vertCount = edgCount = 0;

  // -* 1. process parent vertices
  for (vertex = 0; vertex < 4; vertex++)
    newVertex[vertCount++] = v[vertex];

  // -* 2. process parent edges
  for (edge = 0; edge < 4; edge++) {
    if (e[edge]->isRefined()) {
      e[edge]->getChildren(newVertex[quaVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[quaVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      newVertex[vertCount] = e[edge]->getMidPoint();

      // append new elements in list
      gridLevel.appendVertex(newVertex[vertCount]);
      gridLevel.appendEdge(newEdge[edgCount]);
      gridLevel.appendEdge(newEdge[edgCount+1]);

      edgCount += 2;
      vertCount++;
    }
  }


  // create interior edges

  // first check if volume elment is refined
  if (n) {
    if (n->isRefined()) {
      nFace = n->getCommonFace(this);
      n->getRefElements(nFace,theVertex,theEdge,theChild,theChildFace);

      if (v[0] == theVertex[0])
        offset = 0;
      else if (v[1] == theVertex[0])
        offset = 1;
      else if (v[2] == theVertex[0])
        offset = 2;
      else if (v[3] == theVertex[0])
        offset = 3;
      else {
        // Panic
        Error("can't find common vertex from Quad face");
      }
       
      //index = (4 - offset) % 4;
      index = (5 - offset) % 4;
      newEdge[8] = theEdge[index];

      //index = (5 - offset) % 4;
      //newEdge[11] = theEdge[index];
      index = (6 - offset) % 4;
      newEdge[9] = theEdge[index];
      
      //index = (6 - offset) % 4;
      index = (7 - offset) % 4;
      newEdge[10] = theEdge[index];

      //index = (7 - offset) % 4;
      //newEdge[9] = theEdge[index];
      index = (8 - offset) % 4;
      newEdge[11] = theEdge[index];

      Vertex* vv0 = newEdge[8]->getVertex(0);
      Vertex* vv1 = newEdge[8]->getVertex(1);
      if (!(vv0 == newVertex[4] && vv1 == theVertex[4]) && !(vv1 == newVertex[4] && vv0 == theVertex[4]))
        std::cout << "ERROR: edge wrong vertex no. interpretation at face: " << nFace << std::endl;
      vv0 = newEdge[10]->getVertex(0);
      vv1 = newEdge[10]->getVertex(1);
      if (!(vv0 == newVertex[6] && vv1 == theVertex[4]) && !(vv1 == newVertex[6] && vv0 == theVertex[4]))
        std::cout << "ERROR: edge wrong vertex no. interpretation at face: " << nFace << std::endl;
      vv0 = newEdge[9]->getVertex(0);
      vv1 = newEdge[9]->getVertex(1);
      if (!(vv0 == newVertex[5] && vv1 == theVertex[4]) && !(vv1 == newVertex[5] && vv0 == theVertex[4]))
        std::cout << "ERROR: edge wrong vertex no. interpretation at face: " << nFace << std::endl;
      vv0 = newEdge[11]->getVertex(0);
      vv1 = newEdge[11]->getVertex(1);
      if (!(vv0 == newVertex[7] && vv1 == theVertex[4]) && !(vv1 == newVertex[7] && vv0 == theVertex[4]))
        std::cout << "ERROR: edge wrong vertex no. interpretation at face: " << nFace << std::endl;
      
      // vertex at face center
      newVertex[8] = theVertex[4];
    } // if volume neighbor is refined
    else {
      // create vertex at face center
      double pos[3];
      for (i = 0; i < 3; i++) {
        pos[i] = 0.0;
        for (vertex = 0; vertex < 4; vertex++) {
          pos[i] += (*v[vertex])[i];
        }
        pos[i] *= 0.25;
      }
      Vertex* tmpVt = new Vertex(pos);
      gridLevel.appendVertex(tmpVt);
      newVertex[8] = tmpVt;

      // create and set internal edges
      Edge* e8 = new Edge;
      e8->setVertices(newVertex[4],newVertex[8]);
      Edge* e9 = new Edge;
      e9->setVertices(newVertex[5],newVertex[8]);
      Edge* e10 = new Edge;
      e10->setVertices(newVertex[6],newVertex[8]);
      Edge* e11 = new Edge;
      e11->setVertices(newVertex[7],newVertex[8]);

      newEdge[8]  = e8;
      newEdge[9]  = e9;
      newEdge[10] = e10;
      newEdge[11] = e11;

      gridLevel.appendEdge(e8);
      gridLevel.appendEdge(e9);
      gridLevel.appendEdge(e10);
      gridLevel.appendEdge(e11);


      // if volume not refined init theChild field
      for (i = 0; i < 4; i++)
        theChild[i] = 0;
    }
  } // it has no volumen neighbor
  else {
    // create vertex at face center
    double pos[3];
    for (i = 0; i < 3; i++) {
      pos[i] = 0.0;
      for (vertex = 0; vertex < 4; vertex++) {
        pos[i] += (*v[vertex])[i];
      }
      pos[i] *= 0.25;
    }
    Vertex* tmpVt = new Vertex(pos);
    gridLevel.appendVertex(tmpVt);
    newVertex[8] = tmpVt;

    // create and set internal edges
    Edge* e8 = new Edge;
    e8->setVertices(newVertex[4],newVertex[8]);
    Edge* e9 = new Edge;
    e9->setVertices(newVertex[5],newVertex[8]);
    Edge* e10 = new Edge;
    e10->setVertices(newVertex[6],newVertex[8]);
    Edge* e11 = new Edge;
    e11->setVertices(newVertex[7],newVertex[8]);

    newEdge[8]  = e8;
    newEdge[9]  = e9;
    newEdge[10] = e10;
    newEdge[11] = e11;

    gridLevel.appendEdge(e8);
    gridLevel.appendEdge(e9);
    gridLevel.appendEdge(e10);
    gridLevel.appendEdge(e11);

    // if volume not refined init theChild field
    for (i = 0; i < 4; i++)
      theChild[i] = 0;
  }


  // create children
  child = new Element*[4];

  child[0] = new Quadrangle;
  child[1] = new Quadrangle;
  child[2] = new Quadrangle;
  child[3] = new Quadrangle;

  gridLevel.appendQuadrangle(child[0]);
  gridLevel.appendQuadrangle(child[1]);
  gridLevel.appendQuadrangle(child[2]);
  gridLevel.appendQuadrangle(child[3]);

  // set children parent
  child[0]->setParent(this);
  child[1]->setParent(this);
  child[2]->setParent(this);
  child[3]->setParent(this);


  // assemble children from new elements
  for (ch = 0; ch < 4; ch++) {
    // set vertices
    for (vertex = 0; vertex < 4; vertex++) {
      vertexIndex = quaVertexOfQuaChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 4; edge++) {
      edgeIndex = quaEdgeOfQuaChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
  }

  // connect children with volume
  //index = (4 - offset) % 4;
  //child[0]->setNeighbor(0,theChild[index]);
  index = (4 - offset) % 4;
  child[0]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[0]);

  //index = (5 - offset) % 4;
  //child[3]->setNeighbor(0,theChild[index]);
  index = (5 - offset) % 4;
  child[1]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[1]);

  //index = (6 - offset) % 4;
  //child[2]->setNeighbor(0,theChild[index]);
  index = (6 - offset) % 4;
  child[2]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[2]);

  //index = (7 - offset) % 4;
  //child[1]->setNeighbor(0,theChild[3]);
  index = (7 - offset) % 4;
  child[3]->setNeighbor(0,theChild[index]);
  if (theChild[3])
    theChild[3]->setNeighbor(theChildFace[index],child[3]);

  // set tag
  setRefined();
  
  // set no. of children
  chcount = 4;

  // set value
  for (i = 0; i < 4; i++)
    child[i]->setValue(getValue());
}





// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Quadrangle::coarsen(GridLevel& gridLevel)
{
  size_t noOfChild = 4;
  size_t ch;
  int nFace;
  Element* theNeighbor;

  for (ch = 0; ch < noOfChild; ch++)
  {
    // *- set neighbor's pointer
    theNeighbor = child[ch]->getNeighbor(0);
    if (theNeighbor) {
      nFace = theNeighbor->getCommonFace(child[ch]);
      theNeighbor->setNeighbor(nFace,0);
    }
    gridLevel.removeQuadrangle(child[ch]);
  }

  // Delete the children from memory
  for (ch = 0; ch < noOfChild; ch++)
    delete child[ch];
  // Delete the pointers
  delete [] child;
  // Initialize pointer
  child = 0;

  // set tag
  setRegular();
  // set no. of children
  chcount = 0;

}


// Description
//
void
Quadrangle::close(GridLevel& gridLevel)
{
  size_t refPattern = 0;
  size_t refMask    = 1;
  size_t noOfTriElem = 0;
  size_t noOfQuaElem = 0;

  Element* tris[4];
  Element* quads[2];

  Vertex* vertex[8];

  // Edge refinement pattern
  for (size_t i = 0; i < 4; i++)
  {
    vertex[i] = v[i];
    if (e[i]->isRefined()) {
      refPattern = (refPattern | refMask);
      vertex[4+i] = e[i]->getMidPoint();
    }
    else
      vertex[4+i] = 0;
    refMask = (refMask<<1);
  }


  switch(refPattern)
  {
    case 1:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[4]);
    tris[1]->setVertex(1,vertex[3]);
    tris[1]->setVertex(2,vertex[2]);

    tris[2]->setVertex(0,vertex[4]);
    tris[2]->setVertex(1,vertex[2]);
    tris[2]->setVertex(2,vertex[1]);

    noOfTriElem = 3;
    break;
  case 2:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[5]);

    tris[1]->setVertex(0,vertex[3]);
    tris[1]->setVertex(1,vertex[2]);
    tris[1]->setVertex(2,vertex[5]);

    tris[2]->setVertex(0,vertex[0]);
    tris[2]->setVertex(1,vertex[5]);
    tris[2]->setVertex(2,vertex[1]);

    noOfTriElem = 3;
    break;
  case 3:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;
    tris[3] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[4]);
    tris[1]->setVertex(1,vertex[3]);
    tris[1]->setVertex(2,vertex[5]);

    tris[2]->setVertex(0,vertex[4]);
    tris[2]->setVertex(1,vertex[5]);
    tris[2]->setVertex(2,vertex[1]);

    tris[3]->setVertex(0,vertex[5]);
    tris[3]->setVertex(1,vertex[3]);
    tris[3]->setVertex(2,vertex[2]);

    noOfTriElem = 4;
    break;
  case 4:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[6]);

    tris[1]->setVertex(0,vertex[0]);
    tris[1]->setVertex(1,vertex[6]);
    tris[1]->setVertex(2,vertex[1]);

    tris[2]->setVertex(0,vertex[6]);
    tris[2]->setVertex(1,vertex[2]);
    tris[2]->setVertex(2,vertex[1]);

    noOfTriElem = 3;
    break;
  case 5:
    quads[0] = new Quadrangle;
    quads[1] = new Quadrangle;

    quads[0]->setVertex(0,vertex[0]);
    quads[0]->setVertex(1,vertex[4]);
    quads[0]->setVertex(2,vertex[6]);
    quads[0]->setVertex(3,vertex[3]);

    quads[1]->setVertex(0,vertex[4]);
    quads[1]->setVertex(1,vertex[1]);
    quads[1]->setVertex(2,vertex[2]);
    quads[1]->setVertex(3,vertex[6]);

    noOfQuaElem = 2;
    break;
  case 6:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;
    tris[3] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[6]);

    tris[1]->setVertex(0,vertex[0]);
    tris[1]->setVertex(1,vertex[6]);
    tris[1]->setVertex(2,vertex[5]);

    tris[2]->setVertex(0,vertex[6]);
    tris[2]->setVertex(1,vertex[2]);
    tris[2]->setVertex(2,vertex[5]);

    tris[3]->setVertex(0,vertex[0]);
    tris[3]->setVertex(1,vertex[5]);
    tris[3]->setVertex(2,vertex[1]);

    noOfTriElem = 4;
    break;
  case 7:
    quads[0] = new Quadrangle;

    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[4]);
    tris[0]->setVertex(1,vertex[6]);
    tris[0]->setVertex(2,vertex[5]);

    tris[1]->setVertex(0,vertex[5]);
    tris[1]->setVertex(1,vertex[6]);
    tris[1]->setVertex(2,vertex[2]);

    tris[2]->setVertex(0,vertex[4]);
    tris[2]->setVertex(1,vertex[5]);
    tris[2]->setVertex(2,vertex[1]);

    quads[0]->setVertex(0,vertex[0]);
    quads[0]->setVertex(1,vertex[4]);
    quads[0]->setVertex(2,vertex[6]);
    quads[0]->setVertex(3,vertex[3]);

    noOfTriElem = 3;
    noOfQuaElem = 1;
    break;
  case 8:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[1]);

    tris[1]->setVertex(0,vertex[1]);
    tris[1]->setVertex(1,vertex[7]);
    tris[1]->setVertex(2,vertex[2]);

    tris[2]->setVertex(0,vertex[7]);
    tris[2]->setVertex(1,vertex[3]);
    tris[2]->setVertex(2,vertex[2]);

    noOfTriElem = 3;
    break;
  case 9:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;
    tris[3] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[4]);
    tris[1]->setVertex(1,vertex[7]);
    tris[1]->setVertex(2,vertex[2]);

    tris[2]->setVertex(0,vertex[7]);
    tris[2]->setVertex(1,vertex[3]);
    tris[2]->setVertex(2,vertex[2]);

    tris[3]->setVertex(0,vertex[4]);
    tris[3]->setVertex(1,vertex[2]);
    tris[3]->setVertex(2,vertex[1]);

    noOfTriElem = 4;
    break;
  case 10:
    quads[0] = new Quadrangle;
    quads[1] = new Quadrangle;

    quads[0]->setVertex(0,vertex[0]);
    quads[0]->setVertex(1,vertex[1]);
    quads[0]->setVertex(2,vertex[5]);
    quads[0]->setVertex(3,vertex[7]);

    quads[1]->setVertex(0,vertex[7]);
    quads[1]->setVertex(1,vertex[5]);
    quads[1]->setVertex(2,vertex[2]);
    quads[1]->setVertex(3,vertex[3]);

    noOfQuaElem = 2;
    break;
  case 11:
    quads[0] = new Quadrangle;

    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[4]);
    tris[1]->setVertex(1,vertex[7]);
    tris[1]->setVertex(2,vertex[5]);

    tris[2]->setVertex(0,vertex[4]);
    tris[2]->setVertex(1,vertex[5]);
    tris[2]->setVertex(2,vertex[1]);

    quads[0]->setVertex(0,vertex[7]);
    quads[0]->setVertex(1,vertex[5]);
    quads[0]->setVertex(2,vertex[2]);
    quads[0]->setVertex(3,vertex[3]);

    noOfTriElem = 3;
    noOfQuaElem = 1;
    break;
  case 12:
    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;
    tris[3] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[1]);

    tris[1]->setVertex(0,vertex[1]);
    tris[1]->setVertex(1,vertex[7]);
    tris[1]->setVertex(2,vertex[6]);

    tris[2]->setVertex(0,vertex[7]);
    tris[2]->setVertex(1,vertex[3]);
    tris[2]->setVertex(2,vertex[6]);

    tris[3]->setVertex(0,vertex[1]);
    tris[3]->setVertex(1,vertex[6]);
    tris[3]->setVertex(2,vertex[2]);

    noOfTriElem = 4;
    break;
  case 13:
    quads[0] = new Quadrangle;

    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[4]);
    tris[1]->setVertex(1,vertex[7]);
    tris[1]->setVertex(2,vertex[6]);

    tris[2]->setVertex(0,vertex[7]);
    tris[2]->setVertex(1,vertex[3]);
    tris[2]->setVertex(2,vertex[6]);

    quads[0]->setVertex(0,vertex[4]);
    quads[0]->setVertex(1,vertex[1]);
    quads[0]->setVertex(2,vertex[2]);
    quads[0]->setVertex(3,vertex[6]);

    noOfTriElem = 3;
    noOfQuaElem = 1;
    break;
  case 14:
    quads[0] = new Quadrangle;

    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;

    tris[0]->setVertex(0,vertex[7]);
    tris[0]->setVertex(1,vertex[3]);
    tris[0]->setVertex(2,vertex[6]);

    tris[1]->setVertex(0,vertex[7]);
    tris[1]->setVertex(1,vertex[6]);
    tris[1]->setVertex(2,vertex[5]);

    tris[2]->setVertex(0,vertex[5]);
    tris[2]->setVertex(1,vertex[6]);
    tris[2]->setVertex(2,vertex[2]);

    quads[0]->setVertex(0,vertex[0]);
    quads[0]->setVertex(1,vertex[1]);
    quads[0]->setVertex(2,vertex[5]);
    quads[0]->setVertex(3,vertex[7]);

    noOfTriElem = 2;
    noOfQuaElem = 1;
    break;
  case 15:
    quads[0] = new Quadrangle;

    tris[0] = new Triangle;
    tris[1] = new Triangle;
    tris[2] = new Triangle;
    tris[3] = new Triangle;

    tris[0]->setVertex(0,vertex[0]);
    tris[0]->setVertex(1,vertex[7]);
    tris[0]->setVertex(2,vertex[4]);

    tris[1]->setVertex(0,vertex[7]);
    tris[1]->setVertex(1,vertex[3]);
    tris[1]->setVertex(2,vertex[6]);

    tris[2]->setVertex(0,vertex[6]);
    tris[2]->setVertex(1,vertex[2]);
    tris[2]->setVertex(2,vertex[5]);

    tris[3]->setVertex(0,vertex[5]);
    tris[3]->setVertex(1,vertex[1]);
    tris[3]->setVertex(2,vertex[4]);

    quads[0]->setVertex(0,vertex[4]);
    quads[0]->setVertex(1,vertex[5]);
    quads[0]->setVertex(2,vertex[6]);
    quads[0]->setVertex(3,vertex[7]);

    noOfTriElem = 4;
    noOfQuaElem = 1;
    break;
    default:
      noOfTriElem = 0;
      noOfQuaElem = 0;
      break;
  }

  // Connect hierarchy
  size_t noOfChildren = noOfTriElem + noOfQuaElem;
  child = new Element*[noOfChildren];

  // Append elements in next grid level
  size_t count = 0;
  for (size_t nn = 0; nn < noOfTriElem; nn++)
  {
    tris[nn]->setClosure();
    tris[nn]->setValue(getValue());
    tris[nn]->setParent(this);
    child[count++] = tris[nn];
    gridLevel.appendTriangle(tris[nn]);
  }
  for (size_t qq = 0; qq < noOfQuaElem; qq++)
  {
    quads[qq]->setClosure();
    quads[qq]->setValue(getValue());
    quads[qq]->setParent(this);
    child[count++] = quads[qq];
    gridLevel.appendQuadrangle(quads[qq]);
  }
  
  // set no. of children
  chcount = static_cast<unsigned char>(noOfTriElem + noOfQuaElem);
  // set irregular mark
  setIrregular();
  // set neighbors and create new edges
  connect(gridLevel);
}





// Description
//  computes the edge refinement pattern
int
Quadrangle::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;
  for (int i = 0; i < 4; i++)
  {
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
//	         virtual functions
//
// ____________________________________________________________
//*************************************************************

//=============================================================================
//  TOPOLOGY
//=============================================================================

// Description
//
bool
Quadrangle::queryEdgeRefinement() const
{
  for (int i = 0; i < 4; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Quadrangle::queryChildrenEdgeRefinement() const
{
  for (int i = 0; i < 4; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  return false;
}


// Description
// function for 3D elements
int
Quadrangle::getCommonFace(Element* theNeighbor) const
{
  if (n == theNeighbor)
    return 0;
  else
    return -1;
}



// Description
// common elements for quads as faces of hexa or pyramids
void
Quadrangle::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
{
  if (face == 0) {
    // order determines the orientation
    // towards the volume
    theVertex[0] = v[0];
    theVertex[1] = v[3];
    theVertex[2] = v[2];
    theVertex[3] = v[1];
    // set vertex at barycenter
    theVertex[4] = child[0]->getVertex(2);

    // get Edges at face (internal edges)
    theEdge[0] = child[1]->getEdge(3);    // e8
    theEdge[1] = child[0]->getEdge(2);    // e11
    theEdge[2] = child[3]->getEdge(1);    // e10
    theEdge[3] = child[2]->getEdge(0);    // e9

    theChild[0] = child[0];
    theChild[1] = child[3];
    theChild[2] = child[2];
    theChild[3] = child[1];

    // connect with volume at face 0
    theChildFace[0] = 0;
    theChildFace[1] = 0;
    theChildFace[2] = 0;
    theChildFace[3] = 0;

  }
}


// Description
//
void
Quadrangle::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
void
Quadrangle::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 4; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Quadrangle::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Quadrangle::getVertices(Vertex** vt) const
{
  for (int i = 0; i < 4; i++)
    vt[i] = v[i];
}


// Description
//
void
Quadrangle::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
void
Quadrangle::setEdges(Edge** edg)
{
  for (int i = 0; i < 4; i++)
    e[i] = edg[i];
}



// Description
//
Edge*
Quadrangle::getEdge(int i) const
{
  return e[i];
}


// Description
//
void
Quadrangle::getEdges(Edge** edg) const
{
  for (int i = 0; i < 4; i++)
    edg[i] = e[i];
}


// Description
//
void
Quadrangle::setNeighbor(int face,Element* theNeighbor)
{
  if (face == 0)
    n = theNeighbor;
}


// Description
//
Element*
Quadrangle::getNeighbor(int face) const
{
  if (face == 0)
    return n;
  else
    return 0;
}

//=============================================================================
//  TOPOLOGY
//=============================================================================

// Description
//
int
Quadrangle::type() const
{
  return GRD_QUADRANGLE;
}

// Description
//
int
Quadrangle::noOfVerticesAtFace(int fc) const
{
  return 4;
}


// Description
//
int
Quadrangle::topologicalVertexAtFace(int fc,int vt) const
{
  if (fc == 0) return vt;
  else return -1;
}


// Description
//
Vertex*
Quadrangle::getVertexAtFace(int fc,int vt) const
{
  if (fc == 0) return v[vt];
  else return 0;
  
}

// Description
// 
int
Quadrangle::noOfEdgesAtFace(int fc) const
{
  if (fc == 0) return 4;
  else return -1;  
}

// Description
//
Vertex*
Quadrangle::getVertexAtEdge(int ed,int vt) const
{
  return v[quaVertexOfEdge[ed][vt]];
}


// Description
//
Edge*
Quadrangle::getEdgeAtFace(int fc,int ed) const
{
  if (fc == 0) return e[ed];
  else return 0;
}



// Description
//  Auxiliary function for Quadrangle
bool compareVertices(Vertex** vt,Vertex** vv,int nov)
{
  int count = 0;
  for (int i = 0; i < nov; i++)
  {
    for (int j = 0; j < nov; j++)
    {
      if (vt[i] == vv[j])
      {
        count++;
        break;
      }
    }
  }
  if (count == nov) return true;
  else return false; 
  
}


// Description
//  Connect new irregular elements with
//  the neighbors and create new edges
void
Quadrangle::connect(GridLevel& gridLevel)
{
  // Create List of objects to be connected
  std::list<Element*> nelt;

  // Collect all neighbors from volume
  if (n)
  {
    if (n->getNoOfChildren() > 0)
    {
      for (int ch = 0; ch < n->getNoOfChildren(); ch++)
        nelt.push_back(n->getChild(ch));
    }
    else
    {
      if (!n->isIrregular())
        nelt.push_back(n);
    }
  }  
  // Connect elements in the list
  typedef std::list<Element*>::iterator EI;
  for (int ch = 0; ch < static_cast<int>(chcount); ch++)
  {
    // Quadrangles have only one neighbor: the volume
    Element* neighbor = child[ch]->getNeighbor(0);
    if (neighbor)
      continue;
    else
    {
      Vertex* vt[4];
      Vertex* vv[4];
      int noOfVt = child[ch]->getNoOfVertices();
      for (size_t nn = 0; nn < (size_t)noOfVt; nn++) vt[nn] = child[ch]->getVertex(nn);
      
      // process all neighbors
      for (EI q = nelt.begin(); q != nelt.end(); ++q)
      {
        for (int nn = 0; nn < (*q)->getNoOfFaces(); nn++)
        {
          int noOfVv = (*q)->noOfVerticesAtFace(nn);
          if (noOfVv != noOfVt)
            continue;
            
          for (size_t i = 0; i < (size_t)noOfVv; i++) vv[i] = (*q)->getVertexAtFace(nn,i);
          
          if (compareVertices(vt,vv,noOfVt))
          {
            // neighbour found
            child[ch]->setNeighbor(0,*q);
            (*q)->setNeighbor(nn,child[ch]);
          }
        }
      }
    }
  }
  
  // Check edges
  // Create List of objects to be connected
  typedef std::list<Edge*>::iterator EGI;
  std::list<Edge*> edglt;
  // Get edges from children
  for (int ch = 0; ch < static_cast<int>(chcount); ch++)
  {
    for (int nn = 0; nn < child[ch]->getNoOfEdges(); nn++)
    {
      Edge* tmp = child[ch]->getEdge(nn);
      if (tmp) edglt.push_back(tmp);
    }
  }
  // get edges from neighbors
  for (EI q = nelt.begin(); q != nelt.end(); ++q)
  {
    for (int nn = 0; nn < (*q)->getNoOfEdges(); nn++)
    {
      Edge* tmp = (*q)->getEdge(nn);
      if (tmp) edglt.push_back(tmp);
    } 
  }
  
  // Get the edges of the element
  for (int nn = 0; nn < 4; nn++)
  {
    if (!e[nn]->isRefined())
      edglt.push_back(e[nn]);
    else
    {
      edglt.push_back(e[nn]->getChild(0));
      edglt.push_back(e[nn]->getChild(1));
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
  edglt.clear();
}

//// Description
////
//void
//Quadrangle::getGlobalCoords(const double* lc,double* coord)
//{
//  double u = lc[0];
//  double v = lc[1];
//
//  coord[0] = (((*v[0])[0] - (*v[1])[0] + (*v[2])[0] - (*v[3])[0]) * v - (*v[0])[0] + (*v[1])[0]) * u + ((*v[3])[0] - (*v[0])[0]) * v + (*v[0])[0];
//  coord[1] = (((*v[0])[1] - (*v[1])[1] + (*v[2])[1] - (*v[3])[1]) * v - (*v[0])[1] + (*v[1])[1]) * u + ((*v[3])[1] - (*v[0])[1]) * v + (*v[0])[1];
//  coord[2] = (((*v[0])[2] - (*v[1])[2] + (*v[2])[2] - (*v[3])[2]) * v - (*v[0])[2] + (*v[1])[2]) * u + ((*v[3])[2] - (*v[0])[2]) * v + (*v[0])[2];
//
//}
//
//// Description
////
//void
//Quadrangle::getJacobian(Jacobian& Jac)
//{
//
//  Jac[0][0] = ((*v[0])[0] - (*v[1])[0] + (*v[2])[0] - (*v[3])[0])*0.5 - (*v[0])[0] + (*v[1])[0];
//  Jac[0][1] = ((*v[0])[0] - (*v[1])[0] + (*v[2])[0] - (*v[3])[0])*0.5 - (*v[0])[0] + (*v[3])[0];
//  Jac[0][2] = 0.0;
//
//  Jac[1][0] = ((*v[0])[1] - (*v[1])[1] + (*v[2])[1] - (*v[3])[1])*0.5 - (*v[0])[1] + (*v[1])[1];
//  Jac[1][1] = ((*v[0])[1] - (*v[1])[1] + (*v[2])[1] - (*v[3])[1])*0.5 - (*v[0])[1] + (*v[3])[1];
//  Jac[1][2] = 0.0;
//
//  Jac[2][0] = ((*v[0])[2] - (*v[1])[2] + (*v[2])[2] - (*v[3])[2])*0.5 - (*v[0])[2] + (*v[1])[2];
//  Jac[2][1] = ((*v[0])[2] - (*v[1])[2] + (*v[2])[2] - (*v[3])[2])*0.5 - (*v[0])[2] + (*v[3])[2];
//  Jac[2][2] = 0.0;
//}


} // namespace
