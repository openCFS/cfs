/***************************************************************************
    File        : Quadrangle.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 8 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Quadrangle.h"
#include "GridLevel.h"
#include "ConformingClosure.h"

namespace grd {

// (static obj) initialize memory management here
Pool<Quadrangle> Quadrangle::pool;

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
}



// Description
// if all children are connected to volume elements
// only for 3D elements
bool
Quadrangle::queryChildrenMarks()
{
  int i;
  int counter = 0;

  if (!isRefined())
    return false;
  for (i = 0; i < 4; i++)
  {
    if (child[i]->isMarkedForCoarsening())
      counter++;
  }
  // set back children makrs
  if (counter == 0)
    return false;
  else if (counter > 0 && counter < 4)
  {
    for (i = 0; i < 4; i++)
    {
      if (child[i]->isMarkedForCoarsening())
        child[i]->setRegular();
    }
    return false;
  }
  else if (counter == 4)
  {
    // can't refine if an child's edge is refined
    bool flag = false;
    for (i = 0; i < 4; i++)
    {
      if (child[i]->queryEdgeRefinement()) {
        flag = true;
        break;
      }
    }
    if (flag)
    {
      for (i = 0; i < 4; i++)
      {
        if (child[i]->isMarkedForCoarsening())
          child[i]->setRegular();
      }
      return false;
    }
    else
      return true;
  }
  else
    return false;

}


// Description
//
// refine triangle regular into 4 child triangles
// insert new elements in the lists.
// input: tesselation - contains pointer to the element lists
// output: bool true or false
// -* 1. process parent vertices and get new vertices
// -* 2. process parent edges and get new vertices + new edges
// -* 3. process parent faces and get new edges + neighbors
// *- 4. assemble children from new elements
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
        cerr << "ERROR at Quad: can't find common vertex from Quad face\n";
        exit(1);
      }

      index = (4 - offset) % 4;
      newEdge[8] = theEdge[index];

      index = (5 - offset) % 4;
      newEdge[11] = theEdge[index];

      index = (6 - offset) % 4;
      newEdge[10] = theEdge[index];

      index = (7 - offset) % 4;
      newEdge[9] = theEdge[index];

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
  index = (4 - offset) % 4;
  child[0]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[0]);

  index = (5 - offset) % 4;
  child[3]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[3]);


  index = (6 - offset) % 4;
  child[2]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[2]);

  index = (7 - offset) % 4;
  child[1]->setNeighbor(0,theChild[3]);
  if (theChild[3])
    theChild[3]->setNeighbor(theChildFace[3],child[1]);

  // set tag
  setRefined();

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
  int ch;
  int nFace;
  Element* theNeighbor;

  for (ch = 0; ch < 4; ch++) {
    // *- set neighbor's pointer
    theNeighbor = child[ch]->getNeighbor(0);
    if (theNeighbor) {
      nFace = theNeighbor->getCommonFace(child[ch]);
      theNeighbor->setNeighbor(nFace,0);
    }
    gridLevel.removeQuadrangle(child[ch]);
  }

  for (ch = 0; ch < 4; ch++)
    delete child[ch];
  child = 0;

  // set tag
  setRegular();

}



// Description:
//
void
Quadrangle::close(int refPattern,ConformingClosure& closure)
{
  //int refPattern;
  int noOfElem;
  int id1,id2;

  ConformingClosure::TRI theElement = closure.beginTriangle();
  Vertex* vertex[8];

  // init element buffer
  closure.init();
  for (int i = 0; i < 8; i++)
    vertex[i] = 0;

  //refPattern = getRefinementPattern(vertex);

  switch(refPattern) {
    case 1:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[3]);
    theElement[1]->setVertex(2,vertex[2]);

    theElement[2]->setVertex(0,vertex[4]);
    theElement[2]->setVertex(1,vertex[2]);
    theElement[2]->setVertex(2,vertex[1]);

    noOfElem = 3;
    break;
  case 2:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[5]);

    theElement[1]->setVertex(0,vertex[3]);
    theElement[1]->setVertex(1,vertex[2]);
    theElement[1]->setVertex(2,vertex[5]);

    theElement[2]->setVertex(0,vertex[0]);
    theElement[2]->setVertex(1,vertex[5]);
    theElement[2]->setVertex(2,vertex[1]);

    noOfElem = 3;
    break;
  case 3:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[3]);
    theElement[1]->setVertex(2,vertex[5]);

    theElement[2]->setVertex(0,vertex[4]);
    theElement[2]->setVertex(1,vertex[5]);
    theElement[2]->setVertex(2,vertex[1]);

    theElement[3]->setVertex(0,vertex[5]);
    theElement[3]->setVertex(1,vertex[3]);
    theElement[3]->setVertex(2,vertex[2]);

    noOfElem = 4;
    break;
  case 4:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[6]);

    theElement[1]->setVertex(0,vertex[0]);
    theElement[1]->setVertex(1,vertex[6]);
    theElement[1]->setVertex(2,vertex[1]);

    theElement[2]->setVertex(0,vertex[6]);
    theElement[2]->setVertex(1,vertex[2]);
    theElement[2]->setVertex(2,vertex[1]);

    noOfElem = 3;
    break;
  case 5:
    id1 = vertex[0]->getId();
    id2 = vertex[3]->getId();
    if (id1 > id2) {
      theElement[0]->setVertex(0,vertex[0]);
      theElement[0]->setVertex(1,vertex[3]);
      theElement[0]->setVertex(2,vertex[4]);

      theElement[1]->setVertex(0,vertex[4]);
      theElement[1]->setVertex(1,vertex[3]);
      theElement[1]->setVertex(2,vertex[6]);

      theElement[2]->setVertex(0,vertex[4]);
      theElement[2]->setVertex(1,vertex[6]);
      theElement[2]->setVertex(2,vertex[2]);

      theElement[3]->setVertex(0,vertex[4]);
      theElement[3]->setVertex(1,vertex[2]);
      theElement[3]->setVertex(2,vertex[1]);
    }
    else {
      theElement[0]->setVertex(0,vertex[0]);
      theElement[0]->setVertex(1,vertex[3]);
      theElement[0]->setVertex(2,vertex[6]);

      theElement[1]->setVertex(0,vertex[0]);
      theElement[1]->setVertex(1,vertex[6]);
      theElement[1]->setVertex(2,vertex[4]);

      theElement[2]->setVertex(0,vertex[4]);
      theElement[2]->setVertex(1,vertex[6]);
      theElement[2]->setVertex(2,vertex[1]);

      theElement[3]->setVertex(0,vertex[1]);
      theElement[3]->setVertex(1,vertex[6]);
      theElement[3]->setVertex(2,vertex[2]);
    }

    noOfElem = 4;
    break;
  case 6:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[6]);

    theElement[1]->setVertex(0,vertex[0]);
    theElement[1]->setVertex(1,vertex[6]);
    theElement[1]->setVertex(2,vertex[5]);

    theElement[2]->setVertex(0,vertex[6]);
    theElement[2]->setVertex(1,vertex[2]);
    theElement[2]->setVertex(2,vertex[5]);

    theElement[3]->setVertex(0,vertex[0]);
    theElement[3]->setVertex(1,vertex[5]);
    theElement[3]->setVertex(2,vertex[1]);

    noOfElem = 4;
    break;
  case 7:
    theElement[0]->setVertex(0,vertex[4]);
    theElement[0]->setVertex(1,vertex[6]);
    theElement[0]->setVertex(2,vertex[5]);

    theElement[1]->setVertex(0,vertex[5]);
    theElement[1]->setVertex(1,vertex[6]);
    theElement[1]->setVertex(2,vertex[2]);

    theElement[2]->setVertex(0,vertex[4]);
    theElement[2]->setVertex(1,vertex[5]);
    theElement[2]->setVertex(2,vertex[1]);

    id1 = vertex[0]->getId();
    id2 = vertex[3]->getId();

    if (id1 > id2) {
      theElement[3]->setVertex(0,vertex[0]);
      theElement[3]->setVertex(1,vertex[3]);
      theElement[3]->setVertex(2,vertex[6]);

      theElement[4]->setVertex(0,vertex[0]);
      theElement[4]->setVertex(1,vertex[6]);
      theElement[4]->setVertex(2,vertex[4]);
    }
    else {
      theElement[3]->setVertex(0,vertex[0]);
      theElement[3]->setVertex(1,vertex[3]);
      theElement[3]->setVertex(2,vertex[4]);

      theElement[4]->setVertex(0,vertex[4]);
      theElement[4]->setVertex(1,vertex[3]);
      theElement[4]->setVertex(2,vertex[6]);
    }

    noOfElem = 5;
    break;
  case 8:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[1]);

    theElement[1]->setVertex(0,vertex[1]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[2]);

    theElement[2]->setVertex(0,vertex[7]);
    theElement[2]->setVertex(1,vertex[3]);
    theElement[2]->setVertex(2,vertex[2]);

    noOfElem = 3;
    break;
  case 9:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[2]);

    theElement[2]->setVertex(0,vertex[7]);
    theElement[2]->setVertex(1,vertex[3]);
    theElement[2]->setVertex(2,vertex[2]);

    theElement[3]->setVertex(0,vertex[4]);
    theElement[3]->setVertex(1,vertex[2]);
    theElement[3]->setVertex(2,vertex[1]);

    noOfElem = 4;
    break;
  case 10:
    id1 = vertex[0]->getId();
    id2 = vertex[1]->getId();

    if (id1 > id2) {
      theElement[0]->setVertex(0,vertex[0]);
      theElement[0]->setVertex(1,vertex[7]);
      theElement[0]->setVertex(2,vertex[5]);

      theElement[1]->setVertex(0,vertex[0]);
      theElement[1]->setVertex(1,vertex[5]);
      theElement[1]->setVertex(2,vertex[1]);

      theElement[2]->setVertex(0,vertex[5]);
      theElement[2]->setVertex(1,vertex[7]);
      theElement[2]->setVertex(2,vertex[3]);

      theElement[3]->setVertex(0,vertex[5]);
      theElement[3]->setVertex(1,vertex[3]);
      theElement[3]->setVertex(2,vertex[2]);
    }
    else {
      theElement[0]->setVertex(0,vertex[0]);
      theElement[0]->setVertex(1,vertex[7]);
      theElement[0]->setVertex(2,vertex[1]);

      theElement[1]->setVertex(0,vertex[1]);
      theElement[1]->setVertex(1,vertex[7]);
      theElement[1]->setVertex(2,vertex[5]);

      theElement[2]->setVertex(0,vertex[5]);
      theElement[2]->setVertex(1,vertex[7]);
      theElement[2]->setVertex(2,vertex[2]);

      theElement[3]->setVertex(0,vertex[7]);
      theElement[3]->setVertex(1,vertex[3]);
      theElement[3]->setVertex(2,vertex[2]);
    }

    noOfElem = 4;
    break;
  case 11:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[5]);

    theElement[2]->setVertex(0,vertex[4]);
    theElement[2]->setVertex(1,vertex[5]);
    theElement[2]->setVertex(2,vertex[1]);

    id1 = vertex[2]->getId();
    id2 = vertex[3]->getId();

    if (id1 > id2) {
      theElement[3]->setVertex(0,vertex[7]);
      theElement[3]->setVertex(1,vertex[2]);
      theElement[3]->setVertex(2,vertex[5]);

      theElement[4]->setVertex(0,vertex[7]);
      theElement[4]->setVertex(1,vertex[3]);
      theElement[4]->setVertex(2,vertex[2]);

    }
    else {
      theElement[3]->setVertex(0,vertex[5]);
      theElement[3]->setVertex(1,vertex[7]);
      theElement[3]->setVertex(2,vertex[3]);

      theElement[4]->setVertex(0,vertex[5]);
      theElement[4]->setVertex(1,vertex[3]);
      theElement[4]->setVertex(2,vertex[2]);
    }

    noOfElem = 5;
    break;
  case 12:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[1]);

    theElement[1]->setVertex(0,vertex[1]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[6]);

    theElement[2]->setVertex(0,vertex[7]);
    theElement[2]->setVertex(1,vertex[3]);
    theElement[2]->setVertex(2,vertex[6]);

    theElement[3]->setVertex(0,vertex[1]);
    theElement[3]->setVertex(1,vertex[6]);
    theElement[3]->setVertex(2,vertex[2]);

    noOfElem = 4;
    break;
  case 13:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[6]);

    theElement[2]->setVertex(0,vertex[7]);
    theElement[2]->setVertex(1,vertex[3]);
    theElement[2]->setVertex(2,vertex[6]);

    id1 = vertex[1]->getId();
    id2 = vertex[2]->getId();

    if (id1 > id2) {
      theElement[3]->setVertex(0,vertex[4]);
      theElement[3]->setVertex(1,vertex[6]);
      theElement[3]->setVertex(2,vertex[1]);

      theElement[4]->setVertex(0,vertex[1]);
      theElement[4]->setVertex(1,vertex[6]);
      theElement[4]->setVertex(2,vertex[2]);

    }
    else {
      theElement[3]->setVertex(0,vertex[4]);
      theElement[3]->setVertex(1,vertex[6]);
      theElement[3]->setVertex(2,vertex[2]);

      theElement[4]->setVertex(0,vertex[4]);
      theElement[4]->setVertex(1,vertex[2]);
      theElement[4]->setVertex(2,vertex[1]);
    }

    noOfElem = 5;
    break;
  case 14:
    theElement[0]->setVertex(0,vertex[7]);
    theElement[0]->setVertex(1,vertex[3]);
    theElement[0]->setVertex(2,vertex[6]);

    theElement[1]->setVertex(0,vertex[7]);
    theElement[1]->setVertex(1,vertex[6]);
    theElement[1]->setVertex(2,vertex[5]);

    theElement[2]->setVertex(0,vertex[5]);
    theElement[2]->setVertex(1,vertex[6]);
    theElement[2]->setVertex(2,vertex[2]);

    id1 = vertex[0]->getId();
    id2 = vertex[1]->getId();

    if (id1 > id2) {
      theElement[3]->setVertex(0,vertex[0]);
      theElement[3]->setVertex(1,vertex[7]);
      theElement[3]->setVertex(2,vertex[5]);

      theElement[4]->setVertex(0,vertex[0]);
      theElement[4]->setVertex(1,vertex[5]);
      theElement[4]->setVertex(2,vertex[1]);
    }
    else {
      theElement[3]->setVertex(0,vertex[0]);
      theElement[3]->setVertex(1,vertex[7]);
      theElement[3]->setVertex(2,vertex[1]);

      theElement[4]->setVertex(0,vertex[1]);
      theElement[4]->setVertex(1,vertex[7]);
      theElement[4]->setVertex(2,vertex[5]);
    }

    noOfElem = 5;
    break;
  case 15:
    theElement[0]->setVertex(0,vertex[0]);
    theElement[0]->setVertex(1,vertex[7]);
    theElement[0]->setVertex(2,vertex[4]);

    theElement[1]->setVertex(0,vertex[4]);
    theElement[1]->setVertex(1,vertex[7]);
    theElement[1]->setVertex(2,vertex[6]);

    theElement[2]->setVertex(0,vertex[7]);
    theElement[2]->setVertex(1,vertex[3]);
    theElement[2]->setVertex(2,vertex[6]);

    theElement[3]->setVertex(0,vertex[4]);
    theElement[3]->setVertex(1,vertex[6]);
    theElement[3]->setVertex(2,vertex[5]);

    theElement[4]->setVertex(0,vertex[4]);
    theElement[4]->setVertex(1,vertex[5]);
    theElement[4]->setVertex(2,vertex[1]);

    theElement[5]->setVertex(0,vertex[5]);
    theElement[5]->setVertex(1,vertex[6]);
    theElement[5]->setVertex(2,vertex[2]);

    noOfElem = 6;
    break;
    default:
      noOfElem = 0;
      break;
  }

  closure.setNoOfTriangles(noOfElem);

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
void
Quadrangle::close(ConformingClosure& closure)
{
  int refPattern = getRefinementPattern();
  close(refPattern,closure);
}

// Description
//
bool
Quadrangle::queryChildrenRefinement()
{
  if(isRefined()) {
    for (int i = 0; i < 4; i++)
      if (child[i]->isRefined())
        return true;
  }

  return false;
}


// Description
//
bool
Quadrangle::queryEdgeRefinement()
{
  for (int i = 0; i < 4; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Quadrangle::queryChildrenEdgeRefinement()
{
  for (int i = 0; i < 4; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  return false;
}



// include element topology information
#include "QuadrangleTopology.h"



// Description
// function for 3D elements
int
Quadrangle::getCommonFace(Element* theNeighbor)
{
  if (n == theNeighbor)
    return 0;
  else
    return -1;
}



// Description
// common elements for quads as faces of hexa or pyramids
void
Quadrangle::getRefElements(int face,Vertex** theVertex,
                           Edge** theEdge,Element** theChild,
                           int* theChildFace)
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
Quadrangle::setVertices(Vertex** vt)
{
  for (int i = 0; i < 4; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Quadrangle::getVertex(int i)
{
  return v[i];
}


// Description
//
void
Quadrangle::getVertices(Vertex** vt)
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
Quadrangle::getEdge(int i)
{
  return e[i];
}


// Description
//
void
Quadrangle::getEdges(Edge** edg)
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
Quadrangle::getNeighbor(int face)
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
Quadrangle::type()
{
  return GRD_QUADRANGLE;
}

// Description
//
int
Quadrangle::noOfVertexAtFace(int fc)
{
  return 2;
}

// Description
//
Vertex*
Quadrangle::getVertexAtFace(int fc,int vt)
{
  return v[quaVertexOfFace[fc][vt]];
}

// Description
//
Vertex*
Quadrangle::getVertexAtEdge(int ed,int vt)
{
  return v[quaVertexOfEdge[ed][vt]];
}


// Description
//
Edge*
Quadrangle::getEdgeAtFace(int fc,int ed)
{
  return e[quaEdgeOfFace[fc][ed]];
}


//=============================================================================
//  GEOMETRY
//=============================================================================

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
