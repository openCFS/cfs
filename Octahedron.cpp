/***************************************************************************
    File        : Octahedron.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Octahedron.h"
#include "Tetrahedron.h"
#include "GridLevel.h"
#include "ConformingClosure.h"

namespace grd {

// initialize memory management here
Pool<Octahedron> Octahedron::pool;


// Description:
//
Octahedron::Octahedron( )
{
  int i;
  for (i = 0; i < 7; i++)
    v[i] = 0;

  for (i = 0; i < 12; i++)
    e[i] = 0;

  for (i = 0; i < 8; i++)
    n[i] = 0;
}


// Description
//
void
Octahedron::getTetras(vector<Element*>& tetra)
{
  // set tetra 0
  tetra[0]->setVertex(0,v[0]);
  tetra[0]->setVertex(1,v[2]);
  tetra[0]->setVertex(2,v[1]);
  tetra[0]->setVertex(3,v[6]);

  // set tetra 1
  tetra[1]->setVertex(0,v[0]);
  tetra[1]->setVertex(1,v[3]);
  tetra[1]->setVertex(2,v[2]);
  tetra[1]->setVertex(3,v[6]);

  // set tetra 2
  tetra[2]->setVertex(0,v[0]);
  tetra[2]->setVertex(1,v[4]);
  tetra[2]->setVertex(2,v[3]);
  tetra[2]->setVertex(3,v[6]);

  // set tetra 3
  tetra[3]->setVertex(0,v[0]);
  tetra[3]->setVertex(1,v[1]);
  tetra[3]->setVertex(2,v[4]);
  tetra[3]->setVertex(3,v[6]);


  // set tetra 4
  tetra[4]->setVertex(0,v[5]);
  tetra[4]->setVertex(1,v[1]);
  tetra[4]->setVertex(2,v[2]);
  tetra[4]->setVertex(3,v[6]);

  // set tetra 5
  tetra[5]->setVertex(0,v[5]);
  tetra[5]->setVertex(1,v[2]);
  tetra[5]->setVertex(2,v[3]);
  tetra[5]->setVertex(3,v[6]);

  // set tetra 6
  tetra[6]->setVertex(0,v[5]);
  tetra[6]->setVertex(1,v[3]);
  tetra[6]->setVertex(2,v[4]);
  tetra[6]->setVertex(3,v[6]);

  // set tetra 7
  tetra[7]->setVertex(0,v[5]);
  tetra[7]->setVertex(1,v[4]);
  tetra[7]->setVertex(2,v[1]);
  tetra[7]->setVertex(3,v[6]);
}



// Description
// if all children marked for coarsening
bool
Octahedron::queryChildrenMarks()
{
  int i;
  int counter = 0;

  if (!isRefined())
    return false;
  
  for (i = 0; i < 14; i++)
  {
    if (child[i]->isMarkedForCoarsening())
      counter++;
  }
  // set back children makrs
  if (counter == 0)
    return false;
  if (counter > 0 && counter < 14)
  {
    for (i = 0; i < 14; i++)
    {
      if (child[i]->isMarkedForCoarsening())
        child[i]->setRegular();
    }
    return false;
  }
  else if (counter == 14)
  {
    // an edge of a child is refined?
    bool flag = false;
    for (i = 0; i < 14; i++)
    {
      if (child[i]->queryEdgeRefinement())
      {
        flag = true;
        break;
      }
    }
    if (flag)
    {
      for (i = 0; i < 14; i++)
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


// Description:
//
// refine tetrahedron regular into 4 tetrahedra and one octahedra
// insert new elements in the lists.
// input: tesselation - contains pointer to the element lists
// output: bool true or false
// -* 1. process parent vertices and get new vertices
// -* 2. process parent edges and get new vertices + new edges
// -* 3. process parent faces and get new edges + neighbors
// -* 4. create interior edges from face to barycenter
// -* 5. assemble children from new elements
void
Octahedron::refine(GridLevel& gridLevel)
{
  int ch;
  int i,j,k;
  int vertCount,edgCount;
  int vertex,edge,face,neighbor;
  int vertexIndex,edgeIndex,neighIndex,childIndex;
  int nVertex,nFace;
  int neighFace[32],theChildFace[4];

  Vertex*  newVertex[19];
  Edge*    newEdge[60];
  Element* newNeighbor[46];
  Vertex*  theVertex[3];
  Edge*    theEdge[3];
  Element* theChild[4];


  // -* init counters
  vertCount = edgCount = 0;

  // -* 1. process parent vertices
  for (vertex = 0; vertex < 6; vertex++)
    newVertex[vertCount++] = v[vertex];

  // -* 2. process parent edges
  for (edge = 0; edge < 12; edge++) {
    if (e[edge]->isRefined()) {
      e[edge]->getChildren(newVertex[octVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[octVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      newVertex[vertCount] = e[edge]->getMidPoint();
      gridLevel.appendVertex(newVertex[vertCount]);
      gridLevel.appendEdge(newEdge[edgCount]);
      gridLevel.appendEdge(newEdge[edgCount+1]);

      edgCount += 2;
      vertCount++;
    }
  }

  // add barycenter
  newVertex[vertCount] = v[6];


  // -* 3. process parent faces

  for (face = 0; face < 8; face++) {
    // if neighbor regular refined
    if (n[face]) {
      if (n[face]->isRefined()) {
//        int type = n[face]->type();
//        if (type == GRD_TRIANGLE)
//          int dummy = 0;
        // get neighbor common face
        nFace = n[face]->getCommonFace(this);
        n[face]->getRefElements(nFace,theVertex,theEdge,theChild,theChildFace);

        // for all vertices of neighbor's face
        int index = 2;
        vertex    = 2;
        for (nVertex = 0; nVertex < 3; nVertex++) {
          // find common vertex of parent and neighbor
          vertex = ((index--)%3);
          vertexIndex = octVertexOfFace[face][vertex];
          if (index < 0) index = 2;
          while(v[vertexIndex] != theVertex[nVertex]) {
            vertex      = ((index--)%3);
            vertexIndex = octVertexOfFace[face][vertex];
            if (index < 0) index = 2;
          }
          // set child edge
          edgeIndex = octEdgeFromParentFace[face][vertex];
          newEdge[edgeIndex] = theEdge[nVertex];
          // set child neighbor
          neighIndex = octNeighborOfParentFace[face][vertex];
          newNeighbor[neighIndex] = theChild[nVertex];
          neighFace[neighIndex]   = theChildFace[nVertex];
        }
        //         set last neighbor from face middle
        neighIndex = octNeighborOfParentFace[face][3];
        newNeighbor[neighIndex] = theChild[3];
        neighFace[neighIndex]   = theChildFace[3];
      }
      else {
        // neighbor is not refined
        // create new edges on face
        // set children neighbors at face to NULL
        int v1,v2;
        for (edge = 0; edge < 3; edge++) {
          Edge* tmpEdge = new Edge;
          v1 = octVertexOfEdgeOfParentFace[face][edge][0];
          v2 = octVertexOfEdgeOfParentFace[face][edge][1];
          tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
          newEdge[octEdgeFromParentFace[face][edge]] = tmpEdge;
          gridLevel.appendEdge(tmpEdge);
        }
        for (neighbor = 0; neighbor < 4; neighbor++) {
          neighIndex = octNeighborOfParentFace[face][neighbor];
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
        v1 = octVertexOfEdgeOfParentFace[face][edge][0];
        v2 = octVertexOfEdgeOfParentFace[face][edge][1];
        tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
        newEdge[octEdgeFromParentFace[face][edge]] = tmpEdge;
        gridLevel.appendEdge(tmpEdge);
      }
      for (neighbor = 0; neighbor < 4; neighbor++) {
        neighIndex = octNeighborOfParentFace[face][neighbor];
        newNeighbor[neighIndex] = 0;
      }
    }
  }


  // -* 4. create interior edges from parent edge to barycenter
  for (edge = 0; edge < 12; edge++) {
    Edge* tmpEdge = new Edge;
    int v1,v2;
    v1 = octVertexOfIntEdgeFromParentEdge[edge][0];
    v2 = octVertexOfIntEdgeFromParentEdge[edge][1];
    edgeIndex = octInteriorEdgeFromParentEdge[edge];
    tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
    newEdge[edgeIndex] = tmpEdge;
    gridLevel.appendEdge(tmpEdge);
  }


  // create children
  child = new Element*[14];

  child[0]  = new Octahedron;
  child[1]  = new Octahedron;
  child[2]  = new Octahedron;
  child[3]  = new Octahedron;
  child[4]  = new Octahedron;
  child[5]  = new Octahedron;
  child[6]  = new Tetrahedron;
  child[7]  = new Tetrahedron;
  child[8]  = new Tetrahedron;
  child[9]  = new Tetrahedron;
  child[10] = new Tetrahedron;
  child[11] = new Tetrahedron;
  child[12] = new Tetrahedron;
  child[13] = new Tetrahedron;

  gridLevel.appendOctahedron(child[0]);
  gridLevel.appendOctahedron(child[1]);
  gridLevel.appendOctahedron(child[2]);
  gridLevel.appendOctahedron(child[3]);
  gridLevel.appendOctahedron(child[4]);
  gridLevel.appendOctahedron(child[5]);
  gridLevel.appendTetrahedron(child[6]);
  gridLevel.appendTetrahedron(child[7]);
  gridLevel.appendTetrahedron(child[8]);
  gridLevel.appendTetrahedron(child[9]);
  gridLevel.appendTetrahedron(child[10]);
  gridLevel.appendTetrahedron(child[11]);
  gridLevel.appendTetrahedron(child[12]);
  gridLevel.appendTetrahedron(child[13]);

  // set newNeighbor for interior connections
  newNeighbor[32] = child[0];  // oct 0
  newNeighbor[33] = child[1];  // oct 1
  newNeighbor[34] = child[2];  // oct 2
  newNeighbor[35] = child[3];  // oct 3
  newNeighbor[36] = child[4];  // oct 4
  newNeighbor[37] = child[5];  // oct 5
  newNeighbor[38] = child[6];  // tet 0
  newNeighbor[39] = child[7];  // tet 1
  newNeighbor[40] = child[8];  // tet 2
  newNeighbor[41] = child[9];  // tet 3
  newNeighbor[42] = child[10]; // tet 4
  newNeighbor[43] = child[11]; // tet 5
  newNeighbor[44] = child[12]; // tet 6
  newNeighbor[45] = child[13]; // tet 7

  // set children parent
  child[0]->setParent(this);
  child[1]->setParent(this);
  child[2]->setParent(this);
  child[3]->setParent(this);
  child[4]->setParent(this);
  child[5]->setParent(this);
  child[6]->setParent(this);
  child[7]->setParent(this);
  child[8]->setParent(this);
  child[9]->setParent(this);
  child[10]->setParent(this);
  child[11]->setParent(this);
  child[12]->setParent(this);
  child[13]->setParent(this);

  // assemble octa children from new elements
  for (ch = 0; ch < 6; ch++) {
    // set vertices
    int vertexIndex;
    for (vertex = 0; vertex < 6; vertex++) {
      vertexIndex = octVertexOfOctaChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 12; edge++) {
      edgeIndex = octEdgeOfOctaChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
    //         connect child
    int neighIndex;
    for (face = 0; face < 8; face++) {
      neighIndex = octOctaChildNeighOfFace[ch][face];
      child[ch]->setNeighbor(face,newNeighbor[neighIndex]);

    }
  }

  // assemble tetra children
  for (ch = 0; ch < 8; ch++) {
    for (vertex = 0; vertex < 4; vertex++) {
      vertexIndex = octVertexOfTetraChild[ch][vertex];
      child[ch+6]->setVertex(vertex,newVertex[vertexIndex]);
    }
    for (edge = 0; edge < 6; edge++) {
      edgeIndex = octEdgeOfTetraChild[ch][edge];
      child[ch+6]->setEdge(edge,newEdge[edgeIndex]);
    }
    for (face = 0; face < 4; face++) {
      neighIndex = octTetraChildNeighborOfFace[ch][face];
      child[ch+6]->setNeighbor(face,newNeighbor[neighIndex]);
    }
  }


  // connect external neighbor with children
  for (neighbor = 0; neighbor < 32; neighbor++) {
    if (newNeighbor[neighbor])
    {
      childIndex = octChildNeighOfNeighbor[neighbor];
      newNeighbor[neighbor]->setNeighbor(neighFace[neighbor],child[childIndex]);
    }
  }

  // set the barycenter of the new octahedra
  Octahedron* childOcta;

  for (k = 0; k < 6; k++) {
    Vertex* tmpVertex = new Vertex;
    childOcta = (Octahedron*) child[k];
    childOcta->v[6] = tmpVertex;
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 6; j++) {
        (*(childOcta->v[6]))[i] += (*(childOcta->v[j]))[i];
      }
      (*(childOcta->v[6]))[i] *= 0.166666666666666;
    }
    gridLevel.appendVertex(tmpVertex);
  }

  // set tag
  setRefined();

}



// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Octahedron::coarsen(GridLevel& gridLevel)
{
  int ch;
  int face;
  int nFace;
  Element* theNeighbor;

  // *- octa children
  for (ch = 0; ch < 6; ch++) {
    for (face = 0; face < 8; face++) {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor) {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeOctahedron(child[ch]);
  }

  // *- tets children
  for (ch = 6; ch < 14; ch++) {
    // *- set neighbor's pointer
    //for (face = 0; face < 4; face++) {
    theNeighbor = child[ch]->getNeighbor(3);
    if (theNeighbor) {
      nFace = theNeighbor->getCommonFace(child[ch]);
      theNeighbor->setNeighbor(nFace,0);
    }
      //}
    gridLevel.removeTetrahedron(child[ch]);
  }

  for (ch = 0; ch < 14; ch++)
    delete child[ch];
  child = 0;


  // set tag
  setRegular();
}


// Description
//
void
Octahedron::close(int refPattern,ConformingClosure& closure)
{
  int noOfElem = 0;
  int tet;
  //int refPattern;

  Vertex* aa;
  Vertex* bb;
  Vertex* cc;

  // init element buffer
  closure.init();
  ConformingClosure::TET theElement = closure.beginTetrahedron();

  //refPattern = getRefinementPattern();

  for (tet = 0; tet < 8; tet++) {
    // close tetra
    switch(octSlicePattern[tet][refPattern]) {
    case 0:
      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);

      noOfElem++;
      break;
    case 1:
      aa = e[octEdgeSlicePattern[tet][0]]->getMidPoint();

      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,aa);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,aa);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);

      noOfElem++;
      break;
    case 2:
      aa = e[octEdgeSlicePattern[tet][1]]->getMidPoint();

      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,aa);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,aa);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);

      noOfElem++;
      break;
    case 3:
      aa = e[octEdgeSlicePattern[tet][0]]->getMidPoint();
      bb = e[octEdgeSlicePattern[tet][1]]->getMidPoint();

      if (v[octTetras[tet][0]]->getId() > v[octTetras[tet][1]]->getId()) {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,bb);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
      }
      else {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,bb);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
      }

      noOfElem++;
      break;
    case 4:
      aa = e[octEdgeSlicePattern[tet][2]]->getMidPoint();

      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,aa);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,aa);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);

      noOfElem++;
      break;
    case 5:
      aa = e[octEdgeSlicePattern[tet][0]]->getMidPoint();
      bb = e[octEdgeSlicePattern[tet][2]]->getMidPoint();

      if (v[octTetras[tet][0]]->getId() > v[octTetras[tet][2]]->getId()) {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,bb);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
      }
      else {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,bb);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
      }
      noOfElem++;
      break;
    case 6:
      aa = e[octEdgeSlicePattern[tet][1]]->getMidPoint();
      bb = e[octEdgeSlicePattern[tet][2]]->getMidPoint();

      if (v[octTetras[tet][1]]->getId() > v[octTetras[tet][2]]->getId()) {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,bb);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,aa);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
      }
      else {
        theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
        theElement[noOfElem]->setVertex(1,bb);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,aa);
        theElement[noOfElem]->setVertex(1,bb);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
        noOfElem++;
        theElement[noOfElem]->setVertex(0,bb);
        theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
        theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
        theElement[noOfElem]->setVertex(3,v[6]);
      }

      noOfElem++;
      break;
    case 7:
      aa = e[octEdgeSlicePattern[tet][0]]->getMidPoint();
      bb = e[octEdgeSlicePattern[tet][1]]->getMidPoint();
      cc = e[octEdgeSlicePattern[tet][2]]->getMidPoint();

      theElement[noOfElem]->setVertex(0,v[octTetras[tet][0]]);
      theElement[noOfElem]->setVertex(1,cc);
      theElement[noOfElem]->setVertex(2,bb);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,cc);
      theElement[noOfElem]->setVertex(1,v[octTetras[tet][1]]);
      theElement[noOfElem]->setVertex(2,aa);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,bb);
      theElement[noOfElem]->setVertex(1,aa);
      theElement[noOfElem]->setVertex(2,v[octTetras[tet][2]]);
      theElement[noOfElem]->setVertex(3,v[6]);
      noOfElem++;
      theElement[noOfElem]->setVertex(0,cc);
      theElement[noOfElem]->setVertex(1,aa);
      theElement[noOfElem]->setVertex(2,bb);
      theElement[noOfElem]->setVertex(3,v[6]);

      noOfElem++;
      break;
    default :
      break;
    }
  }

  closure.setNoOfTetrahedra(noOfElem);
}



// ____________________________________________________________
//
//	
//	         functions
//
// ____________________________________________________________
//


//
//
void
Octahedron::close(ConformingClosure& closure)
{
  int refPattern = getRefinementPattern();
  close(refPattern,closure);
}


// Description
// check onle tetrahedra, child at middle of face
bool
Octahedron::queryChildrenRefinement()
{
  if (isRefined())
  {
    for (int i = 6; i < 14; i++)
    {
      if (child[i]->isRefined())
        return true;
    }
    return false;
  }
  else
    return false;
}


// Description
//
bool
Octahedron::queryEdgeRefinement()
{
  for (int i = 0; i < 12; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Octahedron::queryChildrenEdgeRefinement()
{
  int i;
  for (i = 0; i < 12; i++) {
    if (e[i]->queryChildrenRefinement())
      return true;
  }

  // check if neighbor's child is refined
  for (i = 0; i < 8; i++) {
    if (n[i]) {
      if (n[i]->queryChildrenRefinement())
        return true;
    }
  }

  return false;
}



// include topology information
#include "OctahedronTopology.h"


int
Octahedron::getCommonFace(Element* theElem)
{
  for (int face = 0; face < 8; face++)
    if (n[face] == theElem)
      return face;

  return -1;
}



// Description
//
void
Octahedron::getRefElements(int face,Vertex** theVertex,Edge** theEdge,
                           Element** theChild,int* theChildFace)
{
  int vertex,index, childIndex;

  for (vertex = 0; vertex < 3; vertex++) {
    index = octVertexOfFace[face][vertex];
    theVertex[vertex] = v[index];

    childIndex = octChildOfParentFace[face][vertex];
    theChild[vertex] = child[childIndex];
    theChildFace[vertex] = octChildFaceOfParentFace[face][vertex];

    index = octChildEdgeOfParentFace[face][vertex];
    theEdge[vertex] = child[childIndex]->getEdge(index);
  }

  childIndex = octChildOfParentFace[face][3];
  theChild[3] = child[childIndex];
  theChildFace[3] = octChildFaceOfParentFace[face][3];

}



// Description
//
void
Octahedron::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
void
Octahedron::setVertices(Vertex** vt)
{
  for (int i = 0; i < 6; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Octahedron::getVertex(int i)
{
  return v[i];
}


// Description
//
void
Octahedron::getVertices(Vertex** vt)
{
  for (int i = 0; i < 6; i++)
    vt[i] = v[i];
}


// Description
//
void
Octahedron::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
void
Octahedron::setEdges(Edge** edg)
{
  for (int i = 0; i < 12; i++)
    e[i] = edg[i];
}


// Description
//
Edge*
Octahedron::getEdge(int i)
{
  return e[i];
}


// Description
//
void
Octahedron::getEdges(Edge** edg)
{
  for (int i = 0; i < 6; i++)
    edg[i] = e[i];
}


// Description
//
void
Octahedron::setNeighbor(int face,Element* theNeighbor)
{
  n[face] = theNeighbor;
}


// Description
//
Element*
Octahedron::getNeighbor(int face)
{
  return n[face];
}


//
// Topology
//

//
//
int
Octahedron::type()
{
  return GRD_OCTAHEDRON;
}


// Description
//
int
Octahedron::noOfVertexAtFace(int fc)
{
  return 3;
}


// Description
//
Vertex*
Octahedron::getVertexAtFace(int fc,int vt)
{
  return (v[octVertexOfFace[fc][vt]]);
}

//
//
Vertex*
Octahedron::getVertexAtEdge(int ed,int vt)
{
  return (v[octVertexOfEdge[ed][vt]]);
}

//
//
Edge*
Octahedron::getEdgeAtFace(int fc,int ed)
{
  return e[octEdgeOfFace[fc][ed]];
}

} // namespace grd
