// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/***************************************************************************
    File        : Octahedron.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Octahedron.hh"
#include "Tetrahedron.hh"
#include "GridLevel.hh"

namespace grd {

// initialize memory management here
//Pool<Octahedron> Octahedron::pool;


// Description:
//
Octahedron::Octahedron( )
{
  int i;
  for (i = 0; i < 6; i++)
    v[i] = 0;

  for (i = 0; i < 12; i++)
    e[i] = 0;

  for (i = 0; i < 8; i++)
    n[i] = 0;

  // set no. of children
  // chcount = 0; <- set in abstract class Element
}


// Description
//
void
Octahedron::getTetras(std::vector<Element*>& tetra) const
{
  // set tetra 0
  tetra[0]->setVertex(0,v[0]);
  tetra[0]->setVertex(1,v[1]);
  tetra[0]->setVertex(2,v[5]);
  tetra[0]->setVertex(3,v[2]);
  // set neighbors
  tetra[0]->setNeighbor(0,n[4]);
  tetra[0]->setNeighbor(1,tetra[1]);
  tetra[0]->setNeighbor(2,n[0]);
  tetra[0]->setNeighbor(3,tetra[3]);


  // set tetra 1
  tetra[1]->setVertex(0,v[0]);
  tetra[1]->setVertex(1,v[2]);
  tetra[1]->setVertex(2,v[5]);
  tetra[1]->setVertex(3,v[3]);
  // set neighbors
  tetra[1]->setNeighbor(0,n[5]);
  tetra[1]->setNeighbor(1,tetra[2]);
  tetra[1]->setNeighbor(2,n[1]);
  tetra[1]->setNeighbor(3,tetra[0]);

  // set tetra 2
  tetra[2]->setVertex(0,v[0]);
  tetra[2]->setVertex(1,v[3]);
  tetra[2]->setVertex(2,v[5]);
  tetra[2]->setVertex(3,v[4]);
  // set neighbors
  tetra[2]->setNeighbor(0,n[6]);
  tetra[2]->setNeighbor(1,tetra[3]);
  tetra[2]->setNeighbor(2,n[2]);
  tetra[2]->setNeighbor(3,tetra[1]);

  // set tetra 3
  tetra[3]->setVertex(0,v[0]);
  tetra[3]->setVertex(1,v[4]);
  tetra[3]->setVertex(2,v[5]);
  tetra[3]->setVertex(3,v[1]);
  // set neighbors
  tetra[3]->setNeighbor(0,n[7]);
  tetra[3]->setNeighbor(1,tetra[0]);
  tetra[3]->setNeighbor(2,n[3]);
  tetra[3]->setNeighbor(3,tetra[2]);

}


// Description
//  Spit an Octahedron into four
//  tetrahedra and store the new elements
//  in the same grid level
void
Octahedron::split(GridLevel& grid)
{
  int cface;
  // Create new tetrahedra
  Element* tetra[4];
  tetra[0] = new Tetrahedron;
  tetra[1] = new Tetrahedron;
  tetra[2] = new Tetrahedron;
  tetra[3] = new Tetrahedron;
  Edge* medge = new Edge(v[0],v[5]);

  // set tetra 0
  tetra[0]->setVertex(0,v[0]);
  tetra[0]->setVertex(1,v[1]);
  tetra[0]->setVertex(2,v[5]);
  tetra[0]->setVertex(3,v[2]);
  // set neighbors
  tetra[0]->setNeighbor(0,n[4]);
  tetra[0]->setNeighbor(1,tetra[1]);
  tetra[0]->setNeighbor(2,n[0]);
  tetra[0]->setNeighbor(3,tetra[3]);
  // connect neighbors with tets
  cface = n[4]->getCommonFace(this);
  n[4]->setNeighbor(cface,tetra[0]);
  cface = n[0]->getCommonFace(this);
  n[0]->setNeighbor(cface,tetra[0]);
  // set edges
  tetra[0]->setEdge(0,medge);
  tetra[0]->setEdge(1,e[8]);
  tetra[0]->setEdge(2,e[0]);
  tetra[0]->setEdge(3,e[9]);
  tetra[0]->setEdge(4,e[4]);
  tetra[0]->setEdge(5,e[1]);

  // set tetra 1
  tetra[1]->setVertex(0,v[0]);
  tetra[1]->setVertex(1,v[2]);
  tetra[1]->setVertex(2,v[5]);
  tetra[1]->setVertex(3,v[3]);
  // set neighbors
  tetra[1]->setNeighbor(0,n[5]);
  tetra[1]->setNeighbor(1,tetra[2]);
  tetra[1]->setNeighbor(2,n[1]);
  tetra[1]->setNeighbor(3,tetra[0]);
  // connect neighbors with tets
  cface = n[5]->getCommonFace(this);
  n[5]->setNeighbor(cface,tetra[1]);
  cface = n[1]->getCommonFace(this);
  n[1]->setNeighbor(cface,tetra[1]);
  // set edges
  tetra[1]->setEdge(0,medge);
  tetra[1]->setEdge(1,e[9]);
  tetra[1]->setEdge(2,e[1]);
  tetra[1]->setEdge(3,e[10]);
  tetra[1]->setEdge(4,e[5]);
  tetra[1]->setEdge(5,e[2]);
  
  // set tetra 2
  tetra[2]->setVertex(0,v[0]);
  tetra[2]->setVertex(1,v[3]);
  tetra[2]->setVertex(2,v[5]);
  tetra[2]->setVertex(3,v[4]);
  // set neighbors
  tetra[2]->setNeighbor(0,n[6]);
  tetra[2]->setNeighbor(1,tetra[3]);
  tetra[2]->setNeighbor(2,n[2]);
  tetra[2]->setNeighbor(3,tetra[1]);
  // connect neighbors with tets
  cface = n[6]->getCommonFace(this);
  n[6]->setNeighbor(cface,tetra[2]);
  cface = n[2]->getCommonFace(this);
  n[2]->setNeighbor(cface,tetra[2]);
  // set edges
  tetra[2]->setEdge(0,medge);
  tetra[2]->setEdge(1,e[10]);
  tetra[2]->setEdge(2,e[2]);
  tetra[2]->setEdge(3,e[11]);
  tetra[2]->setEdge(4,e[6]);
  tetra[2]->setEdge(5,e[3]);

  // set tetra 3
  tetra[3]->setVertex(0,v[0]);
  tetra[3]->setVertex(1,v[4]);
  tetra[3]->setVertex(2,v[5]);
  tetra[3]->setVertex(3,v[1]);
  // set neighbors
  tetra[3]->setNeighbor(0,n[7]);
  tetra[3]->setNeighbor(1,tetra[0]);
  tetra[3]->setNeighbor(2,n[3]);
  tetra[3]->setNeighbor(3,tetra[2]);
  // connect neighbors with tets
  cface = n[7]->getCommonFace(this);
  n[7]->setNeighbor(cface,tetra[3]);
  cface = n[3]->getCommonFace(this);
  n[3]->setNeighbor(cface,tetra[3]);
  // set edges
  tetra[3]->setEdge(0,medge);
  tetra[3]->setEdge(1,e[11]);
  tetra[3]->setEdge(2,e[3]);
  tetra[3]->setEdge(3,e[8]);
  tetra[3]->setEdge(4,e[7]);
  tetra[3]->setEdge(5,e[0]);

  // Append elements to octahedron
  child = new Element*[4];
  child[0] = tetra[0]; child[0]->setParent(this);
  child[1] = tetra[1]; child[1]->setParent(this);
  child[2] = tetra[2]; child[2]->setParent(this);
  child[3] = tetra[3]; child[3]->setParent(this);
  // Inherit value from parent
  child[0]->setValue(getValue());
  child[1]->setValue(getValue());
  child[2]->setValue(getValue());
  child[3]->setValue(getValue());

  // Insert into grild level
  grid.appendTetrahedron(tetra[0]);
  grid.appendTetrahedron(tetra[1]);
  grid.appendTetrahedron(tetra[2]);
  grid.appendTetrahedron(tetra[3]);
  grid.appendEdge(medge);
}


// Description
// Marks split tetra
// they will be removed later
void
Octahedron::unsplit()
{
  int cface; 
 
  // mark edge for unrefinement
  Edge* medge = child[0]->getEdge(0);
  medge->ref();

  // Set enighborhood
  // tet 0
  n[4] = child[0]->getNeighbor(0);
  cface = n[4]->getCommonFace(child[0]);
  n[4]->setNeighbor(cface,this);
  n[0] = child[0]->getNeighbor(2);
  cface = n[0]->getCommonFace(child[0]);
  n[0]->setNeighbor(cface,this);
  
  // set tetra 1
  n[5] = child[1]->getNeighbor(0);
  cface = n[5]->getCommonFace(child[1]);
  n[5]->setNeighbor(cface,this);
  n[1] = child[1]->getNeighbor(2);
  cface = n[1]->getCommonFace(child[1]);
  n[1]->setNeighbor(cface,this);

  // set tetra 2
  n[6] = child[2]->getNeighbor(0);
  cface = n[6]->getCommonFace(child[2]);
  n[6]->setNeighbor(cface,this);
  n[2] = child[2]->getNeighbor(2);
  cface = n[2]->getCommonFace(child[2]);
  n[2]->setNeighbor(cface,this);

  // set tetra 3
  n[7] = child[3]->getNeighbor(0);
  cface = n[7]->getCommonFace(child[3]);
  n[7]->setNeighbor(cface,this);
  n[3] = child[3]->getNeighbor(2);
  cface = n[3]->getCommonFace(child[3]);
  n[3]->setNeighbor(cface,this);

  // evaluate children marks
  if (child[0]->isMarkedForCoarsening()) markForCoarsening();
  else if (child[0]->isMarkedForRefinement()) markForRefinement();
  if (child[1]->isMarkedForCoarsening()) markForCoarsening();
  else if (child[1]->isMarkedForRefinement()) markForRefinement();
  if (child[2]->isMarkedForCoarsening()) markForCoarsening();
  else if (child[2]->isMarkedForRefinement()) markForRefinement();
  if (child[3]->isMarkedForCoarsening()) markForCoarsening();
  else if (child[3]->isMarkedForRefinement()) markForRefinement();
  
  // Reference elements to be romeved
  child[0]->ref();
  child[1]->ref();
  child[2]->ref();
  child[3]->ref();
  
  // Append elements to octahedron
  delete [] child;
  child = 0;
}


//=============================================================================
//  REFINEMENT
//=============================================================================

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
  int i,j;
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
  Vertex* bc = new Vertex;
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 6; j++)
    {
      (*bc)[i] += (*(v[j]))[i];
    }
    (*bc)[i] *= ONE_SIXTH;
  }
  newVertex[vertCount] = bc;
  gridLevel.appendVertex(bc);

  // -* 3. process parent faces

  for (face = 0; face < 8; face++) {
    // if neighbor regular refined
    if (n[face]) {
      if (n[face]->isRefined()) {
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

  // set tag
  setRefined();
  // Set no. of children
  chcount = 14;

  // set value
  for (ch = 0; ch < 14; ch++) child[ch]->setValue(getValue());
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
  for (ch = 0; ch < 6; ch++)
  {
    for (face = 0; face < 8; face++)
    {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor) {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeOctahedron(child[ch]);
  }

  // *- tets children
  for (ch = 6; ch < 14; ch++)
  {
    // *- set neighbor's pointer
    //for (face = 0; face < 4; face++) {
    theNeighbor = child[ch]->getNeighbor(3);
    if (theNeighbor)
    {
      nFace = theNeighbor->getCommonFace(child[ch]);
      theNeighbor->setNeighbor(nFace,0);
    }
      //}
    gridLevel.removeTetrahedron(child[ch]);
  }

  // Delete the children from memory
  for (ch = 0; ch < 14; ch++)
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
Octahedron::close(GridLevel& gridLevel)
{
  size_t refPattern = 0;
  size_t refMask    = 1;
  size_t noOfElem;
  int i0,i1,i2,i3;

  Vertex* aa;
  Vertex* bb;
  Vertex* cc;
  Vertex* dd;
  Vertex* ee;

  // Element buffer
  Element* theElement[28];

  // init count element buffer
  noOfElem = 0;

  // Edge Refinement Pattern
  for (size_t i = 0; i < 12; i++)
  {
    if (e[i]->isRefined())
    { refPattern = (refPattern | refMask); }

    refMask = (refMask<<1);
  }

  for (int tet = 0; tet < 4; tet++)
  {
    switch(octTetRefPattern[tet][refPattern]) {
// CASE 0
      case 0:
        theElement[noOfElem] = new Tetrahedron;

        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(2,v[5]);
        theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem++;
        break;
// CASE 2
      case 2:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+1]->setVertex(0,v[0]);
        theElement[noOfElem+1]->setVertex(1,aa);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem += 2;
        break;
// CASE 4
      case 4:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[5]);
        theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+1]->setVertex(0,aa);
        theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem += 2;
        break;
// CASE 6
      case 6:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        if ((v[0]->getId()) > v[5]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,bb);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,aa);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,bb);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,aa);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 3;
        break;
// CASE 8
      case 8:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;

        aa = e[octTetEdge[tet][3]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+1]->setVertex(0,v[0]);
        theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,aa);
        noOfElem += 2;
        break;
// CASE 10
      case 10:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        if ((v[octTetVertex[tet][1]]->getId()) > (v[octTetVertex[tet][3]]->getId()))
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,bb);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,aa);
          theElement[noOfElem+1]->setVertex(3,bb);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,bb);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,bb);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,bb);
        }
        noOfElem += 3;
        break;
// CASE 12
      case 12:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[5]);
        theElement[noOfElem]->setVertex(3,bb);

        theElement[noOfElem+1]->setVertex(0,v[0]);
        theElement[noOfElem+1]->setVertex(1,aa);
        theElement[noOfElem+1]->setVertex(2,bb);
        theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+2]->setVertex(0,aa);
        theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+2]->setVertex(2,v[5]);
        theElement[noOfElem+2]->setVertex(3,bb);

        theElement[noOfElem+3]->setVertex(0,aa);
        theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+3]->setVertex(2,bb);
        theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem += 4;
        break;
// CASE 14
      case 14:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;
        theElement[noOfElem+4] = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][3]]->getMidPoint();

        if (i0 > i2)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
        }
        else //if (i2 > i0)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        noOfElem += 5;
        break;
// CASE 16
      case 16:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;

        aa = e[octTetEdge[tet][4]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,aa);
        theElement[noOfElem]->setVertex(2,v[5]);
        theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+1]->setVertex(0,v[0]);
        theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,aa);
        noOfElem += 2;
        break;
// CASE 18
      case 18:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        if ((v[5]->getId()) > v[octTetVertex[tet][3]]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,bb);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,bb);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,bb);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 3;
        break;
// CASE 20
      case 20:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        if ((v[0]->getId()) > v[octTetVertex[tet][3]]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,bb);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,bb);

          theElement[noOfElem+2]->setVertex(0,aa);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,bb);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,aa);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+1]->setVertex(0,aa);
          theElement[noOfElem+1]->setVertex(1,bb);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,aa);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,bb);
        }
        noOfElem += 3;
        break;
// CASE 22
      case 22:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;

        i0 = v[0]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();

        if (i0 > i2 && i0 > i3)
        {
          if (i2 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,cc);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,cc);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);
          }
        }
        else if (i2 > i0 && i2 > i3)
        {
          if (i0 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,cc);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);
          }
        }
        else if (i3 > i0 && i3 > i2)
        {
          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        noOfElem += 4;
        break;
// CASE 24
      case 24:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;

        aa = e[octTetEdge[tet][3]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        if ((v[octTetVertex[tet][1]]->getId()) > v[5]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,aa);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,bb);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,bb);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,aa);

          theElement[noOfElem+2]->setVertex(0,v[0]);
          theElement[noOfElem+2]->setVertex(1,bb);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 3;
        break;
// CASE 26
      case 26:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();

        theElement[noOfElem]->setVertex(0,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(1,cc);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[0]);

        theElement[noOfElem+1]->setVertex(0,aa);
        theElement[noOfElem+1]->setVertex(1,bb);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,v[0]);

        theElement[noOfElem+2]->setVertex(0,cc);
        theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][3]]);
        theElement[noOfElem+2]->setVertex(2,bb);
        theElement[noOfElem+2]->setVertex(3,v[0]);

        theElement[noOfElem+3]->setVertex(0,aa);
        theElement[noOfElem+3]->setVertex(1,cc);
        theElement[noOfElem+3]->setVertex(2,bb);
        theElement[noOfElem+3]->setVertex(3,v[0]);
        noOfElem += 4;
        break;
// CASE 28
      case 28:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;
        theElement[noOfElem+4] = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();

        if (i0 > i3)
        {
          if (i1 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,cc);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,cc);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[5]);
            theElement[noOfElem+4]->setVertex(2,bb);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
        }
        else //if (i3 > i0)
        {
          if (i1 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,aa);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
        }
        noOfElem += 5;
        break;
// CASE 30
      case 30:
        theElement[noOfElem]   = new Tetrahedron;
        theElement[noOfElem+1] = new Tetrahedron;
        theElement[noOfElem+2] = new Tetrahedron;
        theElement[noOfElem+3] = new Tetrahedron;
        theElement[noOfElem+4] = new Tetrahedron;

        i0 = v[0]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][3]]->getMidPoint();
        dd = e[octTetEdge[tet][4]]->getMidPoint();

        if (i0 > i3)
        {
          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,dd);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,dd);
            theElement[noOfElem+4]->setVertex(2,cc);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 5;
          }
          else
          {
            theElement[noOfElem+5]   = new Tetrahedron;

            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,dd);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,dd);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,bb);
            theElement[noOfElem+5]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+5]->setVertex(2,aa);
            theElement[noOfElem+5]->setVertex(3,dd);

            noOfElem += 6;
          }
        }
        else
        {
          theElement[noOfElem+5]   = new Tetrahedron;
          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,dd);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,dd);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,bb);
            theElement[noOfElem+5]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+5]->setVertex(2,aa);
            theElement[noOfElem+5]->setVertex(3,dd);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,dd);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,dd);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,bb);
            theElement[noOfElem+5]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+5]->setVertex(2,aa);
            theElement[noOfElem+5]->setVertex(3,dd);
          }
          noOfElem += 6;
        }
        break;
// CASE 32
      case 32:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;

        aa = e[octTetEdge[tet][5]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(2,v[5]);
        theElement[noOfElem]->setVertex(3,aa);

        theElement[noOfElem+1]->setVertex(0,aa);
        theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(2,v[5]);
        theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem += 2;
        break;
// CASE 34
      case 34:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][5]]->getMidPoint();
        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,bb);

        theElement[noOfElem+1]->setVertex(0,bb);
        theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(2,aa);
        theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

        theElement[noOfElem+2]->setVertex(0,v[0]);
        theElement[noOfElem+2]->setVertex(1,aa);
        theElement[noOfElem+2]->setVertex(2,v[5]);
        theElement[noOfElem+2]->setVertex(3,bb);

        theElement[noOfElem+3]->setVertex(0,bb);
        theElement[noOfElem+3]->setVertex(1,aa);
        theElement[noOfElem+3]->setVertex(2,v[5]);
        theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
        noOfElem += 4;
        break;
// CASE 36
      case 36:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][5]]->getMidPoint();
        if ((v[octTetVertex[tet][1]]->getId()) > v[octTetVertex[tet][3]]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,aa);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,aa);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,bb);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,aa);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,aa);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 3;
        break;
// CASE 38
      case 38:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        else //if (i2 > i0)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        noOfElem += 5;
        break;
// CASE 40
      case 40:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;

        aa = e[octTetEdge[tet][3]]->getMidPoint();
        bb = e[octTetEdge[tet][5]]->getMidPoint();
        if ((v[0]->getId()) > v[5]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,aa);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,aa);
          theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,aa);
        }
        noOfElem += 3;
        break;
// CASE 42
      case 42:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
        }
        else //if (i2 > i0)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
        }
        noOfElem += 5;
        break;
// CASE 44
      case 44:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,aa);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
        }
        else //if (i2 > i0)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,aa);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);
          }
        }
        noOfElem += 5;
        break;
// CASE 46
      case 46:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;
        theElement[noOfElem+5]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][3]]->getMidPoint();
        dd = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,dd);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,cc);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,v[0]);
            theElement[noOfElem+5]->setVertex(1,aa);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,cc);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,dd);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,cc);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,bb);
            theElement[noOfElem+5]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+5]->setVertex(2,aa);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        else {
          if (i1 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,cc);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,cc);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,bb);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,cc);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,dd);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,cc);
            theElement[noOfElem+1]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+4]->setVertex(0,bb);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,bb);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,cc);
          }
        }
        noOfElem += 6;
        break;
// CASE 48
      case 48:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;

        aa = e[octTetEdge[tet][4]]->getMidPoint();
        bb = e[octTetEdge[tet][5]]->getMidPoint();
        if ((v[0]->getId()) > v[octTetVertex[tet][1]]->getId())
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,aa);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,v[0]);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,aa);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,bb);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,aa);

          theElement[noOfElem+2]->setVertex(0,bb);
          theElement[noOfElem+2]->setVertex(1,aa);
          theElement[noOfElem+2]->setVertex(2,v[5]);
          theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 3;
        break;
// CASE 50
      case 50:
      theElement[noOfElem]     = new Tetrahedron;
      theElement[noOfElem+1]   = new Tetrahedron;
      theElement[noOfElem+2]   = new Tetrahedron;
      theElement[noOfElem+3]   = new Tetrahedron;
      theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i1)
        {
          if (i2 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,bb);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[0]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,bb);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,bb);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[0]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
        }
        else //if (i1 > i0)
        {
          if (i2 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,cc);
            theElement[noOfElem+4]->setVertex(1,bb);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);
          }
        }
        noOfElem += 5;
        break;
// CASE 52
      case 52:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        theElement[noOfElem]->setVertex(0,v[0]);
        theElement[noOfElem]->setVertex(1,cc);
        theElement[noOfElem]->setVertex(2,aa);
        theElement[noOfElem]->setVertex(3,v[5]);

        theElement[noOfElem+1]->setVertex(0,aa);
        theElement[noOfElem+1]->setVertex(1,bb);
        theElement[noOfElem+1]->setVertex(2,v[octTetVertex[tet][1]]);
        theElement[noOfElem+1]->setVertex(3,v[5]);

        theElement[noOfElem+2]->setVertex(0,cc);
        theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][3]]);
        theElement[noOfElem+2]->setVertex(2,bb);
        theElement[noOfElem+2]->setVertex(3,v[5]);

        theElement[noOfElem+3]->setVertex(0,aa);
        theElement[noOfElem+3]->setVertex(1,cc);
        theElement[noOfElem+3]->setVertex(2,bb);
        theElement[noOfElem+3]->setVertex(3,v[5]);
        noOfElem += 4;
        break;
// CASE 54
      case 54:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i2 = v[5]->getId();
        i3 = v[octTetVertex[tet][3]]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();
        dd = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          theElement[noOfElem+5]   = new Tetrahedron;

          if (i2 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,dd);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,dd);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,cc);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,dd);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+4]->setVertex(0,v[0]);
            theElement[noOfElem+4]->setVertex(1,aa);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,dd);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,aa);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          noOfElem += 6;
        }
        else
        {
          if (i2 > i3)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,cc);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 5;
          }
          else
          {
            theElement[noOfElem+5]   = new Tetrahedron;

            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,bb);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,bb);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,bb);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,bb);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,cc);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,aa);
            theElement[noOfElem+5]->setVertex(2,v[5]);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 6;
          }
        }
        break;
// CASE 56
      case 56:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();

        aa = e[octTetEdge[tet][3]]->getMidPoint();
        bb = e[octTetEdge[tet][4]]->getMidPoint();
        cc = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i1 && i0 > i2)
        {
          if (i1 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,bb);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[0]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,aa);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,bb);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,bb);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,aa);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[0]);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        else if (i1 > i0 && i1 > i2)
        {
          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,aa);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,aa);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        else if (i2 > i0 && i2 > i1)
        {
          if (i0 > i1)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,bb);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[0]);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,aa);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,cc);
            theElement[noOfElem+1]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,cc);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,aa);

            theElement[noOfElem+3]->setVertex(0,cc);
            theElement[noOfElem+3]->setVertex(1,bb);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,v[octTetVertex[tet][3]]);
          }
        }
        noOfElem += 4;
        break;
// CASE 58
      case 58:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();
        dd = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i1)
        {
          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,cc);
            theElement[noOfElem+4]->setVertex(2,bb);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 5;
          }
          else
          {
            theElement[noOfElem+5]   = new Tetrahedron;

            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,cc);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,aa);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,v[0]);
            theElement[noOfElem+3]->setVertex(1,aa);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,dd);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,bb);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,dd);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 6;
          }
        }
        else
        {
          theElement[noOfElem+5]   = new Tetrahedron;

          if (i0 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,v[0]);
            theElement[noOfElem+2]->setVertex(1,aa);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,bb);

            theElement[noOfElem+3]->setVertex(0,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,aa);
            theElement[noOfElem+3]->setVertex(3,dd);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,cc);
            theElement[noOfElem+4]->setVertex(2,aa);
            theElement[noOfElem+4]->setVertex(3,bb);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem]->setVertex(2,aa);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,aa);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,dd);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,bb);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,dd);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          noOfElem += 6;
        }
        break;
// CASE 60
      case 60:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;

        i0 = v[0]->getId();
        i1 = v[octTetVertex[tet][1]]->getId();
        i2 = v[5]->getId();

        aa = e[octTetEdge[tet][2]]->getMidPoint();
        bb = e[octTetEdge[tet][3]]->getMidPoint();
        cc = e[octTetEdge[tet][4]]->getMidPoint();
        dd = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          theElement[noOfElem+5]   = new Tetrahedron;

          if (i1 > i2)
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,bb);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,v[0]);
            theElement[noOfElem+1]->setVertex(1,aa);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,bb);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,cc);
            theElement[noOfElem+2]->setVertex(2,bb);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,v[5]);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,cc);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);
          }
          noOfElem += 6;
        }
        else
        {
          if (i1 > i2)
          {
            theElement[noOfElem+5]   = new Tetrahedron;

            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,aa);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,bb);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,bb);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,dd);

            theElement[noOfElem+3]->setVertex(0,aa);
            theElement[noOfElem+3]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+3]->setVertex(2,bb);
            theElement[noOfElem+3]->setVertex(3,cc);

            theElement[noOfElem+4]->setVertex(0,aa);
            theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+4]->setVertex(2,v[5]);
            theElement[noOfElem+4]->setVertex(3,bb);

            theElement[noOfElem+5]->setVertex(0,dd);
            theElement[noOfElem+5]->setVertex(1,cc);
            theElement[noOfElem+5]->setVertex(2,bb);
            theElement[noOfElem+5]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 6;
          }
          else
          {
            theElement[noOfElem]->setVertex(0,v[0]);
            theElement[noOfElem]->setVertex(1,aa);
            theElement[noOfElem]->setVertex(2,v[5]);
            theElement[noOfElem]->setVertex(3,dd);

            theElement[noOfElem+1]->setVertex(0,aa);
            theElement[noOfElem+1]->setVertex(1,cc);
            theElement[noOfElem+1]->setVertex(2,v[5]);
            theElement[noOfElem+1]->setVertex(3,dd);

            theElement[noOfElem+2]->setVertex(0,aa);
            theElement[noOfElem+2]->setVertex(1,v[octTetVertex[tet][1]]);
            theElement[noOfElem+2]->setVertex(2,v[5]);
            theElement[noOfElem+2]->setVertex(3,cc);

            theElement[noOfElem+3]->setVertex(0,dd);
            theElement[noOfElem+3]->setVertex(1,cc);
            theElement[noOfElem+3]->setVertex(2,v[5]);
            theElement[noOfElem+3]->setVertex(3,bb);

            theElement[noOfElem+4]->setVertex(0,dd);
            theElement[noOfElem+4]->setVertex(1,cc);
            theElement[noOfElem+4]->setVertex(2,bb);
            theElement[noOfElem+4]->setVertex(3,v[octTetVertex[tet][3]]);

            noOfElem += 5;
          }
        }
        break;
// CASE 62
      case 62:
        theElement[noOfElem]     = new Tetrahedron;
        theElement[noOfElem+1]   = new Tetrahedron;
        theElement[noOfElem+2]   = new Tetrahedron;
        theElement[noOfElem+3]   = new Tetrahedron;
        theElement[noOfElem+4]   = new Tetrahedron;
        theElement[noOfElem+5]   = new Tetrahedron;
        theElement[noOfElem+6]   = new Tetrahedron;

        i0 = v[0]->getId();
        i2 = v[5]->getId();

        aa = e[octTetEdge[tet][1]]->getMidPoint();
        bb = e[octTetEdge[tet][2]]->getMidPoint();
        cc = e[octTetEdge[tet][3]]->getMidPoint();
        dd = e[octTetEdge[tet][4]]->getMidPoint();
        ee = e[octTetEdge[tet][5]]->getMidPoint();

        if (i0 > i2)
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,bb);
          theElement[noOfElem]->setVertex(2,aa);
          theElement[noOfElem]->setVertex(3,ee);

          theElement[noOfElem+1]->setVertex(0,ee);
          theElement[noOfElem+1]->setVertex(1,bb);
          theElement[noOfElem+1]->setVertex(2,aa);
          theElement[noOfElem+1]->setVertex(3,dd);

          theElement[noOfElem+2]->setVertex(0,ee);
          theElement[noOfElem+2]->setVertex(1,dd);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,cc);

          theElement[noOfElem+3]->setVertex(0,v[0]);
          theElement[noOfElem+3]->setVertex(1,aa);
          theElement[noOfElem+3]->setVertex(2,cc);
          theElement[noOfElem+3]->setVertex(3,ee);

          theElement[noOfElem+4]->setVertex(0,bb);
          theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+4]->setVertex(2,aa);
          theElement[noOfElem+4]->setVertex(3,dd);

          theElement[noOfElem+5]->setVertex(0,v[0]);
          theElement[noOfElem+5]->setVertex(1,aa);
          theElement[noOfElem+5]->setVertex(2,v[5]);
          theElement[noOfElem+5]->setVertex(3,cc);

          theElement[noOfElem+6]->setVertex(0,ee);
          theElement[noOfElem+6]->setVertex(1,dd);
          theElement[noOfElem+6]->setVertex(2,cc);
          theElement[noOfElem+6]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        else
        {
          theElement[noOfElem]->setVertex(0,v[0]);
          theElement[noOfElem]->setVertex(1,bb);
          theElement[noOfElem]->setVertex(2,v[5]);
          theElement[noOfElem]->setVertex(3,ee);

          theElement[noOfElem+1]->setVertex(0,bb);
          theElement[noOfElem+1]->setVertex(1,aa);
          theElement[noOfElem+1]->setVertex(2,v[5]);
          theElement[noOfElem+1]->setVertex(3,ee);

          theElement[noOfElem+2]->setVertex(0,ee);
          theElement[noOfElem+2]->setVertex(1,bb);
          theElement[noOfElem+2]->setVertex(2,aa);
          theElement[noOfElem+2]->setVertex(3,dd);

          theElement[noOfElem+3]->setVertex(0,ee);
          theElement[noOfElem+3]->setVertex(1,dd);
          theElement[noOfElem+3]->setVertex(2,aa);
          theElement[noOfElem+3]->setVertex(3,cc);

          theElement[noOfElem+4]->setVertex(0,bb);
          theElement[noOfElem+4]->setVertex(1,v[octTetVertex[tet][1]]);
          theElement[noOfElem+4]->setVertex(2,aa);
          theElement[noOfElem+4]->setVertex(3,dd);

          theElement[noOfElem+5]->setVertex(0,ee);
          theElement[noOfElem+5]->setVertex(1,aa);
          theElement[noOfElem+5]->setVertex(2,v[5]);
          theElement[noOfElem+5]->setVertex(3,cc);

          theElement[noOfElem+6]->setVertex(0,ee);
          theElement[noOfElem+6]->setVertex(1,dd);
          theElement[noOfElem+6]->setVertex(2,cc);
          theElement[noOfElem+6]->setVertex(3,v[octTetVertex[tet][3]]);
        }
        noOfElem += 7;
        break;
      default:
        break;
    }
  }

  // Connect hierarchy
  child = new Element*[noOfElem];
  // Append element into grid level
  for (size_t nn = 0; nn < noOfElem; nn++)
  {
    theElement[nn]->setClosure();
    theElement[nn]->setValue(getValue());
    theElement[nn]->setParent(this);
    child[nn] = theElement[nn];
    gridLevel.appendTetrahedron(theElement[nn]);
  }

  // set no. of children
  chcount = static_cast<unsigned char>(noOfElem);
  // set irregular mark
  setIrregular();
  // connetc
  connect(gridLevel);
} // closure()



// Description
//  computes the edge refinement pattern
int
Octahedron::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;
  for (size_t i = 0; i < 12; i++) {
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
Octahedron::queryEdgeRefinement() const
{
  for (int i = 0; i < 12; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Octahedron::queryChildrenEdgeRefinement() const
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

// Description
//
int
Octahedron::getCommonFace(Element* theElem) const
{
  for (int face = 0; face < 8; face++)
    if (n[face] == theElem)
      return face;

  return -1;
}



// Description
//
void
Octahedron::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
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
Octahedron::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 6; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Octahedron::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Octahedron::getVertices(Vertex** vt) const
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
Octahedron::getEdge(int i) const
{
  return e[i];
}


// Description
//
void
Octahedron::getEdges(Edge** edg) const
{
  for (int i = 0; i < 12; i++)
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
Octahedron::getNeighbor(int face) const
{
  return n[face];
}


//=============================================================================
//  TOPOLOGY
//=============================================================================

//
//
int
Octahedron::type() const
{
  return GRD_OCTAHEDRON;
}


// Description
//
int
Octahedron::noOfVerticesAtFace(int fc) const
{
  return 3;
}


// Description
//
int
Octahedron::topologicalVertexAtFace(int fc,int vt) const
{
  return octVertexOfFace[fc][vt];
}

// Description
//
Vertex*
Octahedron::getVertexAtFace(int fc,int vt) const
{
  return (v[octVertexOfFace[fc][vt]]);
}

// Description
//
int
Octahedron::noOfEdgesAtFace(int fc) const
{
  return 3;
}

// Description
//
Vertex*
Octahedron::getVertexAtEdge(int ed,int vt) const
{
  return (v[octVertexOfEdge[ed][vt]]);
}

//
//
Edge*
Octahedron::getEdgeAtFace(int fc,int ed) const
{
  return e[octEdgeOfFace[fc][ed]];
}

// Description
//
void
Octahedron::connect(GridLevel& gridLevel)
{
  // Create List of objects to be connected
  std::list<Element*> nelt;
  std::list<Element*> nelc;

  // Collect all neighbors
  for (int nn = 0; nn < 8; nn++)
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
  for (int nn = 0; nn < 12; nn++)
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
} // connect


} // namespace grd
