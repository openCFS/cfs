// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/***************************************************************************
    File        : Pyramid.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Nov 29 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Pyramid.hh"
#include "Tetrahedron.hh"
#include "GridLevel.hh"

#include "Utility.hh"

namespace grd {

// initialize memory management here
//Pool<Pyramid> Pyramid::pool;


// Description:
//
Pyramid::Pyramid( )
{
  int i;
  for (i = 0; i < 5; i++)
    v[i] = 0;

  for (i = 0; i < 8; i++)
    e[i] = 0;

  for (i = 0; i < 5; i++)
    n[i] = 0;

  // set no. of children
  // chcount = 0; <- set in abstract class Element
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
Pyramid::refine(GridLevel& gridLevel)
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
      e[edge]->getChildren(newVertex[pyrVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[pyrVertexOfEdge[edge][0]],
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
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 6; j++) {
      (*bc)[i] += (*(v[j]))[i];
    }
    (*(bc))[i] *= ONE_SIXTH;
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
          vertexIndex = pyrVertexOfFace[face][vertex];
          if (index < 0) index = 2;
          while(v[vertexIndex] != theVertex[nVertex]) {
            vertex      = ((index--)%3);
            vertexIndex = pyrVertexOfFace[face][vertex];
            if (index < 0) index = 2;
          }
          // set child edge
          edgeIndex = pyrEdgeFromParentFace[face][vertex];
          newEdge[edgeIndex] = theEdge[nVertex];
          // set child neighbor
          neighIndex = pyrNeighborOfParentFace[face][vertex];
          newNeighbor[neighIndex] = theChild[nVertex];
          neighFace[neighIndex]   = theChildFace[nVertex];
        }
        //         set last neighbor from face middle
        neighIndex = pyrNeighborOfParentFace[face][3];
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
          v1 = pyrVertexOfEdgeOfParentFace[face][edge][0];
          v2 = pyrVertexOfEdgeOfParentFace[face][edge][1];
          tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
          newEdge[pyrEdgeFromParentFace[face][edge]] = tmpEdge;
          gridLevel.appendEdge(tmpEdge);
        }
        for (neighbor = 0; neighbor < 4; neighbor++) {
          neighIndex = pyrNeighborOfParentFace[face][neighbor];
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
        v1 = pyrVertexOfEdgeOfParentFace[face][edge][0];
        v2 = pyrVertexOfEdgeOfParentFace[face][edge][1];
        tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
        newEdge[pyrEdgeFromParentFace[face][edge]] = tmpEdge;
        gridLevel.appendEdge(tmpEdge);
      }
      for (neighbor = 0; neighbor < 4; neighbor++) {
        neighIndex = pyrNeighborOfParentFace[face][neighbor];
        newNeighbor[neighIndex] = 0;
      }
    }
  }


  // -* 4. create interior edges from parent edge to barycenter
  for (edge = 0; edge < 12; edge++) {
    Edge* tmpEdge = new Edge;
    int v1,v2;
    v1 = pyrVertexOfIntEdgeFromParentEdge[edge][0];
    v2 = pyrVertexOfIntEdgeFromParentEdge[edge][1];
    edgeIndex = pyrInteriorEdgeFromParentEdge[edge];
    tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
    newEdge[edgeIndex] = tmpEdge;
    gridLevel.appendEdge(tmpEdge);
  }


  // create children
  child = new Element*[14];

  child[0]  = new Pyramid;
  child[1]  = new Pyramid;
  child[2]  = new Pyramid;
  child[3]  = new Pyramid;
  child[4]  = new Pyramid;
  child[5]  = new Pyramid;
  child[6]  = new Tetrahedron;
  child[7]  = new Tetrahedron;
  child[8]  = new Tetrahedron;
  child[9]  = new Tetrahedron;

  gridLevel.appendPyramid(child[0]);
  gridLevel.appendPyramid(child[1]);
  gridLevel.appendPyramid(child[2]);
  gridLevel.appendPyramid(child[3]);
  gridLevel.appendPyramid(child[4]);
  gridLevel.appendPyramid(child[5]);
  gridLevel.appendTetrahedron(child[6]);
  gridLevel.appendTetrahedron(child[7]);
  gridLevel.appendTetrahedron(child[8]);
  gridLevel.appendTetrahedron(child[9]);

  // set newNeighbor for interior connections
  newNeighbor[20] = child[0];  // pyr 0
  newNeighbor[21] = child[1];  // pyr 1
  newNeighbor[22] = child[2];  // pyr 2
  newNeighbor[23] = child[3];  // pyr 3
  newNeighbor[24] = child[4];  // pyr 4
  newNeighbor[25] = child[5];  // pyr 5
  newNeighbor[26] = child[6];  // tet 0
  newNeighbor[27] = child[7];  // tet 1
  newNeighbor[28] = child[8];  // tet 2
  newNeighbor[29] = child[9];  // tet 3

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

  // assemble pyra children from new elements
  for (ch = 0; ch < 6; ch++) {
    // set vertices
    int vertexIndex;
    for (vertex = 0; vertex < 5; vertex++) {
      vertexIndex = pyrVertexOfPyramidChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 8; edge++) {
      edgeIndex = pyrEdgeOfPyramidChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
    //         connect child
    int neighIndex;
    for (face = 0; face < 5; face++) {
      neighIndex = pyrPyramidChildNeighOfFace[ch][face];
      child[ch]->setNeighbor(face,newNeighbor[neighIndex]);

    }
  }

  // assemble tetra children
  for (ch = 0; ch < 4; ch++) {
    for (vertex = 0; vertex < 4; vertex++) {
      vertexIndex = pyrVertexOfTetraChild[ch][vertex];
      child[ch+6]->setVertex(vertex,newVertex[vertexIndex]);
    }
    for (edge = 0; edge < 6; edge++) {
      edgeIndex = pyrEdgeOfTetraChild[ch][edge];
      child[ch+6]->setEdge(edge,newEdge[edgeIndex]);
    }
    for (face = 0; face < 4; face++) {
      neighIndex = pyrTetraChildNeighborOfFace[ch][face];
      child[ch+6]->setNeighbor(face,newNeighbor[neighIndex]);
    }
  }


  // connect external neighbor with children
  for (neighbor = 0; neighbor < 20; neighbor++) {
    if (newNeighbor[neighbor])
    {
      childIndex = pyrChildNeighOfNeighbor[neighbor];
      newNeighbor[neighbor]->setNeighbor(neighFace[neighbor],child[childIndex]);
    }
  }

  // set value
  for (ch = 0; ch < 10; ch++) child[ch]->setValue(getValue());

  // set tag
  setRefined();
  // set no. of children
  chcount = 10;

}



// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Pyramid::coarsen(GridLevel& gridLevel)
{
  int ch;
  int face;
  int nFace;
  Element* theNeighbor;

  // *- children, 6 pyramids
  for (ch = 0; ch < 6; ch++)
  {
    // *- set neighbor's pointer
    for (face = 0; face < 5; face++)
    {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor) {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeHexahedron(child[ch]);
  }

  // *- children, 4 tetrahedra
  for (ch = 6; ch < 10; ch++)
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
    gridLevel.removeHexahedron(child[ch]);
  }

  // Delete the children from memory
  for (ch = 0; ch < 10; ch++)
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
// Compute the conoforming closure
// Store the children element in the
// next grid level
void
Pyramid::close(GridLevel& gridLevel)
{
  int refPattern = 0;
  int refMask    = 1;
  size_t noOfTet = 0;
  size_t noOfPyr = 0;


  // Grid elements
  Element* tet[20];
  Element* pyr[20];
  Vertex* aa;
  Vertex* bb;
  Vertex* cc;
  Vertex* dd;
  Vertex* ee;
  Vertex* ff;
  Vertex* gg;
  Vertex* hh;
  
  // Declare tag and mask for case 15
  size_t tag15 = 0;
  size_t mask15[4] = {1,2,4,8};

  // Get edge refinement pattern
  for (size_t i = 0; i < 8; i++)
  {
    if (e[i]->isRefined()) refPattern = (refPattern | refMask);
    refMask = (refMask<<1);
  }

  // Set case and rotation
  size_t rot = pyrIsometryCases[refPattern][0];
  size_t ref = pyrIsometryCases[refPattern][1];

  // Set vertices and edges
  Vertex* vt[5];
  Edge*   eg[8];
  for (size_t n = 0; n < 5; n++)  vt[n] = v[pyrIsometryVertex[rot][n]];
  for (size_t n = 0; n < 8; n++)  eg[n] = e[pyrIsometryEdge[rot][n]];

  // Consider the isometry
  switch (ref)
  {
    case 1:
      pyr[0] = new Pyramid;
      tet[0] = new Tetrahedron;
      tet[1] = new Tetrahedron;
      // Get vertex from refined edges
      aa = eg[0]->getMidPoint();
      // Set pyramids
      pyr[0]->setVertex(0,aa);
      pyr[0]->setVertex(1,vt[1]);
      pyr[0]->setVertex(2,vt[2]);
      pyr[0]->setVertex(3,vt[3]);
      pyr[0]->setVertex(4,vt[4]);
      // Set tetrahedra
      tet[0]->setVertex(0,aa);
      tet[0]->setVertex(1,vt[2]);
      tet[0]->setVertex(2,vt[3]);
      tet[0]->setVertex(3,vt[0]);

      tet[1]->setVertex(0,aa);
      tet[1]->setVertex(1,vt[3]);
      tet[1]->setVertex(2,vt[4]);
      tet[1]->setVertex(3,vt[0]);
      // Set new number of elements
      noOfPyr = 1;
      noOfTet = 2;
      break;
    case 3:
      pyr[0] = new Pyramid;
      tet[0] = new Tetrahedron;
      tet[1] = new Tetrahedron;
      tet[2] = new Tetrahedron;
      // Get vertices from refined edges
      aa = eg[0]->getMidPoint();
      bb = eg[1]->getMidPoint();
      // Cases
      if (vt[1]->getId() > vt[2]->getId())
      {
        pyr[0]->setVertex(0,bb);
        pyr[0]->setVertex(1,vt[1]);
        pyr[0]->setVertex(2,vt[2]);
        pyr[0]->setVertex(3,vt[3]);
        pyr[0]->setVertex(4,vt[4]);

        tet[0]->setVertex(0,aa);
        tet[0]->setVertex(0,bb);
        tet[0]->setVertex(0,vt[4]);
        tet[0]->setVertex(0,vt[0]);

        tet[1]->setVertex(0,bb);
        tet[1]->setVertex(0,vt[3]);
        tet[1]->setVertex(0,vt[4]);
        tet[1]->setVertex(0,vt[0]);

        tet[2]->setVertex(0,vt[1]);
        tet[2]->setVertex(0,bb);
        tet[2]->setVertex(0,vt[4]);
        tet[2]->setVertex(0,aa);
      }
      else
      {
        pyr[0]->setVertex(0,aa);
        pyr[0]->setVertex(1,vt[1]);
        pyr[0]->setVertex(2,vt[2]);
        pyr[0]->setVertex(3,vt[3]);
        pyr[0]->setVertex(4,vt[4]);

        tet[0]->setVertex(0,aa);
        tet[0]->setVertex(0,bb);
        tet[0]->setVertex(0,vt[3]);
        tet[0]->setVertex(0,vt[0]);

        tet[1]->setVertex(0,bb);
        tet[1]->setVertex(0,vt[3]);
        tet[1]->setVertex(0,vt[4]);
        tet[1]->setVertex(0,vt[0]);

        tet[2]->setVertex(0,vt[3]);
        tet[2]->setVertex(0,aa);
        tet[2]->setVertex(0,vt[2]);
        tet[2]->setVertex(0,bb);
      }
      noOfPyr = 1;
      noOfTet = 3;
      break;
    case 5:
      pyr[0] = new Pyramid;
      tet[0] = new Tetrahedron;
      tet[1] = new Tetrahedron;
      tet[2] = new Tetrahedron;
      tet[3] = new Tetrahedron;

      aa = eg[0]->getMidPoint();
      bb = eg[2]->getMidPoint();

      pyr[0]->setVertex(0,bb);
      pyr[0]->setVertex(1,vt[1]);
      pyr[0]->setVertex(2,vt[2]);
      pyr[0]->setVertex(3,vt[3]);
      pyr[0]->setVertex(4,vt[4]);

      tet[0]->setVertex(0,aa);
      tet[0]->setVertex(1,vt[2]);
      tet[0]->setVertex(2,bb);
      tet[0]->setVertex(3,vt[0]);

      tet[1]->setVertex(0,aa);
      tet[1]->setVertex(1,bb);
      tet[1]->setVertex(2,vt[4]);
      tet[1]->setVertex(3,vt[0]);

      tet[2]->setVertex(0,vt[1]);
      tet[2]->setVertex(1,vt[2]);
      tet[2]->setVertex(2,bb);
      tet[2]->setVertex(3,aa);

      tet[3]->setVertex(0,vt[1]);
      tet[3]->setVertex(1,bb);
      tet[3]->setVertex(2,vt[4]);
      tet[3]->setVertex(3,aa);

      noOfPyr = 1;
      noOfTet = 4;
      break;
    case 7:
      aa = eg[0]->getMidPoint();
      bb = eg[2]->getMidPoint();
      if (v[1]->getId() > v[2]->getId())
      {
        if (v[2]->getId() > v[3]->getId())
        {
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;

          pyr[0]->setVertex(0,cc);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,cc);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,aa);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,vt[4]);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,vt[1]);
          tet[2]->setVertex(1,bb);
          tet[2]->setVertex(2,cc);
          tet[2]->setVertex(3,aa);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,cc);
          tet[3]->setVertex(2,vt[4]);
          tet[3]->setVertex(3,aa);

          tet[4]->setVertex(0,vt[1]);
          tet[4]->setVertex(1,vt[2]);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,bb);

          noOfTet = 5;
          noOfPyr = 1;
        }
        else
        {
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;

          pyr[0]->setVertex(0,bb);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,vt[4]);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,vt[4]);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,bb);
          tet[2]->setVertex(1,vt[3]);
          tet[2]->setVertex(2,vt[4]);
          tet[2]->setVertex(3,cc);

          tet[3]->setVertex(0,bb);
          tet[3]->setVertex(1,vt[4]);
          tet[3]->setVertex(2,vt[1]);
          tet[3]->setVertex(3,aa);

          noOfTet = 4;
          noOfPyr = 1;
        }
      }
      else
      { 
        if (v[2]->getId() > v[3]->getId())
        {
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;

          pyr[0]->setVertex(0,cc);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,cc);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,aa);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,vt[4]);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,vt[1]);
          tet[2]->setVertex(1,cc);
          tet[2]->setVertex(2,vt[4]);
          tet[2]->setVertex(3,aa);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,vt[2]);
          tet[3]->setVertex(2,cc);
          tet[3]->setVertex(3,aa);

          tet[4]->setVertex(0,aa);
          tet[4]->setVertex(1,vt[2]);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,bb);

          noOfTet = 5;
          noOfPyr = 1;
        }
        else
        {
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;

          pyr[0]->setVertex(0,aa);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,cc);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,aa);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,vt[4]);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,vt[3]);
          tet[2]->setVertex(2,vt[4]);
          tet[2]->setVertex(3,cc);

          tet[3]->setVertex(0,aa);
          tet[3]->setVertex(1,bb);
          tet[3]->setVertex(2,vt[3]);
          tet[3]->setVertex(3,cc);

          tet[4]->setVertex(0,aa);
          tet[4]->setVertex(1,vt[2]);
          tet[4]->setVertex(2,vt[3]);
          tet[4]->setVertex(3,bb);

          noOfTet = 5;
          noOfPyr = 1;
        }
      }
      break;
    case 15:
      aa = eg[0]->getMidPoint();
      bb = eg[1]->getMidPoint();
      cc = eg[2]->getMidPoint();
      dd = eg[3]->getMidPoint();

      if (vt[0]->getId() > vt[1]->getId()) tag15 |= mask15[0];
      if (vt[1]->getId() > vt[2]->getId()) tag15 |= mask15[1];
      if (vt[2]->getId() > vt[3]->getId()) tag15 |= mask15[2];
      if (vt[0]->getId() > vt[3]->getId()) tag15 |= mask15[3];
      switch(tag15)
      {
        case 15:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,dd);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,bb);
          tet[3]->setVertex(2,vt[2]);
          tet[3]->setVertex(3,dd);
          
          tet[4]->setVertex(0,bb);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[2]);
          
          tet[5]->setVertex(0,vt[2]);
          tet[5]->setVertex(1,cc);
          tet[5]->setVertex(2,vt[3]);
          tet[5]->setVertex(3,dd);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 13:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,bb);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,bb);
          tet[3]->setVertex(2,vt[4]);
          tet[3]->setVertex(3,dd);

          tet[4]->setVertex(0,vt[3]);
          tet[4]->setVertex(1,vt[4]);
          tet[4]->setVertex(2,bb);
          tet[4]->setVertex(3,dd);
          
          tet[5]->setVertex(0,bb);
          tet[5]->setVertex(1,dd);
          tet[5]->setVertex(2,cc);
          tet[5]->setVertex(3,vt[3]);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 11:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,cc);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,bb);
          tet[3]->setVertex(2,cc);
          tet[3]->setVertex(3,dd);

          tet[4]->setVertex(0,vt[1]);
          tet[4]->setVertex(1,cc);
          tet[4]->setVertex(2,vt[4]);
          tet[4]->setVertex(3,dd);

          tet[5]->setVertex(0,vt[1]);
          tet[5]->setVertex(1,bb);
          tet[5]->setVertex(2,vt[2]);
          tet[5]->setVertex(3,cc);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 9:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,bb);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,vt[3]);
          tet[3]->setVertex(2,dd);
          tet[3]->setVertex(3,bb);

          tet[4]->setVertex(0,vt[4]);
          tet[4]->setVertex(1,cc);
          tet[4]->setVertex(2,dd);
          tet[4]->setVertex(3,bb);

          tet[5]->setVertex(0,vt[4]);
          tet[5]->setVertex(1,vt[3]);
          tet[5]->setVertex(2,cc);
          tet[5]->setVertex(3,bb);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 14:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,dd);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,vt[1]);
          tet[2]->setVertex(1,aa);
          tet[2]->setVertex(2,vt[2]);
          tet[2]->setVertex(3,dd);

          tet[3]->setVertex(0,aa);
          tet[3]->setVertex(1,dd);
          tet[3]->setVertex(2,bb);
          tet[3]->setVertex(3,vt[2]);

          tet[4]->setVertex(0,bb);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[2]);

          tet[5]->setVertex(0,vt[1]);
          tet[5]->setVertex(1,cc);
          tet[5]->setVertex(2,vt[2]);
          tet[5]->setVertex(3,dd);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 3:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,cc);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,cc);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,aa);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,cc);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,aa);
          tet[3]->setVertex(1,cc);
          tet[3]->setVertex(2,vt[1]);
          tet[3]->setVertex(3,vt[4]);

          tet[4]->setVertex(0,aa);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[4]);

          tet[5]->setVertex(0,vt[1]);
          tet[5]->setVertex(1,bb);
          tet[5]->setVertex(2,vt[2]);
          tet[5]->setVertex(3,cc);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 10:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,cc);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,cc);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,aa);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,cc);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[2]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,aa);
          tet[3]->setVertex(2,vt[2]);
          tet[3]->setVertex(3,cc);

          tet[4]->setVertex(0,aa);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[1]);

          tet[5]->setVertex(0,vt[1]);
          tet[5]->setVertex(1,dd);
          tet[5]->setVertex(2,vt[2]);
          tet[5]->setVertex(3,cc);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 1:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,bb);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,dd);
          tet[3]->setVertex(2,bb);
          tet[3]->setVertex(3,vt[4]);

          tet[4]->setVertex(0,bb);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[4]);

          tet[5]->setVertex(0,vt[3]);
          tet[5]->setVertex(1,cc);
          tet[5]->setVertex(2,vt[4]);
          tet[5]->setVertex(3,bb);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        case 12:
          pyr[0] = new Pyramid;
          tet[0] = new Tetrahedron;
          tet[1] = new Tetrahedron;
          tet[2] = new Tetrahedron;
          tet[3] = new Tetrahedron;
          tet[4] = new Tetrahedron;
          tet[5] = new Tetrahedron;

          pyr[0]->setVertex(0,bb);
          pyr[0]->setVertex(1,vt[1]);
          pyr[0]->setVertex(2,vt[2]);
          pyr[0]->setVertex(3,vt[3]);
          pyr[0]->setVertex(4,vt[4]);

          tet[0]->setVertex(0,aa);
          tet[0]->setVertex(1,bb);
          tet[0]->setVertex(2,dd);
          tet[0]->setVertex(3,vt[0]);

          tet[1]->setVertex(0,bb);
          tet[1]->setVertex(1,cc);
          tet[1]->setVertex(2,dd);
          tet[1]->setVertex(3,vt[0]);

          tet[2]->setVertex(0,aa);
          tet[2]->setVertex(1,dd);
          tet[2]->setVertex(2,bb);
          tet[2]->setVertex(3,vt[1]);

          tet[3]->setVertex(0,vt[1]);
          tet[3]->setVertex(1,dd);
          tet[3]->setVertex(2,bb);
          tet[3]->setVertex(3,vt[4]);

          tet[4]->setVertex(0,bb);
          tet[4]->setVertex(1,dd);
          tet[4]->setVertex(2,cc);
          tet[4]->setVertex(3,vt[4]);

          tet[5]->setVertex(0,vt[3]);
          tet[5]->setVertex(1,cc);
          tet[5]->setVertex(2,vt[4]);
          tet[5]->setVertex(3,bb);

          noOfPyr = 1;
          noOfTet = 6;
          break;
        default:
          break;
      } // switch tag15
    default:
      break;
  }

  // Connect hierarchy
  size_t noOfChildren = noOfTet + noOfPyr;
  child = new Element*[noOfChildren];

  // Add elements to list
  size_t count = 0;
  for (size_t nn = 0; nn < noOfTet; nn++)
  {
    tet[nn]->setClosure();
    tet[nn]->setValue(getValue());
    tet[nn]->setParent(this);
    child[count++] = tet[nn];
    gridLevel.appendTetrahedron(tet[nn]);
  }
  for (size_t nn = 0; nn < noOfPyr; nn++)
  {
    pyr[nn]->setClosure();
    pyr[nn]->setValue(getValue());
    pyr[nn]->setParent(this);
    child[count++] = pyr[nn];
    gridLevel.appendPyramid(pyr[nn]);
  }
  
  // set no. of children
  chcount = static_cast<unsigned char>(noOfTet + noOfPyr);
  // set irregular mark
  setIrregular();
}



// Description
//  computes the edge refinemt pattern
int
Pyramid::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;

  for (size_t i = 0; i < 8; i++)
  {
    if (e[i]->isRefined())
    {
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
Pyramid::queryEdgeRefinement() const
{
  for (int i = 0; i < 8; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Pyramid::queryChildrenEdgeRefinement() const
{
  int i;
  for (i = 0; i < 8; i++) 
  {
    if (e[i]->queryChildrenRefinement())
      return true;
  }

  // check if neighbor's child is refined
  for (i = 0; i < 5; i++) 
  {
    if (n[i]) 
    {
      if (n[i]->queryChildrenRefinement())
        return true;
    }
  }
  return false;
}


int
Pyramid::getCommonFace(Element* theElem) const
{
  for (int face = 0; face < 5; face++)
    if (n[face] == theElem)
      return face;

  return -1;
}



// Description
//
void
Pyramid::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
{
  int vertex,index, childIndex;

  if (face < 4)
  {
    for (vertex = 0; vertex < 3; vertex++) {
      index = pyrVertexOfFace[face][vertex];
      theVertex[vertex] = v[index];

      childIndex = pyrChildOfParentFace[face][vertex];
      theChild[vertex] = child[childIndex];
      theChildFace[vertex] = pyrChildFaceOfParentFace[face][vertex];

      index = pyrChildEdgeOfParentFace[face][vertex];
      theEdge[vertex] = child[childIndex]->getEdge(index);
    }

    childIndex = pyrChildOfParentFace[face][3];
    theChild[3] = child[childIndex];
    theChildFace[3] = pyrChildFaceOfParentFace[face][3];
  }
  else
  {
    for (vertex = 0; vertex < 4; vertex++) {
      index = pyrVertexOfFace[face][vertex];
      theVertex[vertex] = v[index];

      childIndex = pyrChildOfParentFace[face][vertex];
      theChild[vertex] = child[childIndex];
      theChildFace[vertex] = pyrChildFaceOfParentFace[face][vertex];

      index = pyrChildEdgeOfParentFace[face][vertex];
      theEdge[vertex] = child[childIndex]->getEdge(index);
    }

    // get Vertex at midpoint of quadrilateral face
    theVertex[4] = child[1]->getVertex(3);
  }
}



// Description
//
void
Pyramid::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
void
Pyramid::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 5; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Pyramid::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Pyramid::getVertices(Vertex** vt) const
{
  for (int i = 0; i < 5; i++)
    vt[i] = v[i];
}


// Description
//
void
Pyramid::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
void
Pyramid::setEdges(Edge** edg)
{
  for (int i = 0; i < 8; i++)
    e[i] = edg[i];
}


// Description
//
Edge*
Pyramid::getEdge(int i) const
{
  return e[i];
}


// Description
//
void
Pyramid::getEdges(Edge** edg) const
{
  for (int i = 0; i < 8; i++)
    edg[i] = e[i];
}


// Description
//
void
Pyramid::setNeighbor(int face,Element* theNeighbor)
{
  n[face] = theNeighbor;
}


// Description
//
Element*
Pyramid::getNeighbor(int face) const
{
  return n[face];
}


//=============================================================================
//  TOPOLOGY
//=============================================================================

//
//
int
Pyramid::type() const
{
  return GRD_PYRAMID;
}



// Description
//
int
Pyramid::noOfVerticesAtFace(int fc) const
{
  return pyrNoOfVertexOfFace[fc];
}

// Description
//
int
Pyramid::topologicalVertexAtFace(int fc,int vt) const
{
  return pyrVertexOfFace[fc][vt];
}

// Description
//
Vertex*
Pyramid::getVertexAtFace(int fc,int vt) const
{
  return (v[pyrVertexOfFace[fc][vt]]);
}

// Description
//
int
Pyramid::noOfEdgesAtFace(int fc) const
{
  return pyrNoOfEdgeOfFace[fc];
}

// Description
//
Vertex*
Pyramid::getVertexAtEdge(int ed,int vt) const
{
  return (v[pyrVertexOfEdge[ed][vt]]);
}

//
//
Edge*
Pyramid::getEdgeAtFace(int fc,int ed) const
{
  return e[pyrEdgeOfFace[fc][ed]];
}



} // namespace grd
