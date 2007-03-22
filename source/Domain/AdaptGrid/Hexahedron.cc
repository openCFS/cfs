// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : Hexahedron.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Hexahedron.hh"
#include "Tetrahedron.hh"
#include "Pyramid.hh"
#include "GridLevel.hh"


namespace grd {

// Definition and initialisation of static class data
//Pool<Hexahedron> Hexahedron::pool;


// Description:
//
Hexahedron::Hexahedron( )
{
  int i;
  for (i = 0; i < 8; i++)
    v[i] = 0;

  for (i = 0; i < 12; i++)
    e[i] = 0;

  for (i = 0; i < 6; i++)
    n[i] = 0;

 // set no. of children
 // chcount = 0; <- set in abstract class Element
}

//=============================================================================
//  REFINEMENT
//=============================================================================

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
Hexahedron::refine(GridLevel& gridLevel)
{
  int ch;
  int i;
  int vertCount,edgCount;
  int vertex,edge,face,neighbor;
  int vertexIndex,edgeIndex,neighIndex,childIndex;
  int nVertex,nFace;
  int neighFace[32],theChildFace[4];

  Vertex*  newVertex[27];
  Edge*    newEdge[54];
  Element* newNeighbor[32];
  Vertex*  theVertex[5]; // include vertex from face
  Edge*    theEdge[4];
  Element* theChild[4];

  // -* init counters
  vertCount = edgCount = 0;

  // -* 1. process parent vertices
  for (vertex = 0; vertex < 8; vertex++)
    newVertex[vertCount++] = v[vertex];

  // -* 2. process parent edges
  for (edge = 0; edge < 12; edge++)
  {
    if (e[edge]->isRefined())
    {
      e[edge]->getChildren(newVertex[hexVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else
    {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[hexVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      newVertex[vertCount] = e[edge]->getMidPoint();
      gridLevel.appendVertex(newVertex[vertCount]);
      gridLevel.appendEdge(newEdge[edgCount]);
      gridLevel.appendEdge(newEdge[edgCount+1]);

      edgCount += 2;
      vertCount++;
    }
  }


  // -* 3. process parent faces
  for (face = 0; face < 6; face++)
  {
    // if neighbor regular refined
    if (n[face])
    {
      if (n[face]->isRefined())
      {
        // get neighbor's common face
        nFace = n[face]->getCommonFace(this);
        n[face]->getRefElements(nFace,theVertex,theEdge,theChild,theChildFace);

        // for all vertices of neighbor's face
        int index = 3;
        vertex    = 3;
        for (nVertex = 0; nVertex < 4; nVertex++)
        {
          // find common vertex of parent and neighbor
          vertex = ((index--)%4);
          vertexIndex = hexVertexOfFace[face][vertex];
          if (index < 0) index = 3;
          while(v[vertexIndex] != theVertex[nVertex])
          {
            vertex      = ((index--)%4);
            vertexIndex = hexVertexOfFace[face][vertex];
            if (index < 0) index = 3;
          }
          // set child edge
          edgeIndex = hexEdgeFromParentFace[face][vertex];
          newEdge[edgeIndex] = theEdge[nVertex];
          // set child neighbor
          neighIndex = hexNeighborOfParentFace[face][vertex];
          newNeighbor[neighIndex] = theChild[nVertex];
          neighFace[neighIndex]   = theChildFace[nVertex];
        }
        // set last vertex from face middle
        vertexIndex = hexVertexFromParentFace[face];
        newVertex[vertexIndex] = theVertex[4];
      } // if n[face]->isRefined()
      else {
        // neighbor is not refined
        // create new vertex and edges on face
        // set children neighbors at face to NULL

        Vertex* tmpVertex = new Vertex;
        // set position
        for (i = 0; i < 3; i++)
        {
          for (vertex = 0; vertex < 4; vertex++)
          {
            vertexIndex = hexVertexOfFace[face][vertex];
            (*tmpVertex)[i] += (*v[vertexIndex])[i];
          }
          (*tmpVertex)[i] *= 0.25;
        }
        vertexIndex = hexVertexFromParentFace[face];
        newVertex[vertexIndex] = tmpVertex;
        gridLevel.appendVertex(tmpVertex);

        int v1,v2;
        for (edge = 0; edge < 4; edge++)
        {
          Edge* tmpEdge = new Edge;
          v1 = hexVertexOfEdgeOfParentFace[face][edge][0];
          v2 = hexVertexOfEdgeOfParentFace[face][edge][1];
          tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
          newEdge[hexEdgeFromParentFace[face][edge]] = tmpEdge;
          gridLevel.appendEdge(tmpEdge);
        }
        for (neighbor = 0; neighbor < 4; neighbor++)
        {
          neighIndex = hexNeighborOfParentFace[face][neighbor];
          newNeighbor[neighIndex] = 0;
        }
      }  // else
    } // if n[face]
    else {
      // neighbor is not set
      // create new vertex and edges on face
      // set children neighbors at face to NULL

      Vertex* tmpVertex = new Vertex;
      // set position
      for (i = 0; i < 3; i++) {
        for (vertex = 0; vertex < 4; vertex++) {
          vertexIndex = hexVertexOfFace[face][vertex];
          (*tmpVertex)[i] += (*v[vertexIndex])[i];
        }
        (*tmpVertex)[i] *= 0.25;
      }
      vertexIndex = hexVertexFromParentFace[face];
      newVertex[vertexIndex] = tmpVertex;
      gridLevel.appendVertex(tmpVertex);

      int v1,v2;
      for (edge = 0; edge < 4; edge++)
      {
        Edge* tmpEdge = new Edge;
        v1 = hexVertexOfEdgeOfParentFace[face][edge][0];
        v2 = hexVertexOfEdgeOfParentFace[face][edge][1];
        tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
        newEdge[hexEdgeFromParentFace[face][edge]] = tmpEdge;
        gridLevel.appendEdge(tmpEdge);
      }
      for (neighbor = 0; neighbor < 4; neighbor++)
      {
        neighIndex = hexNeighborOfParentFace[face][neighbor];
        newNeighbor[neighIndex] = 0;
      }
    }
  }


  // *- create elements from volume
  // *- barycenter
  Vertex* bc = new Vertex;
  for (i = 0; i < 3; i++)
  {
    for (vertex = 0; vertex < 8; vertex++)
    {
      (*bc)[i] += (*v[vertex])[i];
    }
    (*bc)[i] *= 0.125;
  }
  gridLevel.appendVertex(bc);

  newVertex[26] = bc;

  // -* create interior edges from parent face vertices
  for (edge = 0; edge < 6; edge++)
  {
    Edge* tmpEdge = new Edge;
    int v1,v2;
    v1 = hexVertexOfIntEdge[edge][0];
    v2 = hexVertexOfIntEdge[edge][1];
    tmpEdge->setVertices(newVertex[v1],newVertex[v2]);
    newEdge[48+edge] = tmpEdge;
    gridLevel.appendEdge(tmpEdge);
  }


  // create children
  child = new Element*[8];

  child[0]  = (Element* ) new Hexahedron;
  child[1]  = (Element* ) new Hexahedron;
  child[2]  = (Element* ) new Hexahedron;
  child[3]  = (Element* ) new Hexahedron;
  child[4]  = (Element* ) new Hexahedron;
  child[5]  = (Element* ) new Hexahedron;
  child[6]  = (Element* ) new Hexahedron;
  child[7]  = (Element* ) new Hexahedron;

  gridLevel.appendHexahedron(child[0]);
  gridLevel.appendHexahedron(child[1]);
  gridLevel.appendHexahedron(child[2]);
  gridLevel.appendHexahedron(child[3]);
  gridLevel.appendHexahedron(child[4]);
  gridLevel.appendHexahedron(child[5]);
  gridLevel.appendHexahedron(child[6]);
  gridLevel.appendHexahedron(child[7]);

  // set newNeighbor for interior connections
  newNeighbor[24] = child[0];  // hex 0
  newNeighbor[25] = child[1];  // hex 1
  newNeighbor[26] = child[2];  // hex 2
  newNeighbor[27] = child[3];  // hex 3
  newNeighbor[28] = child[4];  // hex 4
  newNeighbor[29] = child[5];  // hex 5
  newNeighbor[30] = child[6];  // hex 6
  newNeighbor[31] = child[7];  // hex 7

  // set children parent
  child[0]->setParent(this);
  child[1]->setParent(this);
  child[2]->setParent(this);
  child[3]->setParent(this);
  child[4]->setParent(this);
  child[5]->setParent(this);
  child[6]->setParent(this);
  child[7]->setParent(this);


  // assemble hexa children from new elements
  for (ch = 0; ch < 8; ch++)
  {
    // set vertices
    int vertexIndex;
    for (vertex = 0; vertex < 8; vertex++)
    {
      vertexIndex = hexVertexOfHexaChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 12; edge++)
    {
      edgeIndex = hexEdgeOfChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
    //   connect child
    int neighIndex;
    for (face = 0; face < 6; face++)
    {
      neighIndex = hexHexaChildNeighOfFace[ch][face];
      child[ch]->setNeighbor(face,newNeighbor[neighIndex]);
    }
  } // for loop


  // connect external neighbor with children
  for (neighbor = 0; neighbor < 24; neighbor++)
  {
    if (newNeighbor[neighbor])
    {
      childIndex = hexChildNeighOfNeighbor[neighbor];
      newNeighbor[neighbor]->setNeighbor(neighFace[neighbor],child[childIndex]);
    }
  }

  // set tag
  setRefined();
  // set no. of children
  chcount = 8;
  
  // set element value
  for (ch = 0; ch < 8; ch++) child[ch]->setValue(getValue());
} // refine()


// Description
// all children are marked for coarsening
// remove children from element list, and
// set neighbor's pointers to NULL.
void
Hexahedron::coarsen(GridLevel& gridLevel)
{
  int ch;
  int face;
  int nFace;
  Element* theNeighbor;

  // *- children
  for (ch = 0; ch < 8; ch++)
  {
    // *- set neighbor's pointer
    for (face = 0; face < 6; face++)
    {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor)
      {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeHexahedron(child[ch]);
  }

  // Delete the children from memory
  for (ch = 0; ch < 8; ch++)
    delete child[ch];
  // Delete the pointers
  delete [] child;
  // Initialize pointer
  child = 0;

  // set tag
  setRegular();
  // set no. of children
  chcount = 0;

} // coarse()


// Description
//
void
Hexahedron::close(GridLevel& gridLevel)
{
  // Inidices
  size_t face;
  size_t index;
  size_t noOfTets = 0;
  size_t noOfPyrs = 0;
  size_t noOfHexs = 0;

  // Masks
  size_t ALL_FACES_REFINED = 63;
  size_t faceMask[6]  = {1,2,4,8,16,32};
  size_t edgeMask[12] = {1,2,4,8,16,32,64,256,512,1048,4096,8192};

  // Auxiliary vertices
  Vertex* vertex[27];

  // Pointers to closure elements
  Element* tets[24];
  Element* pyrs[24];
  Element* hexs[8];

  // Get vertices
  for (index = 0; index < 8; index++) vertex[index] = v[index];

  // Compute Edge refinement pattern
  size_t edgeRefPattern = 0;
  for (index = 0; index < 12; index++)
  {
    Edge* edge = e[index];
    if (e[index]->isRefined())
    {
      edgeRefPattern = (edgeRefPattern | edgeMask[index]);
      vertex[8+index] = edge->getMidPoint();
    }
    else
    {
      vertex[8+index] = 0;
    }
  }

  // Compute Face refinement pattern
  size_t faceRefPattern = 0;
  for (face = 0; face < 6; face++)
  {
    if (n[face])
    {
      if (n[face]->isRefined())
      {
        vertex[20+face] = ((Hexahedron*) n[face])->getFaceVertex(this);
        faceRefPattern = (faceRefPattern | faceMask[face]);
      }
    }
  }

  // Set barycenter
  vertex[26] = new Vertex;
  for (index = 0; index < 3; index++)
  {
    for (size_t vt = 0; vt < 8; vt++)
    {
        (*vertex[26])[index] += (*v[vt])[index];
    }
    (*vertex[26])[index] *= 0.125;
  }
  // store vertex in list at gird level
  gridLevel.appendVertex(vertex[26]); 
  
  
  // If all faces refined, subdivide element
  // into hexahedra
  if (faceRefPattern == ALL_FACES_REFINED)
  {
    for (size_t ch = 0; ch < 8; ch++)
    {
      hexs[ch] = new Hexahedron;
      // set vertices
      for (size_t vt = 0; vt < 8; vt++)
      {
        size_t vertexIndex = hexVertexOfHexaChild[ch][vt];
        hexs[ch]->setVertex(vt,vertex[vertexIndex]);
      }
    }
  }
  else // Refinement: process face wise
  {
    for (face = 0; face < 6; face++)
    {
      // Face is refined, create pyramids
      if ((faceRefPattern & faceMask[face]))
      {
        // Create new pyramids
        pyrs[noOfPyrs]   = new Pyramid;
        pyrs[noOfPyrs+1] = new Pyramid;
        pyrs[noOfPyrs+2] = new Pyramid;
        pyrs[noOfPyrs+3] = new Pyramid;

        for (size_t pyr = 0; pyr < 4; pyr++)
        {
          // Pyramid has five vertices
          for (index = 0; index < 5; index++)
          {
            pyrs[noOfPyrs+pyr]->setVertex(index,vertex[hexSlicePyr[face][pyr][index]]);
          }
        }
        // Increse number or pyramids
        noOfPyrs += 4;
      } // Face is refined, create pyramids
      else // Face is not refined, refine quad and create tets and pyrs
      {
        // Refinement case
        switch(hexSlicePattern[face][edgeRefPattern])
        {
          case 0: // face and edges not refined, create pyr
            pyrs[noOfPyrs] = new Pyramid;
            pyrs[noOfPyrs]->setVertex(0,vertex[hexVertexOfFace[face][0]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexVertexOfFace[face][3]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexVertexOfFace[face][2]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexVertexOfFace[face][1]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfPyrs++;
            break;
          case 1:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            noOfTets += 3;
            break;
          case 2:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            noOfTets += 3;
            break;
          case 3:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;
            tets[noOfTets+3] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            tets[noOfTets+3]->setVertex(0,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+3]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+3]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+3]->setVertex(3,vertex[26]);

            noOfTets += 4;
            break;
          case 4:
            tets[noOfTets] = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            noOfTets += 3;
            break;
          case 5:
            pyrs[noOfPyrs]   = new Pyramid;
            pyrs[noOfPyrs+1] = new Pyramid;

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][4]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][3]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            pyrs[noOfPyrs+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            pyrs[noOfPyrs+1]->setVertex(1,vertex[hexSliceTetra[face][1]]);
            pyrs[noOfPyrs+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            pyrs[noOfPyrs+1]->setVertex(3,vertex[hexSliceTetra[face][6]]);
            pyrs[noOfPyrs+1]->setVertex(4,vertex[26]);

            noOfPyrs += 2;
            break;
          case 6:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;
            tets[noOfTets+3] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            tets[noOfTets+3]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets+3]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+3]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+3]->setVertex(3,vertex[26]);

            noOfTets += 4;
            break;
          case 7:
            pyrs[noOfPyrs] = new Pyramid;

            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][4]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][3]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfTets += 3;
            noOfPyrs += 1;
            break;
          case 8:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            noOfTets += 3;
            break;
          case 9:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;
            tets[noOfTets+3] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            tets[noOfTets+3]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+3]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+3]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+3]->setVertex(3,vertex[26]);

            noOfTets += 4;
            break;
          case 10:
            pyrs[noOfPyrs]   = new Pyramid;
            pyrs[noOfPyrs+1] = new Pyramid;

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][1]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][7]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            pyrs[noOfPyrs+1]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            pyrs[noOfPyrs+1]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            pyrs[noOfPyrs+1]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            pyrs[noOfPyrs+1]->setVertex(3,vertex[hexSliceTetra[face][3]]);
            pyrs[noOfPyrs+1]->setVertex(4,vertex[26]);

            noOfPyrs += 2;
            break;
          case 11:
            pyrs[noOfPyrs] = new Pyramid;

            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][3]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfTets += 3;
            noOfPyrs += 1;
            break;
          case 12:
            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;
            tets[noOfTets+3] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            tets[noOfTets+3]->setVertex(0,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+3]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+3]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+3]->setVertex(3,vertex[26]);

            noOfTets += 4;
            break;
          case 13:
            pyrs[noOfPyrs] = new Pyramid;

            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][1]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][6]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfTets += 3;
            noOfPyrs += 1;
            break;
          case 14:
            pyrs[noOfPyrs] = new Pyramid;

            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][1]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][7]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfTets += 2;
            noOfPyrs += 1;
            break;
          case 15:
            pyrs[noOfPyrs] = new Pyramid;

            tets[noOfTets]   = new Tetrahedron;
            tets[noOfTets+1] = new Tetrahedron;
            tets[noOfTets+2] = new Tetrahedron;
            tets[noOfTets+3] = new Tetrahedron;

            tets[noOfTets]->setVertex(0,vertex[hexSliceTetra[face][0]]);
            tets[noOfTets]->setVertex(1,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets]->setVertex(3,vertex[26]);

            tets[noOfTets+1]->setVertex(0,vertex[hexSliceTetra[face][7]]);
            tets[noOfTets+1]->setVertex(1,vertex[hexSliceTetra[face][3]]);
            tets[noOfTets+1]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+1]->setVertex(3,vertex[26]);

            tets[noOfTets+2]->setVertex(0,vertex[hexSliceTetra[face][6]]);
            tets[noOfTets+2]->setVertex(1,vertex[hexSliceTetra[face][2]]);
            tets[noOfTets+2]->setVertex(2,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+2]->setVertex(3,vertex[26]);

            tets[noOfTets+3]->setVertex(0,vertex[hexSliceTetra[face][5]]);
            tets[noOfTets+3]->setVertex(1,vertex[hexSliceTetra[face][1]]);
            tets[noOfTets+3]->setVertex(2,vertex[hexSliceTetra[face][4]]);
            tets[noOfTets+3]->setVertex(3,vertex[26]);

            pyrs[noOfPyrs]->setVertex(0,vertex[hexSliceTetra[face][4]]);
            pyrs[noOfPyrs]->setVertex(1,vertex[hexSliceTetra[face][5]]);
            pyrs[noOfPyrs]->setVertex(2,vertex[hexSliceTetra[face][6]]);
            pyrs[noOfPyrs]->setVertex(3,vertex[hexSliceTetra[face][7]]);
            pyrs[noOfPyrs]->setVertex(4,vertex[26]);

            noOfTets += 4;
            noOfTets += 1;
            break;
          default:
            break;
        } // switch(hexSlicePattern[face][edgeRefPattern])
      } // Face is not refined, refine quad and create tets and pyrs
    } // Main loop, index is face
  } // if all face refined

  // Connect hierarchy
  size_t noOfChildren = noOfTets + noOfPyrs + noOfHexs;
  child = new Element*[noOfChildren];

  // Add elements to list
  size_t count = 0;
  for (size_t nn = 0; nn < noOfTets; nn++)
  {
    tets[nn]->setClosure();
    tets[nn]->setValue(getValue());
    tets[nn]->setParent(this);
    child[count++] = tets[nn];
    gridLevel.appendTetrahedron(tets[nn]);
  }
  for (size_t nn = 0; nn < noOfPyrs; nn++)
  {
    pyrs[nn]->setClosure();
    pyrs[nn]->setValue(getValue());
    pyrs[nn]->setParent(this);
    child[count++] = pyrs[nn];
    gridLevel.appendPyramid(pyrs[nn]);
  }
  for (size_t nn = 0; nn < noOfHexs; nn++)
  {
    hexs[nn]->setClosure();
    hexs[nn]->setValue(getValue());
    hexs[nn]->setParent(this);
    child[count++] = hexs[nn];
    gridLevel.appendHexahedron(hexs[nn]);
  }
  
  // set no. of children
  chcount = static_cast<unsigned char>(noOfTets + noOfPyrs + noOfHexs);
  // set irregular mark
  setIrregular();
} // closure()



//// Description
//// trilinear mapping
//void
//Hexahedron::getGlobalCoords(const double* lc,double* coord)
//{
//  int i;
//  double f0,f1,f2,f3,f4,f5,f6,f7;
//  double* x[8];
//
//  for(i = 0; i < 8; i++)
//    x[i] = v[i]->getPosition();
//
//  for(i = 0; i < 3; i++) {
//    f0 = (-x[0][i] - x[6][i] - x[2][i] - x[4][i] + x[5][i] + x[1][i] + x[3][i] + x[7][i]);
//    f1 = x[0][i] - x[5][i] + x[6][i] - x[7][i];
//    f2 = x[0][i] + x[2][i] - x[1][i] - x[7][i];
//    f3 = - x[0][i] + x[7][i];
//    f4 = x[0][i] + x[4][i] - x[5][i] - x[1][i];
//    f5 = - x[0][i] + x[5][i];
//    f6 = -x[0][i] + x[1][i];
//    f7 = x[0][i];
//
//    coord[i] = (((f0*lc[0] + f1)*lc[1] + f2*lc[0] + f3)*lc[2] + (f4*lc[0] + f5)*lc[1] + f6* lc[0] + f7);
//  }
//
//}



// ____________________________________________________________
//
//
//           functions
//
// ____________________________________________________________
//

//=============================================================================
//  TOPOLOGY
//=============================================================================

//
// Description
//
int
Hexahedron::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;

  for (size_t i = 0; i < 12; i++)
  {
    if (e[i]->isRefined())
    {
      refPattern = (refPattern | refMask);
    }
    refMask = (refMask<<1);
  }
  return refPattern;
}



// Description
//
bool
Hexahedron::queryEdgeRefinement() const
{
  int i;
  for (i = 0; i < 12; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Hexahedron::queryChildrenEdgeRefinement() const
{
  int i;
  for (i = 0; i < 12; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  return false;
}


int
Hexahedron::getCommonFace(Element* theNeighbor) const
{
  int face;
  for (face = 0; face < 8; face++) {
    if (n[face] == theNeighbor)
      return face;
  }
  return -1;
}



// Description
//
void
Hexahedron::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
{
  int vertex,index, childIndex;

  for (vertex = 0; vertex < 4; vertex++) {
    index = hexVertexOfFace[face][vertex];
    theVertex[vertex] = v[index];

    childIndex = hexChildOfParentFace[face][vertex];
    theChild[vertex] = child[childIndex];
    theChildFace[vertex] = hexChildFaceOfParentFace[face][vertex];

    index = hexChildEdgeOfParentFace[face][vertex];
    theEdge[vertex] = child[childIndex]->getEdge(index);
  }

  childIndex = hexChildOfParentFace[face][0];
  theVertex[4] = child[childIndex]->getVertex(hexChildVertexOfParentFace[face]);

}


//


// Description
//
void
Hexahedron::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
Vertex*
Hexahedron::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Hexahedron::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
Edge*
Hexahedron::getEdge(int i) const
{
  return e[i];

}


// Description
//
void
Hexahedron::setNeighbor(int face,Element* theNeighbor)
{
  n[face] = theNeighbor;
}


// Description
//
Element*
Hexahedron::getNeighbor(int face) const
{
  return n[face];
}



//
//
void
Hexahedron::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 8; i++)
  {
    v[i] = vt[i];
  }
}

//
//
void
Hexahedron::getVertices(Vertex** vt) const
{
  for (int i = 0; i < 8; i++)
  {
    vt[i] = v[i];
  }
}


//
//
void
Hexahedron::setEdges(Edge** edg)
{
  for (int i = 0; i < 12; i++)
  {
    e[i] = edg[i];
  }
}


//
//
void
Hexahedron::getEdges(Edge** edg) const
{
  for (int i = 0; i < 12; i++)
  {
    edg[i] = e[i];
  }
}


//=============================================================================
//  TOPOLOGY
//=============================================================================

// Description
//
int
Hexahedron::type() const
{
  return GRD_HEXAHEDRON;
}

// Description
//
int
Hexahedron::noOfVerticesAtFace(int fc) const
{
  return 4;
}


// Description
//
int
Hexahedron::topologicalVertexAtFace(int fc,int vt) const
{
  return hexVertexOfFace[fc][vt];
}


// Description
//
Vertex*
Hexahedron::getVertexAtFace(int fc,int vt) const
{
  return v[hexVertexOfFace[fc][vt]];
}

// Description
//
int
Hexahedron::noOfEdgesAtFace(int fc) const
{
  return 4;
}

// Description
//
Vertex*
Hexahedron::getVertexAtEdge(int edg,int vt) const
{
  return v[hexVertexOfEdge[edg][vt]];
}

Edge*
Hexahedron::getEdgeAtFace(int fc,int edg) const
{
  return e[hexEdgeOfFace[fc][edg]];
}


} // namespace grd
