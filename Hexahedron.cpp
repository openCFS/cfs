// - C++ -
/***************************************************************************
    File        : Hexahedron.cpp
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Feb 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Hexahedron.h"
#include "GridLevel.h"

using namespace std;
namespace grd {

// Definition and initialisation of static class data
Pool<Hexahedron> Hexahedron::pool;


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

  c = 0;
}


// Description
// if all children marked for coarsening
bool
Hexahedron::queryChildrenMarks()
{
  int i;
  int counter = 0;

  if (!isRefined())
    return false;

  for (i = 0; i < 8; i++)
  {
    if (child[i]->isMarkedForCoarsening())
      counter++;
  }
  // set back children makrs
  if (counter == 0)
    return false;
  else if (counter > 0 && counter < 8)
  {
    for (i = 0; i < 8; i++)
    {
      if (child[i]->isMarkedForCoarsening())
        child[i]->setRegular();
    }
    return false;
  }
  else if (counter == 8)
  {
    // can't refine if a child's edge is refined
    bool flag = false;
    for (i = 0; i < 8; i++)
    {
      if (child[i]->queryEdgeRefinement()) {
        flag = true;
        break;
      }
    }
    if (flag)
    {
      for (i = 0; i < 8; i++)
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
  for (edge = 0; edge < 12; edge++) {
    if (e[edge]->isRefined())
    {
      e[edge]->getChildren(newVertex[hexVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
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
  if (c == 0)
  {
    c = new Vertex;
    for (i = 0; i < 3; i++) {
      for (vertex = 0; vertex < 8; vertex++) {
        (*c)[i] += (*v[vertex])[i];
      }
      (*c)[i] *= 0.125;
    }
    gridLevel.appendVertex(c);
  }
  newVertex[26] = c;

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
}


// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Hexahedron::coarsen(GridLevel& gridLevel)
{
  int ch;
  int face;
  int nFace;
  Element* theNeighbor;

  // *- children
  for (ch = 0; ch < 8; ch++) {
    // *- set neighbor's pointer
    for (face = 0; face < 4; face++) {
      theNeighbor = child[ch]->getNeighbor(face);
      if (theNeighbor) {
        nFace = theNeighbor->getCommonFace(child[ch]);
        theNeighbor->setNeighbor(nFace,0);
      }
    }
    gridLevel.removeHexahedron(child[ch]);
  }


  for (ch = 0; ch < 5; ch++)
    delete child[ch];
  child = 0;

  // set tag
  setRegular();

}


//
//
void
Hexahedron::close(int refP,ConformingClosure& closure)
{
  return;
}


// Description
//  close element and append vertex
//  into GridLevel vertex list
void
Hexahedron::close(GridLevel& gridLevel)
{
  Element::close();
  if (c == 0)
  {
    c = new Vertex;
    for (int i = 0; i < 3; i++) {
      for (int vertex = 0; vertex < 8; vertex++) {
        (*c)[i] += (*v[vertex])[i];
      }
      (*c)[i] *= 0.125;
    }
    gridLevel.appendVertex(c);
  }
}



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


// Description
//
void
Hexahedron::close(ConformingClosure& closure)
{
  int refPattern = getRefinementPattern();
  close(refPattern,closure);
}



// Description
//
bool
Hexahedron::queryChildrenRefinement()
{
  if(isRefined()) {
    for (int i = 0; i < 8; i++) {
      if (child[i]->isRefined())
        return true;
    }
  }

  return false;
}

// Description
//
bool
Hexahedron::queryEdgeRefinement()
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
Hexahedron::queryChildrenEdgeRefinement()
{
  int i;
  for (i = 0; i < 12; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  return false;
}


int
Hexahedron::getCommonFace(Element* theNeighbor)
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
Hexahedron::getRefElements(int face,Vertex** theVertex,
                           Edge** theEdge,Element** theChild,
                           int* theChildFace)
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
Hexahedron::getVertex(int i)
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
Hexahedron::getEdge(int i)
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
Hexahedron::getNeighbor(int face)
{
  return n[face];
}



//
//
void
Hexahedron::setVertices(Vertex** vt)
{
  for (int i = 0; i < 8; i++)
  {
    v[i] = vt[i];
  }
}

//
//
void
Hexahedron::getVertices(Vertex** vt)
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
Hexahedron::getEdges(Edge** edg)
{
  for (int i = 0; i < 12; i++)
  {
    edg[i] = e[i];
  }
}


// Topology
int
Hexahedron::type()
{
  return GRD_HEXAHEDRON;
}

int
Hexahedron::noOfVertexAtFace(int fc)
{
  return 4;
}

Vertex*
Hexahedron::getVertexAtFace(int fc,int vt)
{
  return v[hexVertexOfFace[fc][vt]];
}

Vertex*
Hexahedron::getVertexAtEdge(int edg,int vt)
{
  return v[hexVertexOfEdge[edg][vt]];
}

Edge*
Hexahedron::getEdgeAtFace(int fc,int edg)
{
  return e[hexEdgeOfFace[fc][edg]];
}



} // namespace grd