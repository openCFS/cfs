/***************************************************************************
    File        : Triangle.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Triangle.h"
#include "GridLevel.h"
#include "ConformingClosure.h"

namespace grd {

// initialize memory management here
Pool<Triangle> Triangle::pool;


// Description:
//
Triangle::Triangle( )
{
  int i;
  // Vertices and edges
  for (i = 0; i < 3; i++) {
    v[i] = 0;
    e[i] = 0;
  }
  // Pointer to volume element
  n = 0;

  // sboundary surface
  skeleton = 0;
}


// Description
// if all children are connected to volume elements
// only for 3D elements
bool
Triangle::queryChildrenMarks()
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
Triangle::refine(GridLevel& gridLevel)
{
  int i,ch;
  int vertCount,edgCount;
  int vertex,edge;
  int vertexIndex;

  int index;
  int offset;
  int nFace;

  int theChildFace[4];

  Vertex*  newVertex[6];
  Edge*    newEdge[9];
  Vertex*  theVertex[3];
  Edge*    theEdge[3];
  Element* theChild[4];

  // -* init indices
  offset = 0;
  index = -1;
  nFace = -1;

  // -* init counters
  vertCount = edgCount = 0;

  // -* 1. process parent vertices
  for (vertex = 0; vertex < 3; vertex++)
    newVertex[vertCount++] = v[vertex];

  // -* 2. process parent edges
  for (edge = 0; edge < 3; edge++) {
    if (e[edge]->isRefined()) {
      e[edge]->getChildren(newVertex[triVertexOfEdge[edge][0]],
                           &newEdge[edgCount],&newEdge[edgCount+1]);
      edgCount += 2;
      newVertex[vertCount++] = e[edge]->getMidPoint();
    }
    else {
      e[edge]->refine();
      e[edge]->getChildren(newVertex[triVertexOfEdge[edge][0]],
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
      else {
        std::cerr << "ERROR: triangle can't find offset\n";
        exit(1);
      }

      index = (3 - offset) % 3;
      newEdge[6] = theEdge[index];

      index = (4 - offset) % 3;
      newEdge[7] = theEdge[index];

      index = (5 - offset) % 3;
      newEdge[8] = theEdge[index];

    }
    else {
      Edge* e6 = new Edge;
      e6->setVertices(newVertex[5],newVertex[4]);
      Edge* e7 = new Edge;
      e7->setVertices(newVertex[3],newVertex[5]);
      Edge* e8 = new Edge;
      e8->setVertices(newVertex[4],newVertex[3]);

      newEdge[6] = e6;
      newEdge[7] = e7;
      newEdge[8] = e8;

      gridLevel.appendEdge(e6);
      gridLevel.appendEdge(e7);
      gridLevel.appendEdge(e8);

      // if volume not refined init theChild field
      for (i = 0; i < 4; i++)
              theChild[i] = 0;
    }
  }
  else {
    Edge* e6 = new Edge;
    e6->setVertices(newVertex[5],newVertex[4]);
    Edge* e7 = new Edge;
    e7->setVertices(newVertex[3],newVertex[5]);
    Edge* e8 = new Edge;
    e8->setVertices(newVertex[4],newVertex[3]);

    newEdge[6] = e6;
    newEdge[7] = e7;
    newEdge[8] = e8;

    gridLevel.appendEdge(e6);
    gridLevel.appendEdge(e7);
    gridLevel.appendEdge(e8);

    // if volume not refined init theCild field
    for (i = 0; i < 4; i++)
      theChild[i] = 0;
  }


  // create children
  child = new Element*[4];

  child[0] = new Triangle;
  child[1] = new Triangle;
  child[2] = new Triangle;
  child[3] = new Triangle;

  gridLevel.appendTriangle(child[0]);
  gridLevel.appendTriangle(child[1]);
  gridLevel.appendTriangle(child[2]);
  gridLevel.appendTriangle(child[3]);

  // set children parent
  child[0]->setParent(this);
  child[1]->setParent(this);
  child[2]->setParent(this);
  child[3]->setParent(this);


  // assemble children from new elements
  for (ch = 0; ch < 4; ch++) {
    // set vertices
    for (vertex = 0; vertex < 4; vertex++) {
      vertexIndex = triVertexOfTriChild[ch][vertex];
      child[ch]->setVertex(vertex,newVertex[vertexIndex]);
    }
    // set edges
    int edgeIndex;
    for (edge = 0; edge < 3; edge++) {
      edgeIndex = triEdgeOfTriChild[ch][edge];
      child[ch]->setEdge(edge,newEdge[edgeIndex]);
    }
  }

  // connect children with volume
  index = (3 - offset) % 3;
  child[0]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[0]);

  index = (4 - offset) % 3;
  child[1]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[1]);


  index = (5 - offset) % 3;
  child[2]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[2]);

  child[3]->setNeighbor(0,theChild[3]);
  if (theChild[3])
    theChild[3]->setNeighbor(theChildFace[3],child[3]);

  // set tag
  setRefined();

  // set value
  for (i = 0; i < 4; i++)
    child[i]->setValue(getValue());

  // set new vertex position
  if (skeleton)
  {
    setVertexPosition();
    // set reference to child vertex pos
    ((Triangle*)child[0])->skeleton = skeleton->child[0];
    ((Triangle*)child[1])->skeleton = skeleton->child[1];
    ((Triangle*)child[2])->skeleton = skeleton->child[2];
    ((Triangle*)child[3])->skeleton = skeleton->child[3];
  }
}



// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Triangle::coarsen(GridLevel& gridLevel)
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
    gridLevel.removeTriangle(child[ch]);
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
Triangle::close(int refPattern,ConformingClosure& closure)
{
  //int refPattern;
  int noOfElem;

  ConformingClosure::TRI theElement = closure.beginTriangle();
  Vertex* aa;
  Vertex* bb;
  Vertex* cc;

  // init element buffer
  closure.init();

  //refPattern = getRefinementPattern();

  switch(refPattern) {
    case 0:
      noOfElem = 0;
      break;
    case 1:
      aa = e[0]->getMidPoint();
      theElement[0]->setVertex(0,v[0]);
      theElement[0]->setVertex(1,v[1]);
      theElement[0]->setVertex(2,aa);

      theElement[1]->setVertex(0,v[0]);
      theElement[1]->setVertex(1,aa);
      theElement[1]->setVertex(2,v[2]);

      noOfElem = 2;
      break;
    case 2:
      aa = e[1]->getMidPoint();
      theElement[0]->setVertex(0,v[0]);
      theElement[0]->setVertex(1,v[1]);
      theElement[0]->setVertex(2,aa);

      theElement[1]->setVertex(0,aa);
      theElement[1]->setVertex(1,v[1]);
      theElement[1]->setVertex(2,v[2]);

      noOfElem = 2;
      break;
    case 4:
      aa = e[2]->getMidPoint();
      theElement[0]->setVertex(0,v[0]);
      theElement[0]->setVertex(1,aa);
      theElement[0]->setVertex(2,v[2]);

      theElement[1]->setVertex(0,aa);
      theElement[1]->setVertex(1,v[1]);
      theElement[1]->setVertex(2,v[2]);

      noOfElem = 2;
      break;
    case 3:
      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      if ((v[0]->getId()) > (v[1]->getId())) {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,aa);
        theElement[0]->setVertex(2,bb);

        theElement[1]->setVertex(0,v[0]);
        theElement[1]->setVertex(1,v[1]);
        theElement[1]->setVertex(2,aa);

        theElement[2]->setVertex(0,bb);
        theElement[2]->setVertex(1,aa);
        theElement[2]->setVertex(2,v[2]);
      }
      else {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,v[1]);
        theElement[0]->setVertex(2,bb);

        theElement[1]->setVertex(0,bb);
        theElement[1]->setVertex(1,v[1]);
        theElement[1]->setVertex(2,aa);

        theElement[2]->setVertex(0,bb);
        theElement[2]->setVertex(1,aa);
        theElement[2]->setVertex(2,v[2]);
      }
      noOfElem = 3;
      break;
    case 5:
      aa = e[0]->getMidPoint();
      bb = e[2]->getMidPoint();
      if ((v[0]->getId()) > (v[2]->getId())) {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,bb);
        theElement[0]->setVertex(2,aa);

        theElement[1]->setVertex(0,bb);
        theElement[1]->setVertex(1,v[1]);
        theElement[1]->setVertex(2,aa);

        theElement[2]->setVertex(0,v[0]);
        theElement[2]->setVertex(1,aa);
        theElement[2]->setVertex(2,v[2]);
      }
      else {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,bb);
        theElement[0]->setVertex(2,v[2]);

        theElement[1]->setVertex(0,bb);
        theElement[1]->setVertex(1,v[1]);
        theElement[1]->setVertex(2,aa);

        theElement[2]->setVertex(0,bb);
        theElement[2]->setVertex(1,aa);
        theElement[2]->setVertex(2,v[2]);
      }
      noOfElem = 3;
      break;
    case 6:
      aa = e[1]->getMidPoint();
      bb = e[2]->getMidPoint();
      if ((v[1]->getId()) > (v[2]->getId())) {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,bb);
        theElement[0]->setVertex(2,aa);

        theElement[1]->setVertex(0,bb);
        theElement[1]->setVertex(1,v[1]);
        theElement[1]->setVertex(2,aa);

        theElement[2]->setVertex(0,aa);
        theElement[2]->setVertex(1,v[1]);
        theElement[2]->setVertex(2,v[2]);
      }
      else {
        theElement[0]->setVertex(0,v[0]);
        theElement[0]->setVertex(1,bb);
        theElement[0]->setVertex(2,aa);

        theElement[1]->setVertex(0,aa);
        theElement[1]->setVertex(1,bb);
        theElement[1]->setVertex(2,v[2]);

        theElement[2]->setVertex(0,bb);
        theElement[2]->setVertex(1,v[1]);
        theElement[2]->setVertex(2,v[2]);
      }
      noOfElem = 3;
      break;
    case 7:
      aa = e[0]->getMidPoint();
      bb = e[1]->getMidPoint();
      cc = e[2]->getMidPoint();

      theElement[0]->setVertex(0,v[0]);
      theElement[0]->setVertex(1,cc);
      theElement[0]->setVertex(2,bb);

      theElement[1]->setVertex(0,cc);
      theElement[1]->setVertex(1,v[1]);
      theElement[1]->setVertex(2,aa);

      theElement[2]->setVertex(0,bb);
      theElement[2]->setVertex(1,aa);
      theElement[2]->setVertex(2,v[2]);

      theElement[3]->setVertex(0,aa);
      theElement[3]->setVertex(1,bb);
      theElement[3]->setVertex(2,cc);

      noOfElem = 4;
      break;
    default:
      noOfElem = 0;
      break;
  }

  closure.setNoOfTriangles(noOfElem);

}


//
//
void
Triangle::setVertexPosition()
{
  // for each child
  for (int ch = 0; ch < 4; ch++)
  {
    Triangle* tri = (Triangle*) child[ch];
    TriSkeleton* s = skeleton->child[ch];
    // for each vertex
    for (int vt = 0; vt < 3; vt++)
    {
      double* skPos = s->pos[vt];
      double* pos   = (tri->getVertex(vt))->getPosition();
      // set 3D coordinates
      for (int n = 0; n < 3; n++)
      {
        pos[n] = skPos[n];
      }
    }
  }
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
void
Triangle::close(ConformingClosure& closure)
{
  int refPattern = getRefinementPattern();
  close(refPattern,closure);
}


// Description
//
bool
Triangle::queryChildrenRefinement()
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
Triangle::queryEdgeRefinement()
{
  for (int i = 0; i < 3; i++)
    if (e[i]->isRefined())
      return true;

  return false;
}


// Description
//
bool
Triangle::queryChildrenEdgeRefinement()
{
  for (int i = 0; i < 3; i++)
    if (e[i]->queryChildrenRefinement())
      return true;

  return false;
}


// include element topology information
#include "TriangleTopology.h"



// Description
// function for 3D elements
int
Triangle::getCommonFace(Element* theNeighbor)
{
  if (n == theNeighbor)
    return 0;
  else
    return -1;
}



// Description
// common elements for tris as faces of tetras or prisms or pyramids
void
Triangle::getRefElements(int face,Vertex** theVertex,
                         Edge** theEdge,Element** theChild,
                         int* theChildFace)
{
  if (face == 0) {
    // the order determines the orientation
    // towards the volume
    theVertex[0] = v[0];
    theVertex[1] = v[2];
    theVertex[2] = v[1];

    theEdge[0] = child[0]->getEdge(0);
    theEdge[1] = child[2]->getEdge(2);
    theEdge[2] = child[1]->getEdge(1);

    theChild[0] = child[0];
    theChild[1] = child[2];
    theChild[2] = child[1];
    theChild[3] = child[3];

    theChildFace[0] = 0;
    theChildFace[1] = 0;
    theChildFace[2] = 0;
    theChildFace[3] = 0;

  }
}


// Description
//
void
Triangle::setVertex(int i,Vertex* vt)
{
  v[i] = vt;
}


// Description
//
void
Triangle::setVertices(Vertex** vt)
{
  for (int i = 0; i < 3; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Triangle::getVertex(int i)
{
  return v[i];
}


// Description
//
void
Triangle::getVertices(Vertex** vt)
{
  for (int i = 0; i < 3; i++)
    vt[i] = v[i];
}


// Description
//
void
Triangle::setEdge(int i,Edge* edg)
{
  e[i] = edg;
}


// Description
//
void
Triangle::setEdges(Edge** edg)
{
  for (int i = 0; i < 3; i++)
    e[i] = edg[i];
}



// Description
//
Edge*
Triangle::getEdge(int i)
{
  return e[i];
}


// Description
//
void
Triangle::getEdges(Edge** edg)
{
  for (int i = 0; i < 3; i++)
    edg[i] = e[i];
}


// Description
//
void
Triangle::setNeighbor(int face,Element* theNeighbor)
{
  if (face == 0)
    n = theNeighbor;
}


// Description
//
Element*
Triangle::getNeighbor(int face)
{
  if (face == 0)
    return n;
  else
    return 0;
}

// Topology

//
//
int
Triangle::type()
{
  return GRD_TRIANGLE;
}


// Description
//
int
Triangle::noOfVertexAtFace(int fc)
{
  return 2;
}

// Description
//
Vertex*
Triangle::getVertexAtFace(int fc,int vt)
{
  return (v[triVertexOfFace[fc][vt]]);
}


// Description
//  get vertex from edge by internal topology number
Vertex*
Triangle::getVertexAtEdge(int ed,int vt)
{
  return (v[triVertexOfEdge[ed][vt]]);
}


// Description
//  get edge from face by toplogy number
Edge*
Triangle::getEdgeAtFace(int fc,int ed)
{
  return e[triEdgeOfFace[fc][ed]];
}


//=============================================================================
//  GEOMETRY
//=============================================================================

//// Description
////
//void
//Triangle::getGlobalCoords(const double* lc,double* coords)
//{
//  double* pos[3];
//  // set local coordinates
//  double l1 = lc[0];
//  double l2 = lc[1];
//  double l0 = 1.0 - l1 - l2;
//  // get vertex position
//  pos[0] = v[0]->getPosition();
//  pos[1] = v[1]->getPosition();
//  pos[2] = v[2]->getPosition();
//  // set global coordinates
//  coords[0] = l0*pos[0][0] + l1*pos[1][0] + l2*pos[2][0];
//  coords[1] = l0*pos[0][1] + l1*pos[1][1] + l2*pos[2][1];
//  coords[2] = l0*pos[0][2] + l1*pos[1][2] + l2*pos[2][2];
//}
//
//
//// Description
////
//void
//Triangle::getJacobian(Jacobian& Jac)
//{
//  Jac(0,0) = (*v[1])[0] - (*v[0])[0];
//  Jac(0,1) = (*v[2])[0] - (*v[0])[0];
//  Jac(0,2) = 0.0;
//
//  Jac(1,0) = (*v[1])[1] - (*v[0])[1];
//  Jac(1,1) = (*v[2])[1] - (*v[0])[1];
//  Jac(1,2) = 0.0;
//
//  Jac(2,0) = 0.0;
//  Jac(2,1) = 0.0;
//  Jac(2,2) = 1.0;
//
//}



} // namespace grd
