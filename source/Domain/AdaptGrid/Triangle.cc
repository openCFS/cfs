/***************************************************************************
    File        : Triangle.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/


#include "Triangle.hh"
#include "GridLevel.hh"

#include "Utility.hh"

namespace grd {

// initialize memory management here
//Pool<Triangle> Triangle::pool;


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
  
  // no. of children
  // chcount = 0; <- set in abstract class Element

}


//=============================================================================
//  REFINEMENT
//=============================================================================


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
      else
      {
        Error("triangle can't find offset");
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
    for (vertex = 0; vertex < 3; vertex++) {
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
  //else
    //std::cout << "Error: 1 NULL neighbor" << std::endl;
    
  index = (4 - offset) % 3;
  child[1]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[1]);
  //else
    //std::cout << "Error: 2 NULL neighbor" << std::endl;

  index = (5 - offset) % 3;
  child[2]->setNeighbor(0,theChild[index]);
  if (theChild[index])
    theChild[index]->setNeighbor(theChildFace[index],child[2]);
  //else
    //std::cout << "Error: 3 NULL neighbor" << std::endl;

  child[3]->setNeighbor(0,theChild[3]);
  if (theChild[3])
    theChild[3]->setNeighbor(theChildFace[3],child[3]);
  //else
    //std::cout << "Error: 4 NULL neighbor" << std::endl;
    
  // set tag
  setRefined();
  
  // set no. of children
  chcount = 4;

  // set value
  for (i = 0; i < 4; i++) child[i]->setValue(getValue());
}



// Description
// all children are marked for coarsening
// remove children from element list, and
// set pointers to neighbors to NULL.
void
Triangle::coarsen(GridLevel& gridLevel)
{
  int nFace;
  Element* theNeighbor;

  for (size_t ch = 0; ch < 4; ch++)
  {
    // *- set neighbor's pointer
    theNeighbor = child[ch]->getNeighbor(0);
    if (theNeighbor) {
      nFace = theNeighbor->getCommonFace(child[ch]);
      theNeighbor->setNeighbor(nFace,0);
    }
    gridLevel.removeTriangle(child[ch]);
  }

  // Remove the children from memory
  for (size_t ch = 0; ch < 4; ch++)
    delete child[ch];

  // Delete the pointers
  delete [] child;
  // Initialize pointer
  child = 0;

  // set tag
  setRegular();
  // set no. of children
  chcount = 0;
} // corasen()


// Description
//
void
Triangle::close(GridLevel& gridLevel)
{
  size_t refPattern = 0;
  size_t refMask    = 1;
  size_t noOfElem;

  Vertex* aa;
  Vertex* bb;
  Vertex* cc;

  // Element buffer
  Element* theElement[4];


  // Get refinement pattern
  for (size_t i = 0; i < 3; i++)
  {
    if (e[i]->isRefined())
    { refPattern = (refPattern | refMask); }
    refMask = (refMask<<1);
  }

  switch(refPattern)
  {
    case 0:
      noOfElem = 0;
      break;
    case 1:
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;
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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;

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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;

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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;
      theElement[2] = new Triangle;

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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;
      theElement[2] = new Triangle;

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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;
      theElement[2] = new Triangle;

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
      theElement[0] = new Triangle;
      theElement[1] = new Triangle;
      theElement[2] = new Triangle;
      theElement[3] = new Triangle;

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

  // Connect hierarchy
  child = new Element*[noOfElem];

  // Append new elements to next level grid
  for (size_t nn = 0; nn < noOfElem; nn++)
  {
    theElement[nn]->setClosure();
    theElement[nn]->setValue(getValue());
    theElement[nn]->setParent(this);
    child[nn] = theElement[nn];
    gridLevel.appendTriangle(theElement[nn]);
  }
  
  // set no. of children
  chcount = static_cast<unsigned char>(noOfElem);
  // set irregular mark
  setIrregular();
  // set neighbors and create new edges
  connect(gridLevel);
}




// Description
//  Computes the edge refinement pattern
int
Triangle::getRefinementPattern() const
{
  int refMask = 1;
  int refPattern = 0;
  for (int i = 0; i < 3; i++)
  {
    if (e[i]->isRefined())
    { refPattern = (refPattern | refMask); }
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
inline bool
Triangle::queryEdgeRefinement() const
{
  for (int i = 0; i < 3; i++)
  {
    if (e[i]->isRefined())
      return true;
  }
  return false;
}


// Description
//
inline bool
Triangle::queryChildrenEdgeRefinement() const
{
  for (int i = 0; i < 3; i++)
  {
    if (e[i]->queryChildrenRefinement())
      return true;
  }
  return false;
}

// Description
// function for 3D elements
int
Triangle::getCommonFace(Element* theNeighbor) const
{
  if (n == theNeighbor)
    return 0;
  else
    return -1;
}



// Description
// common elements for tris as faces of tetras or prisms or pyramids
void
Triangle::getRefElements(int face,Vertex** theVertex,Edge** theEdge,Element** theChild,int* theChildFace) const
{
  if (face == 0) {
    // the order determines the orientation
    // towards the volume
    theVertex[0] = v[0];
    theVertex[1] = v[1];
    theVertex[2] = v[2];

    theEdge[0] = child[0]->getEdge(0);
    theEdge[1] = child[1]->getEdge(1);
    theEdge[2] = child[2]->getEdge(2);

    theChild[0] = child[0];
    theChild[1] = child[1];
    theChild[2] = child[2];
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
Triangle::setVertices(Vertex* const* vt)
{
  for (int i = 0; i < 3; i++)
    v[i] = vt[i];
}


// Description
//
Vertex*
Triangle::getVertex(int i) const
{
  return v[i];
}


// Description
//
void
Triangle::getVertices(Vertex** vt) const
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
Triangle::getEdge(int i) const
{
  return e[i];
}


// Description
//
void
Triangle::getEdges(Edge** edg) const
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
Triangle::getNeighbor(int face) const
{
  if (face == 0)
    return n;
  else
    return 0;
}

//=============================================================================
//  TOPOLOGY
//=============================================================================

//
//
int
Triangle::type() const
{
  return GRD_TRIANGLE;
}


// Description
//
int
Triangle::noOfVerticesAtFace(int fc) const
{
  return 3;
}

// Description
//
int
Triangle::topologicalVertexAtFace(int fc,int vt) const
{
  if (fc == 0) return vt;
  else return -1;
}


// Description
//
Vertex*
Triangle::getVertexAtFace(int fc,int vt) const
{
  if (fc == 0) return (v[vt]);
  else return 0;
  
}

// Description
//
int
Triangle::noOfEdgesAtFace(int fc) const
{
  return 3;
}

// Description
//  get vertex from edge by internal topology number
Vertex*
Triangle::getVertexAtEdge(int ed,int vt) const
{
  return (v[triVertexOfEdge[ed][vt]]);
}


// Description
//  get edge from face by toplogy number
Edge*
Triangle::getEdgeAtFace(int fc,int ed) const
{
  if (fc == 0) return e[ed];
  else return 0;
  
}


// Description
//
void
Triangle::connect(GridLevel& gridLevel)
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
    // Triangles have only one neighbor 
    Element* neighbor = child[ch]->getNeighbor(0);
    if (neighbor)
      continue;
    else
    {
      Vertex* vt[3];
      vt[0] = child[ch]->getVertex(0);
      vt[1] = child[ch]->getVertex(1);
      vt[2] = child[ch]->getVertex(2);
      
      // process all neighbors
      for (EI q = nelt.begin(); q != nelt.end(); ++q)
      {
        for (int nn = 0; nn < (*q)->getNoOfFaces(); nn++)
        {
          Vertex* vv[3];
          vv[0] = (*q)->getVertexAtFace(nn,0);
          vv[1] = (*q)->getVertexAtFace(nn,1);
          vv[2] = (*q)->getVertexAtFace(nn,2);
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
    for (int nn = 0; nn < 3; nn++)
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
  for (int nn = 0; nn < 3; nn++)
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
    for (int i = 0; i < 3; i++)
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


// // Private methods
// void
// Triangle::connect(GridLevel& gridLevel)
// {
//   // get all existing edges
//   std::list<Edge*> elt;
//   typedef std::list<Edge*>::iterator EI;
//   for (int nn = 0; nn < 3; nn++)
//   {
//     if (e[nn]->isRefined())
//     {
//       elt.push_back(e[nn]->getChild(0));
//       elt.push_back(e[nn]->getChild(1));
//     }
//     else
//     {
//       elt.push_back(e[nn]);
//     }
//   }
//   // connect this edges
//   for (int nn = 0; nn < static_cast<int>(chcount); nn++)
//   {
//     for (int i = 0; i < 3; i++)
//     {
//       Vertex* vt0 = child[nn]->getVertex(triVertexOfEdge[i][0]);
//       Vertex* vt1 = child[nn]->getVertex(triVertexOfEdge[i][1]);
//       for (EI p = elt.begin(); p != elt.end(); ++p)
//       {
//         Vertex* ev0 = (*p)->getVertex(0);
//         Vertex* ev1 = (*p)->getVertex(1);
//         if ((vt0 == ev0 && vt1 == ev1) || (vt1 == ev0 && vt0 == ev1))
//         {
//           child[nn]->setEdge(i,*p);
//           elt.erase(p);
//           break;
//         }
//         
//       }
//     }
//   }
//   
//   // if list size > 0
//   // we have a catastrophic error
//   if (elt.size() > 0)
//   {
//     std::cerr << "ERROR: connection of triangles failed, edge list has more than zero elements" << std::endl;
//     exit(1);
//   }
//   
//   // create new edges
//   for (int nn = 0; nn < static_cast<int>(chcount); nn++)
//   {
//     for (int i = 0; i < 3; i++)
//     {
//       Edge* edge = child[nn]->getEdge(i);
//       if (edge == 0)
//       {
//         Vertex* vt0 = child[nn]->getVertex(triVertexOfEdge[i][0]);
//         Vertex* vt1 = child[nn]->getVertex(triVertexOfEdge[i][1]);
//         bool flag = false;
//         for (EI p = elt.begin(); p != elt.end(); ++p)
//         {
//           Vertex* ev0 = (*p)->getVertex(0);
//           Vertex* ev1 = (*p)->getVertex(1);
//           if ((vt0 == ev0 && vt1 == ev1) || (vt1 == ev0 && vt0 == ev1))
//           {
//             child[nn]->setEdge(i,*p);
//             flag = true;
//             break;
//           }
//         } // for EI p
//         
//         // check if edge found
//         // if not, create a new edge
//         if (!flag)
//         {
//           Edge* tmp = new Edge(vt0,vt1);
//           child[nn]->setEdge(i,tmp);
//           elt.push_back(tmp);
//           gridLevel.appendEdge(tmp);
//         }
//       } // if (edge == 0)
//     } // for i
//   } // for nn chcount
//   
//   // clear list
//   elt.clear();
// }

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
