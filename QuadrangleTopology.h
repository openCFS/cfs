// - C++ -
/***************************************************************************
    File        : QuadrangleTopology.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Feb 7 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef QUADRANGLETOPOLOGY_H
#define QUADRANGLETOPOLOGY_H

namespace grd {

//
// Topological description of quadrangle
// __________________________________________________________________
//
// Description:
/*
    v3  e2   v2
     --------
		|        |
  e3|        |e1
    |        |
    |        |
     --------
    v0  e0  v1

*/



// __________________________________________________________________
//
// Tetrahedron refinement
// __________________________________________________________________

// refinement into 4 quadrangle: regular ref.of quads
//




//_____________________________________________________________________
//
// Topology

// vertex indices at each face. Indices are in that order to
// get the outer normal
const int quaVertexOfFace[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};

// vertex indices at each edge
// don't change order; is for refinement
const int quaVertexOfEdge[4][2] = {{0,1}, {1,2}, {2,3}, {3,0}};

// edge indices at each face, outer normal order
const int quaEdgeOfFace[4][1] = {{0},{1},{2},{3}};


//_____________________________________________________________________
//
// Refinement
//
// 1. refinement: 4 elements
// 2. new vertices: 6
// 3. new edges:    9
// 4. neighbors over parent element faces: 6
// 5. interior neighbors: 4
// 6. total number of neighbors: 10
//
//
// Numbering of children
/*
   v3     e2     v2
    --------------
   |      |       |
   |      |       |
   | ch3  | ch2   |
e3	------|------   e1
   |      |       |
   |      |       |
   | ch0  | ch1   |
    --------------
   v0     e0      v1
*/

// child at parent face by vertex
//const int quaChildOfParentFace[4][2] = {{0,1},{1,2},{2,3}, {3,0}};

// child face at parent face
//const int quaChildFaceOfParentFace[4][2] = {{0,0},{1,1},{2,2}, {3,3}};


// child edge of child at face
//const int quaChildEdgeOfParentFace[4][2] = {{0,0},{1,1},{2,2}, {4,4}};


// the new vertices, edges and neighbors are incrementally numbred
// assemble children from this elements
// Vertices of qua children
const int quaVertexOfQuaChild[4][4] = {{0,4,8,7},{4,1,5,8},
                                        {8,5,2,6},{7,8,6,3}};

// Edges of quad children
const int quaEdgeOfQuaChild[4][4] = {{0,8,11,7},{1,2,9,8},
                                      {9,3,4,10},{11,10,5,6}};


// _____________________________________________________________________
//
// information necessary to get vertices, edges and children
// neighbors from refined neighbor element. Take vertex of face defintion
// into account for orientation and consistency


// neighbor over parent face number according to outer normal at
// parent face and to parent face order
//const int quaNeighborOfParentFace[3][2] = {{0,1},{2,3},{4,5}};


// tetra child external connection or connection over parent face
//const int quaNeighOfFaceOfQuaChild[4][3] = {{9,3,4},{0,9,5},
//                                             {1,2,9},{6,7,8}};


// child neighbor of external neighbor
//const int quaChildNeighOfNeighbor[6] = {1,2,2,0,0,1};


// _____________________________________________________________________
//
// Topology for getting neighbor connection from volume element
// This is used for boundary refinement from volume refinement

// vertex from tetrahedron faces
//const int quaVertexFromQuadrFace[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

// edge from tetrahedron faces
//const int triEdgeFromQuadrFace[4][3] = {{3,4,1},{3,0,5},{4,5,2},{1,2,0}};

// vertex from octahedron faces
//const int triVertexFromFace[8][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,1},
//				      {5,2,1},{5,3,2},{5,4,3},{5,1,4}};

// edge from octahedron faces
//const int triEdgeFromOctaFace[8][3] = {{4,1,0},{5,2,1},{6,3,2},{7,0,3},
//				    {4,8,9},{5,9,10},{6,10,11},{7,11,8}};

} // namespace grd
#endif // QUADRANGLETOPOLOGY_H
