// - C++ -
/***************************************************************************
    File        : TriangleTopology.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TRIANGLETOPOLOGY_H
#define TRIANGLETOPOLOGY_H

namespace grd {

//
// Topological description of triangle
// __________________________________________________________________
//
// Description:
/*
       v2
       / \
    e1/   \e0
     /     \
    /       \
 v0 ---------v1
       e2
*/



// __________________________________________________________________
//
// Triangle refinement
// __________________________________________________________________

// refinement into 4 triangle: regular ref.of tris
//
// Vertices of children: 6 Vertices
//      from parent	from parent edges
// 	0	0	0	3
// 	1	1	1	4
// 	2	2	2	5
//
// Edges of children: from parent edges 3x2
// 	from parent edg		form parent interior
//	0	0,1		pi	6,7,8
//	1	2,3
//	2	4,5




//_____________________________________________________________________
//
// Topology

// vertex indices at each face. Indices are in that order to
// get the outer normal, triangles has one face (itself)
const int triVertexOfFace[3] = {0,1,2};

// vertex indices at each edge
// don't change order; is for refinement
const int triVertexOfEdge[3][2] = {{1,2}, {2,0}, {0,1}};

// edge indices at each face, outer normal order
// triangles has one face (itself)
const int triEdgeOfFace[3] = {0,1,2};



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
        ^
       /2\
    e1/___\e0
     /\ 3/ \
    /0 \/ 1 \
 v0 ---------v1
       e2
*/

// child at parent face by vertex
//const int triChildOfParentFace[3][2] = {{1,2},{2,0},{0,1}};

// child face at parent face
//const int triChildFaceOfParentFace[3][2] = {{0,0},{1,1},{2,2}};


// child edge of child at face
//const int triChildEdgeOfParentFace[3][2] = {{0,0},{1,1},{2,2}};


// the new vertices, edges and neighbors are incrementally numbred
// assemble children from this elements

// Vertices of tri children
const int triVertexOfTriChild[4][3] = {{0,5,4},{5,1,3},
                                        {4,3,2},{3,4,5}};

// Edges of tetra children
const int triEdgeOfTriChild[4][3] = {{6,3,4},{0,7,5},
                                      {1,2,8},{6,7,8}};



// _____________________________________________________________________
//
// information necessary to get vertices, edges and children
// neighbors from refined neighbor element. Take vertex of face defintion
// into account for orientation and consistency


// neighbor over parent face number according to outer normal at
// parent face and to parent face order
//const int triNeighborOfParentFace[3][2] = {{0,1},{2,3},{4,5}};


// tetra child external connection or connection over parent face
//const int triNeighOfFaceOfTriChild[4][3] = {{9,3,4},{0,9,5},
//                                             {1,2,9},{6,7,8}};


// child neighbor of external neighbor
//const int triChildNeighOfNeighbor[6] = {1,2,2,0,0,1};


// _____________________________________________________________________
//
// Topology for getting neighbor connection from volume element
// This is used for boundary refinement from volume refinement

// vertex from tetrahedron faces
//const int triVertexFromTetraFace[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};

// edge from tetrahedron faces
//const int triEdgeFromTetraFace[4][3] = {{3,4,1},{3,0,5},{4,5,2},{1,2,0}};

// vertex from octahedron faces
//const int triVertexFromOctaFace[8][3] = {{0,1,2},{0,2,3},{0,3,4},{0,4,1},
//                                          {5,2,1},{5,3,2},{5,4,3},{5,1,4}};

// edge from octahedron faces
//const int triEdgeFromOctaFace[8][3] = {{4,1,0},{5,2,1},{6,3,2},{7,0,3},
//                                        {4,8,9},{5,9,10},{6,10,11},{7,11,8}};

} // namespace grd
#endif // TRIANGLETOPOLOGY_H
