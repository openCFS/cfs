// - C++ -
/***************************************************************************
    File        : TetrahedronTopology.h
    Description :

 ---------------------------------------------------------------------------
    Begin       : Fri Dec 7 2001
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef TETRAHEDRONTOPOLOGY_H
#define TETRAHEDRONTOPOLOGY_H

namespace grd {

//
// Topological description of tetrahedron
// __________________________________________________________________
//
// Description:
/*
       3
      /|\
     / | \
    /  |  \
   /   |   \
 0 ----|---- 2
   \   |   /
    \  |  /
     \ | /
      \|/
       1

       3         Faces are numbered with the scheme of face opposite
      /|\
     / | \
    / 1|  \
   / 2 | 0 \
 0 ----|---- 2
   \   |3  /
    \  |  /
     \ | /
      \|/
       1



       3         Edges are numbered according to the followin scheme:
      /|\
   5 / | \3           Edge    Orientation     Edge Number
    /  |  \           0-1         +              2
   /  0|   \          1-2         +              1
 0 ----|---- 2        0-2         +              0
   \   |4  /          2-3         +              3
  2 \  |  /1          1-3         +              4
     \ | /            0-3         +              5
      \|/
       1
*/
// __________________________________________________________________


// __________________________________________________________________
//
// Tetrahedron refinement
// __________________________________________________________________

// refinement into 4 tetrahedra and one octahedron
//
// Vertices of children: 11 Vertices
//      from parent	from parent edges	from parent center
// 	0	0	0	4		10
// 	1	1	1	5
// 	2	2	2	6
// 	3	3	3	7
//			4	8
//			5	9
// Edges of children: from parent edges 6x2, from faces 4x3 : 24
// 	from parent edg		form parent faces
//	0	0,1		0	12,13,14
//	1	2,3		1	15,16,17
//	2	4,5		2	18,19,20
//	3	6,7		3	21,22,23
//	4	8,9		
//	5	10,11		
//



//_____________________________________________________________________
//
// Topology

// vertex indices at each face. Indices are in that order to
// get the outer normal
const int tetVertexOfFace[4][3] = {{1,2,3}, {0,3,2}, {0,1,3}, {0,2,1}};

// vertex indices at each edge
// don't change order; is for refinement
const int tetVertexOfEdge[6][2] = {{0,2},{2,1},{1,0},{2,3},{1,3},{0,3}};

// edge indices at each face, outer normal order
const int tetEdgeOfFace[4][3] = {{1,3,4},{0,5,3},{2,4,5},{0,1,2}};


//_____________________________________________________________________
//
// Refinement
//
// 1. refinement: 5 elements
//		  4 tetras
//		  1 octa
// 2. new vertices: 11
// 3. new edges:    24
// 4. neighbors over parent element faces: 16
// 5. interior neighbors: 5 (tetras + octa)
// 6. total number of neighbors: 	   21
//
//
// Numbering of children
// tetrahedron child is numbered according to the parent vertex it contains
// child[0] contains parent vertex v[0]
// child[1] contains parent vertex v[1]
// child[2] contains parent vertex v[0]
// child[3] contains parent vertex v[0]
//
// octahedron child is the last one
// child[4] = octahedron

// new edges from parent faces
const int tetEdgeFromParentFace[4][3] = {{12,13,14},{15,16,17},
                                          {18,19,20},{21,22,23}};

// child vertex of child edge, coupled with
// tetEdgeFromParentFace
const int tetVertexOfEdgeOfParentFace[4][3][2] = {{{8,5},{5,7},{7,8}},
                                                    {{4,9},{9,7},{7,4}},
                                                    {{9,6},{6,8},{8,9}},
                                                    {{6,4},{4,5},{5,6}}};

// child at parent face by vertex
const int tetChildOfParentFace[4][4] = {{1,2,3,4},{0,3,2,4},
					 {0,1,3,4},{0,2,1,4}};
// child face at parent face
const int tetChildFaceOfParentFace[4][4] = {{0,0,0,1},{1,1,1,6},
					     {2,2,2,3},{3,3,3,4}};


// child edge of child at face
const int tetChildEdgeOfParentFace[4][3] = {{3,4,1},{3,0,5},
					     {4,5,2},{1,2,0}};


// the new vertices, edges and neighbors are incrementally numbred
// assemble children from this elements

// Vertices of tetra children
const int tetVertexOfTetraChild[4][4] = {{0,6,4,9},{6,1,5,8},
                                          {4,5,2,7},{9,8,7,3}};

// Vertices of octa includes barycenter
const int tetVertexOfOctaChild[7] = {8,6,5,7,9,4,10};


// Edges of tetra children
const int tetEdgeOfTetraChild[4][6] = {{0,21,5,15,18,10},
                                        {23,3,4,12,8,19},
                                        {1,2,22,6,13,17},
                                        {16,14,20,7,9,11}};

// Edges of octa child
const int tetEdgeOfOctaChild[12] = {19,12,14,20,23,13,16,18,21,22,17,15};


// _____________________________________________________________________
//
// information necessary to get vertices, edges and children
// neighbors from refined neighbor element. Take vertex of face defintion
// into account for orientation and consistency


// neighbor over parent face number according to outer normal at
// parent face and to parent face order
const int tetNeighborOfParentFace[4][4] = {{0,1,2,3},{4,5,6,7},
                                            {8,9,10,11},{12,13,14,15}};


// tetra child external connection or connection over parent face
const int tetNeighOfFaceOfTetraChild[4][4] = {{20,4,8,12},{0,20,9,14},
                                               {1,6,20,13},{2,5,10,20}};


// Octa chil neighbors
const int tetNeighOfFaceOfOctaChild[8] = {17,3,19,11,15,18,7,16};


// child neighbor of external neighbor
const int tetChildNeighOfNeighbor[16] = {1,2,3,4,0,3,2,4,
                                          0,1,3,4,0,2,1,4};

// _____________________________________________________________________
//
//	Further topological information
// 	not actually used for the refinement algorithm
//

// number of vertices for each face
// const int tetNoOfVerticesOfFace[4] = {3,3,3,3};

// number of edges for each face
// const int tetNoOfEdgesOfFace[4] = {3,3,3,3};

//  vertex index of first vertex at edge
// const int tetFirstVertexOfEdge[6] = {0,2,1,2,1,0};

//  vertices of a given edge
//const int tetVerticesOfEdge[6][2] = {{0,1},{1,2},{2,0},
//                                      {0,3},{1,3},{2,3}};

// // edge index between two vertices
// const int tetEdgeWithVertices[4][4] = {{-1,2,0,5},{2,-1,1,4},{0,1,-1,3},
// 				     {5,4,3,-1}};
// // face indices at each edge
// const int tetFaceWithEdge[6][2] = {{1,3},{0,3},{2,3},{0,1},{0,2},{1,2}};

// // new vertices from parent vertices
// const int tetVertexFromParentVertex[4] = {0,1,2,3};

// // new vertices from parent edge midpoints
// const int tetVertexFromParentEdge[6] = {4,5,6,7,8,9};

// // new vertex from parent barycenter
// const int tetVertexFromParentBarycenter[1] = {10};

// // new edges from parent edges
// const int tetEdgeFromParentEdge[6][2] = {{0,1},{2,3},{4,5},
// 					  {6,7},{8,9},{10,11}};
// // Children at Face	
// const int tetChildOfFace[4][4] = {{1,2,3,4},{0,3,2,4},{0,1,3,4},
// 				   {0,2,1,4}};

// // edge of child opposite to vertex of child at face
// // remark: use with vertex index
// const int tetEdgeOfFace[4][3] = {{3,4,1},{3,0,5},{4,5,2},{1,2,0}};


// // same parent-child vertex of child at face
// const int tetVertexOfChildOfFace[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

// _____________________________________________________________________


// _____________________________________________________________________

// 	TETRAHEDRON CLOSURE


// map to build refinement pattern into refinement case
// case correspondence
//	case	pattern			pattern id
//	Empty	null pattern		0
//	A 	one edge		1
//	B	two edges on face	2
//	C	two oposite edges	3
//	D	three edges on face	4
//	E	retriangulation		5
//	Full	full pattern		6
//
// Pattern address is computed bit-wise. If an edge is refined, set
// corresponding bit
//
// Remark: vector generated with the program: GenerateRefPappernMap

// const int refPatternMap[64] = {0,1,1,2,1,2,2,4,
// 				1,2,2,5,3,5,5,5,
// 				1,3,2,5,2,5,5,5,
// 				2,5,4,5,5,5,5,5,
// 				1,2,3,5,2,5,5,5,
// 				2,4,5,5,5,5,5,5,
// 				2,5,5,5,4,5,5,5,
// 				5,5,5,5,5,5,5,6};





// ____________________________________________________________
//
//	
//	         Tables for slice tetras
//
// ____________________________________________________________
//

// change edge numbering
const int edgesTetra[6][2] = {{0,1},{1,2},{2,0},
                               {0,3},{1,3},{2,3}};

// Description
// table of cases by slicing a tetra
typedef int EDGE_LIST_TETRA;
typedef struct {
  EDGE_LIST_TETRA edges[7];
} TRIANGLE_CASES_TETRA;

const TRIANGLE_CASES_TETRA triCasesTetra[] = {
  {{-1, -1, -1, -1, -1, -1, -1}},
  {{ 0, 3, 2, -1, -1, -1, -1}},
  {{ 0, 1, 4, -1, -1, -1, -1}},
  {{ 3, 2, 4, 4, 2, 1, -1}},
  {{ 1, 2, 5, -1, -1, -1, -1}},
  {{ 3, 5, 1, 3, 1, 0, -1}},
  {{ 0, 2, 5, 0, 5, 4, -1}},
  {{ 3, 5, 4, -1, -1, -1, -1}},
  {{ 3, 4, 5, -1, -1, -1, -1}},
  {{ 0, 4, 5, 0, 5, 2, -1}},
  {{ 0, 5, 3, 0, 1, 5, -1}},
  {{ 5, 2, 1, -1, -1, -1, -1}},
  {{ 3, 4, 1, 3, 1, 2, -1}},
  {{ 0, 4, 1, -1, -1, -1, -1}},
  {{ 0, 2, 3, -1, -1, -1, -1}},
  {{-1, -1, -1, -1, -1, -1, -1}}
};

const int CASE_MASK_TETRA[4] = {1,2,4,8};



// ____________________________________________________________
//
//	
//	         Tables for numerical integration
//
// ____________________________________________________________
//


//const double LumpedTetraOrder2[10][3] =
//{
//  {0.0,0.0,0.0},
//  {1.0,0.0,0.0},
//  {0.0,1.0,0.0},
//  {0.0,0.0,1.0},
//  {0.5,0.0,0.0},
//  {0.5,0.5,0.0},
//  {0.0,0.5,0.0},
//  {0.0,0.5,0.5},
//  {0.5,0.0,0.5},
//  {0.0,0.0,0.5}
//};
//
//
//
//const double GaussTetraOrder3[5][3] =
//{
//  {0.250000000000000,0.250000000000000,0.250000000000000},
//  {0.166666666666667,0.166666666666667,0.166666666666667},
//  {0.500000000000000,0.166666666666667,0.166666666666667},
//  {0.166666666666667,0.500000000000000,0.166666666666667},
//  {0.166666666666667,0.166666666666667,0.500000000000000}
//};
//
//
//const double GaussTetraOrder5[15][3] =
//{
//  {0.25000000000000000,0.25000000000000000,0.25000000000000000},
//  {0.09197107805272303,0.09197107805272303,0.09197107805272303},
//  {0.72408676584183096,0.09197107805272303,0.09197107805272303},
//  {0.09197107805272303,0.72408676584183096,0.09197107805272303},
//  {0.09197107805272303,0.09197107805272303,0.72408676584183096},
//  {0.44364916731037080,0.05635083268962915,0.05635083268962915},
//  {0.05635083268962915,0.44364916731037080,0.05635083268962915},
//  {0.05635083268962915,0.05635083268962915,0.44364916731037080},
//  {0.05635083268962915,0.44364916731037080,0.44364916731037080},
//  {0.44364916731037080,0.05635083268962915,0.44364916731037080},
//  {0.44364916731037080,0.44364916731037080,0.05635083268962915},
//  {0.31979362782962989,0.31979362782962989,0.31979362782962989},
//  {0.04061911651111023,0.31979362782962989,0.31979362782962989},
//  {0.31979362782962989,0.04061911651111023,0.31979362782962989},
//  {0.31979362782962989,0.31979362782962989,0.04061911651111023}
//};


} // namespace grd

#endif // TETRAHEDRONTOPOLOGY_H
