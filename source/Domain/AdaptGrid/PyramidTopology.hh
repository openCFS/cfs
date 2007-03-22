// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// - C++ -
/***************************************************************************
    File        : PyramidTopology.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Thu Dec 12 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef PYRAMIDTOPOLOGY_H
#define PYRAMIDTOPOLOGY_H

namespace grd {

// __________________________________________________________________
//
// Topological description of a pyramid
// __________________________________________________________________


// vertex indices at each face. Indices are in that order to
// get the outer normal
const int pyrNoOfVertexOfFace[5] = {3,3,3,3,4};
const int pyrVertexOfFace[5][4]  = {{0,1,2,-1},{0,2,3,-1},{0,3,4,-1},{0,4,1,-1},{1,4,3,2}};

// vertex indices at each edge
const int pyrVertexOfEdge[8][2] = {{0,1},{0,2},{0,3},{0,4},
                                   {1,2},{2,3},{3,4},{4,1}};

// edge indices at each face
const int pyrNoOfEdgeOfFace[5] = {3,3,3,3,4};
const int pyrEdgeOfFace[5][4]  = {{0,4,1,-1},{1,5,2,-1},{2,6,3,-1},{3,7,0,-1},{7,6,5,4}};



//_____________________________________________________________________
//
// Refinement
//
// 1. refinement: 10 elements
//      6 pyras
//      4 tetras
// 2. new vertices: 14
//  from parent vert:  5
//  from parent edgs:  8
//  from parent face:  1
// 3. new edges:    36
//  from parent edges: 16
//  from parent tri faces: 12
//  from parent qua faces:  4
//  from parent interior:   4 (interior edges)
// 4. neighbors over parent element faces: 20 (faces x 4 = 5 x 4)
// 5. interior neighbors: 10  (tetras + pyras)
// 6. total number of neighbors: 30
//
// Numbering of children
// pyramid child is numbered according to the parent vertex it contains
// the interior pyramid becomes the last number
// child[0] contains parent vertex v[0]
// child[1] contains parent vertex v[1]
// child[2] contains parent vertex v[2]
// child[3] contains parent vertex v[3]
// child[4] contains parent vertex v[4]
// child[5] interior pyramid
//
// tetrahedron child is numberd according to the parent face it contains
// child[6] contains parent face n[0]
// child[7] contains parent face n[1]
// child[8] contains parent face n[2]
// child[9] contains parent face n[3]


// new edges from parent faces
const int pyrNoOfEdgeFromParentFace[5] = {3,3,3,3,4};
const int pyrEdgeFromParentFace[5][4]  = {{16,17,18,-1},{19,20,21,-1},{22,23,24,-1},
                                          {25,26,27,-1},{28,29,30,31}};

// child vertex of child edge at parent face
// coupled with pyrEdgeFromParentFace
const int pyrVertexOfEdgeOfParentFace[5][4][2] = {{{6,5},{5,9},{9,6},{-1,-1}},
                                                  {{7,6},{6,10},{10,7},{-1,-1}},
                                                  {{8,7},{7,11},{11,8},{-1,-1}},
                                                  {{5,8},{8,12},{12,5},{-1,-1}},
                                                  {{13,9},{13,10},{13,11},{13,12}}};


// child at parent face by vertex
const int pyrNoOfChildOfParentFace[5] = {4,4,4,4,4};
const int pyrChildOfParentFace[5][4]  = {{0,1,2,6},{0,2,3,7},
                                         {0,3,4,8},{0,4,1,9},
                                         {1,4,3,2}};

// child face at parent face
const int pyrChildFaceOfParentFace[5][4] = {{0,0,0,3},{1,1,1,3},
                                            {2,2,2,3},{3,3,3,3},
                                            {4,4,4,4}};

// edge of child at parent face
const int pyrChildEdgeOfParentFace[5][4] = {{4,1,0,-1},{5,2,1,-1},
                                            {6,3,2,-1},{7,0,3,-1},
                                            {5,6,7,4}};


// interior edge from parent edge
const int pyrInteriorEdgeFromParentEdge[4] = {32,33,34,35};

// vertex of interior edge from face
const int pyrVertexOfIntEdgeFromParentEdge[4][2] = {{5,13},{6,13},{7,13},{8,13}};


// vertices of pyramid children
const int pyrVertexOfPyramidChild[6][5] = {{0,5,6,7,8},{5,1,9,13,12},
                                           {6,9,2,10,13},{7,13,10,3,11},
                                           {8,12,13,11,4},{13,7,6,5,8}};

// edge of pyramid children
const int pyrEdgeOfPyramidChild[6][8] = {{0,2,4,6,16,19,22,25},{1,17,32,27,8,28,29,15},
										{18,3,20,33,9,10,31,28},{34,21,5,23,11,12,30},
										{26,35,24,7,29,30,13,14},{34,33,32,35,19,16,25,22}};

// vertex of tetra children
const int pyrVertexOfTetraChild[4][4] = {{5,6,9,13},{6,7,10,13},
                                         {7,8,11,13},{8,5,12,13}};

// edge of tetra child
const int pyrEdgeOfTetraChild[4][6] = {{17,18,16,28,33,32},
                                       {20,21,19,31,34,33},
                                       {23,14,22,30,35,34},
                                       {26,27,25,29,32,35}};


// _____________________________________________________________________
//
// information necessary to get vertices, edges and children
// neighbors from refined neighbor element. Take vertex of face
// defintion into account for orientation and consistency


// XXXXXXXXXX here has to be checked!!!!***********************

// neighbor over parent face
const int pyrNeighborOfParentFace[5][4] = {{0,1,2,3},{4,5,6,7},
                                           {8,9,10,11},{12,13,14,15},
                                           {16,17,18,19}};


// child pyramid connenction
const int pyrPyramidChildNeighOfFace[6][5] = {{0,4,8,12,25},
                                              {1,26,29,14,16},
                                              {2,5,27,26,17},
                                              {27,6,8,28,18},
                                              {29,28,10,13,19},
                                              {28,27,26,29,20}};


// child tetrahedron connection
const int pyrTetraChildNeighborOfFace[4][4] = {{22,21,25,3},{23,22,25,7},
                                               {24,23,25,11},{21,24,25,15}};

// ************ to be done from here ***********
// child neighbor of external neighbor
// connect external neighbors with new children
const int pyrChildNeighOfNeighbor[20] = {0,1,2,6,
					                               0,2,3,7,
					                               0,3,4,8,
					                               0,4,1,9,
                                         1,2,3,4};



//_____________________________________________________________________
//
//	Further topologycal information
//	not used by refinement algorithm

// ... no further topologycal information stored


//_____________________________________________________________________
// Irregular refinement

// The symmetries given here are linear isometries
// The isometries 0 to 3 are the identity and rotation, i.e. with det = 1
// The isometries 5 to 7 are obtained from the previous one by reflexion
// i.e. they have det = -1
const int pyrIsometryVertex[8][5] = {{0,1,2,3,4},{0,2,3,4,1},{0,3,4,1,2},{0,4,1,2,3},
                                     {0,4,3,2,1},{0,2,1,4,3},{0,1,4,3,2},{0,3,2,1,4}};
const int pyrIsometryEdge[8][8]   = {{0,1,2,3,4,5,6,7},{1,2,3,0,5,6,7,4},{2,3,0,1,6,7,4,5},{3,0,1,2,7,4,5,6},
                                     {3,2,1,0,6,5,4,7},{1,0,3,2,4,7,6,5},{0,3,2,1,7,6,5,4},{2,1,0,3,5,4,7,6}};
// Equivalence cases and rotations
// First number is rotation
// Second number is equivalence case
const size_t pyrNoOfSymmetries = 8;
const size_t pyrIsometryCases[256][2] = {
{0,0},{0,1},{1,1},{0,3},{2,1},{0,5},{1,3},{0,7},
{3,1},{3,3},{1,5},{3,7},{2,3},{2,7},{1,7},{0,15},
{0,16},{0,17},{5,17},{0,19},{0,20},{0,21},{0,22},{0,23},
{5,20},{5,22},{5,21},{5,23},{0,28},{0,29},{5,29},{0,31},
{1,16},{7,20},{1,17},{7,22},{7,17},{7,21},{1,19},{7,23},
{1,20},{1,28},{1,21},{1,29},{1,22},{7,29},{1,23},{1,31},
{0,48},{0,49},{0,50},{0,51},{7,49},{0,53},{7,51},{0,55},
{0,56},{0,57},{0,58},{0,59},{7,57},{0,61},{7,59},{0,63},
{2,16},{2,20},{4,20},{2,28},{2,17},{2,21},{4,22},{2,29},
{4,17},{2,22},{4,21},{4,29},{2,19},{2,23},{4,23},{2,31},
{0,80},{0,81},{5,81},{0,83},{2,81},{0,85},{0,86},{0,87},
{4,81},{2,86},{4,85},{5,87},{2,83},{2,87},{4,87},{0,95},
{1,48},{1,56},{1,49},{1,57},{1,50},{1,58},{1,51},{1,59},
{4,49},{4,57},{1,53},{1,61},{4,51},{4,59},{1,55},{1,63},
{0,112},{0,113},{0,114},{0,115},{4,114},{0,117},{0,118},{0,119},
{4,113},{0,121},{4,117},{0,123},{4,115},{4,123},{4,119},{0,127},
{3,16},{6,17},{3,20},{3,22},{6,20},{6,21},{3,28},{6,29},
{3,17},{3,19},{3,21},{3,23},{6,22},{6,23},{3,29},{3,31},
{3,48},{3,50},{5,49},{5,51},{3,56},{3,58},{5,57},{5,59},
{3,49},{3,51},{3,53},{3,55},{3,57},{3,59},{3,61},{3,63},
{1,80},{6,81},{1,81},{3,86},{7,81},{6,85},{1,83},{7,87},
{3,81},{3,83},{1,85},{3,87},{1,86},{6,87},{1,87},{1,95},
{3,112},{3,114},{7,114},{3,118},{7,113},{7,117},{7,115},{7,119},
{3,113},{3,115},{3,117},{3,119},{3,121},{3,123},{7,123},{3,127},
{2,48},{6,49},{2,56},{6,57},{2,49},{2,53},{2,57},{2,61},
{2,50},{6,51},{2,58},{6,59},{2,51},{2,55},{2,59},{2,63},
{2,112},{5,114},{5,113},{5,115},{2,113},{2,117},{2,121},{5,123},
{2,114},{2,118},{5,117},{5,119},{2,115},{2,119},{2,123},{2,127},
{1,112},{6,113},{1,113},{1,121},{1,114},{6,117},{1,115},{1,123},
{6,114},{6,115},{1,117},{6,123},{1,118},{6,119},{1,119},{1,127},
{0,240},{0,241},{1,241},{0,243},{2,241},{0,245},{1,243},{0,247},
{3,241},{3,243},{1,245},{3,247},{2,243},{2,247},{1,247},{0,255},
};


} // namespace grd
#endif // PYRAMIDTOPOLOGY_H
