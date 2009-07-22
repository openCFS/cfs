// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "elemshapes.hh"

namespace CoupledField {


  void ElemShapes::Initialize() {
  
  // ************************************************************************
  //  POINT
  // ************************************************************************
  
  // ************************************************************************
  // LINE2
  // ************************************************************************
  {
    ElemShape s;
    s.dim = 1;
    s.numVertices = 2;
    s.numNodes = 2;
    s.numEdges = 1;
    s.numFaces = 0;
    Double midPoint[1] = {0.0};
    Double nodeCoords[] = 
    { 
     -1.0, // #1
      1.0  // #2 
    };
    UInt edgeIndices[] = 
    { 
     1, 2 // #1
    };
    SetElemInfo( s, midPoint, nodeCoords, edgeIndices ); 
    shapes[Elem::LINE2] = s;
  }
  
  // ************************************************************************
  // LINE3
  // ************************************************************************

  // ************************************************************************
  // TRIA3
  // ************************************************************************
  
  // ************************************************************************
  // TRIA6
  // ************************************************************************
  
  // ************************************************************************
  // QUAD4
  // ************************************************************************
  {
    ElemShape s;
    s.dim = 2;
    s.numVertices = 4;
    s.numNodes = 4;
    s.numEdges = 4;
    s.numFaces = 1;
    Double midPoint[2] = {0.0, 0.0};
    Double nodeCoords[] = 
    { 
     -1.0, -1.0, // #1
      1.0, -1.0, // #2
      1.0,  1.0, // #3
     -1.0,  1.0  // #4 
    };
    UInt edgeIndices[] = 
    { 
     1, 2, // #1
     2, 3, // #2
     4, 3, // #3
     1, 4  // #4 
    };
    UInt numFaceNodes[] = 
    {
     4 // #1
    };
    UInt faceIndices[] = 
    {
     1, 2, 3, 4 // #1
    };

    SetElemInfo( s, midPoint, nodeCoords, edgeIndices, 
                 numFaceNodes, faceIndices ); 
    shapes[Elem::QUAD4] = s;
  }

  // ************************************************************************
  // QUAD8
  // ************************************************************************
  
  // ************************************************************************
  // QUAD9
  // ************************************************************************
  
  // ************************************************************************
  // TET4 
  // ************************************************************************
  
  // ************************************************************************
  // TET10
  // ************************************************************************
  
  // ************************************************************************
  // HEXA8
  // ************************************************************************
  {
    ElemShape s;
    s.dim = 3;
    s.numVertices = 8;
    s.numNodes = 8;
    s.numEdges = 12;
    s.numFaces = 6;
    Double midPoint[3] = {0.0, 0.0, 0.0};
    Double nodeCoords[] = 
    { -1.0, -1.0, -1.0, // #1
       1.0, -1.0, -1.0, // #2
       1.0,  1.0, -1.0, // #3
      -1.0,  1.0, -1.0, // #4
      -1.0, -1.0,  1.0, // #5
       1.0, -1.0,  1.0, // #6
       1.0,  1.0,  1.0, // #7
      -1.0,  1.0,  1.0  // #8
    };
    UInt edgeIndices[] = 
    { 1, 2, // #1
      2, 3, // #2
      4, 3, // #3
      1, 4, // #4
      1, 5, // #5
      2, 6, // #6
      3, 7, // #7
      4, 8, // #8
      5, 6, // #9
      6, 7, // #10
      8, 7, // #11
      5, 8  // #12
    };
    UInt numFaceNodes[] = 
    {
     4, // #1
     4, // #2
     4, // #3
     4, // #4
     4, // #5
     4  // #6
    };
    UInt faceIndices[] = 
    {
     1, 2, 3, 4, // #1
     1, 2, 6, 5, // #2
     2, 3, 7, 6, // #3
     4, 3, 7, 8, // #4
     1, 4, 8, 5, // #5
     5, 6, 7, 8  // #6
    };
    SetElemInfo( s, midPoint, nodeCoords, edgeIndices,
                 numFaceNodes, faceIndices ); 
    shapes[Elem::HEXA8] = s;
  }
  // ************************************************************************
  // HEXA20
  // ************************************************************************
  
  // ************************************************************************
  // HEXA27
  // ************************************************************************
  
  // ************************************************************************
  // PYRA5
  // ************************************************************************
  
  // ************************************************************************
  // PYRA13
  // ************************************************************************
  
  // ************************************************************************
  // WEDGE6
  // ************************************************************************
  
  // ************************************************************************
  // WEDGE15
  // ************************************************************************

  }
  
  void ElemShapes::SetElemInfo( ElemShape& shape, Double midPoint[], 
                                Double nodeCoords[], UInt edgeIndices[], 
                                UInt numFaceNodes[], UInt faceIndices[] ) {

    // midPoint
    UInt dim = shape.dim;
    shape.midPointCoord.Resize( dim );
    for( UInt iDim = 0; iDim < dim; iDim++ ) {
      shape.midPointCoord[iDim] = midPoint[iDim];
    }

    // nodal coordinates
    shape.nodeCoords.Resize( shape.numNodes );
    for( UInt iNode = 0; iNode < shape.numNodes; iNode++ ) {
      shape.nodeCoords[iNode].Resize( dim);
      for( UInt iDim = 0; iDim < dim; iDim++ ) {
        (shape.nodeCoords[iNode])[iDim] = nodeCoords[iNode*dim+iDim];
      }
    }

    // edgeIndices
    shape.edges.Resize( shape.numEdges );
    for( UInt iEdge = 0; iEdge < shape.numEdges; iEdge++ ) {
      shape.edges[iEdge].Resize( 2 );
      shape.edges[iEdge][0] = edgeIndices[2*iEdge];
      shape.edges[iEdge][1] = edgeIndices[2*iEdge+1];
    }

    // faceIndices
    shape.faces.Resize( shape.numFaces );
    UInt pos = 0;
    for( UInt iFace = 0; iFace < shape.numFaces; iFace++ ) {
      shape.faces[iFace].Resize( numFaceNodes[iFace] );
      for( UInt iNode = 0; iNode < numFaceNodes[iFace]; iNode++ ) {
        shape.faces[iFace][iNode] = faceIndices[pos++];
      }
    }
  }


}
