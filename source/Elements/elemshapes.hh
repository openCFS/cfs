// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_ELEMSHAPES_HH
#define FILE_CFS_ELEMSHAPES_HH

#include <map>

#include "General/environment.hh"
#include "Domain/elem.hh"
#include "MatVec/vector.hh"

namespace CoupledField {

//! Description of geometry of reference elements
//! ToDO: Move array numElemNodes, elemDims_, elemQuadratic_ and 
//!       definition of feType from elem.hh to here
struct ElemShape {

  //! dimension of element
  UInt dim;

  //! number of vertices (corners)
  UInt numVertices;

  //! number of nodes
  UInt numNodes;

  //! number of edges
  UInt numEdges;

  //! number of faces
  UInt numFaces;

  //! coordinate of element midpoint
  Vector<Double> midPointCoord;

  //! coordinates of nodes
  StdVector<StdVector<Double> > nodeCoords; 

  //! indices (node numbers) of edges for each edge 
  StdVector<StdVector<UInt> > edges;

  //! indices (node numbers) of faces
  StdVector<StdVector<UInt> > faces;
};

//! Collection of geometric shapes of all element shapes
struct ElemShapes {

  //! Collection of reference element shape
  static StdVector<ElemShape> shapes;

  //! Fill information
  static void Initialize();

private:

  //! Helper method for setting element shape information (1D)
  static void SetElemInfo( ElemShape& shape, Double midPoint[], Double nodeCoords[],
                              UInt edgeIndices[], UInt numFaceNodes[] = NULL, UInt faceIndices[] = NULL );
};


} // namespace Coupledfield
#endif
