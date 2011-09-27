// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "elem.hh"
#include "Elements/basefe.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "General/Enum.hh"

namespace CoupledField {

// static declaration
std::map<Elem::FEType,ElemShape> Elem::shapes;

  Elem::Elem() : 
    elemNum(0),
    type(ET_UNDEF),
    regionId( NO_REGION_ID )
  {}

//

  std::string Elem::ToString() const {
    std::ostringstream os;
    os << "elemNum=" << elemNum << " region=" << regionId;
    return os.str();
  }
  
 Elem & Elem::operator=(const Elem& t) {
   if (this!=&t) {
     elemNum = t.elemNum;
     type = t.type;
     regionId=t.regionId;
     connect = t.connect;
     edges = t.edges;
     faces = t.faces;
     faceFlags = t.faceFlags;
   }
   return *this;
 }
  void Elem::CorrectConnectivity() 
  {
    UInt dummy;

    // This funtion rearanges the connectivity of the element so
    // that the orientation is in such a way that the Jacobion determinant
    // will always be positive. 
    switch(type) 
    {
    case ET_TRIA6:
      dummy = connect[4];
      connect[4] = connect[5];
      connect[5] = dummy;
    case ET_TRIA3:
      dummy = connect[0];
      connect[0] = connect[1];
      connect[1] = dummy;
      break;
    case ET_QUAD8:
      dummy = connect[4];
      connect[4] = connect[7];
      connect[7] = dummy;
      dummy = connect[5];
      connect[5] = connect[6];
      connect[6] = dummy;
    case ET_QUAD4:
      dummy = connect[1];
      connect[1] = connect[3];
      connect[3] = dummy;
      break;
    default:
      EXCEPTION("Connectivity for " << feType.ToString(type) << " element "
                << elemNum << " in region "
                << domain->GetGrid()->GetRegion().ToString(regionId)
                << " is not properly oriented!" );
      break;
    }
  }
  
  void SetElemInfo( ElemShape& shape, Double midPoint[], 
                      Double nodeCoords[], UInt edgeVertices[],
                      UInt numEdgeNodes[], UInt edgeNodes[],
                      UInt numEdgesPerDir[],  UInt locDirEdges[],
                      UInt numFaceVertices[] = NULL, UInt faceVertices[] = NULL,
                      UInt numFaceNodes[] = NULL, UInt faceNodes[] = NULL,
                      UInt numFacesPerDir[] = NULL, UInt locDirFaces[][2] = NULL) {

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

      // edgeVertices
      shape.edgeVertices.Resize( shape.numEdges );
      for( UInt iEdge = 0; iEdge < shape.numEdges; iEdge++ ) {
        shape.edgeVertices[iEdge].Resize( 2 );
        shape.edgeVertices[iEdge][0] = edgeVertices[2*iEdge];
        shape.edgeVertices[iEdge][1] = edgeVertices[2*iEdge+1];
      }
      
      // edgeNodes
      UInt pos = 0;
      shape.edgeNodes.Resize( shape.numEdges );
      for( UInt iEdge = 0; iEdge < shape.numEdges; iEdge++ ) {
        shape.edgeNodes[iEdge].Resize( numEdgeNodes[iEdge] );
        for( UInt iNode = 0; iNode < numEdgeNodes[iEdge]; iNode++ ) {
          shape.edgeNodes[iEdge][iNode] = edgeNodes[pos++];
        }
      }
      
      // locDirEdges
      pos = 0;
      shape.locDirEdges.Resize( shape.dim );
      for( UInt iDim = 0; iDim < shape.dim; iDim++ ) {
        shape.locDirEdges[iDim].Resize( numEdgesPerDir[iDim]);
        for( UInt iEdge = 0; iEdge < numEdgesPerDir[iDim]; iEdge++ )
          shape.locDirEdges[iDim][iEdge] = locDirEdges[pos++];
      }

      // faceVertices
      shape.faceVertices.Resize( shape.numFaces );
      pos = 0;
      for( UInt iFace = 0; iFace < shape.numFaces; iFace++ ) {
        shape.faceVertices[iFace].Resize( numFaceVertices[iFace] );
        for( UInt iNode = 0; iNode < numFaceVertices[iFace]; iNode++ ) {
          shape.faceVertices[iFace][iNode] = faceVertices[pos++];
        }
      }
      
      // faceNodes
      shape.faceNodes.Resize( shape.numFaces );
      pos = 0;
      for( UInt iFace = 0; iFace < shape.numFaces; iFace++ ) {
        shape.faceNodes[iFace].Resize( numFaceNodes[iFace] );
        for( UInt iNode = 0; iNode < numFaceNodes[iFace]; iNode++ ) {
          shape.faceNodes[iFace][iNode] = faceNodes[pos++];
        }
      }
      
      // locDirFaces
      if( shape.numFaces > 0 ) {
        pos = 0;
        shape.locDirFaces.Resize( shape.dim );
        for( UInt iDim = 0; iDim < shape.dim; iDim++ ) {
          shape.locDirFaces[iDim].Resize( numFacesPerDir[iDim]);
          for( UInt iFace = 0; iFace < numFacesPerDir[iDim]; iFace++ ) {
            shape.locDirFaces[iDim][iFace] = 
                std::pair<UInt,UInt>(locDirFaces[pos][0], 
                                     locDirFaces[pos][1]);
            pos++;
          }
        }
      }
    }


  ElemShape::ElemShape():
        dim(0),
        order(0),
        numVertices(0),
        numNodes(0),
        numEdges(0),
        numFaces(0) {
    // nothing to do here
  }

    void ElemShape::Initialize() {
      
    // ************************************************************************
    //  POINT
    // ************************************************************************
    
    // ************************************************************************
    // LINE2
    // ************************************************************************
    {
      ElemShape s;
      s.dim = 1;
      s.order = 1;
      s.numVertices = 2;
      s.numNodes = 2;
      s.numNodes = 2;
      s.numEdges = 1;
      s.numFaces = 0;
      
      Double midPoint[1] = {0.0};
      Double nodeCoords[] = 
      { 
       -1.0, // #1
        1.0  // #2 
      };
      UInt edgeVertices[] = 
      { 
       1, 2 // #1
      };
      UInt numEdgeNodes[] = 
      { 
       2 // #1
      };
      UInt numEdgesPerDir[] =
      {
       1 // #edges in xi-dir
      };
      UInt locDirEdges[] =
      {
       1 // xi
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeVertices, numEdgesPerDir, locDirEdges, 
                   NULL, NULL, NULL, NULL );
      Elem::shapes[Elem::ET_LINE2] = s;
    }
    
    // ************************************************************************
    // LINE3
    // ************************************************************************

    // ************************************************************************
    // TRIA3
    // ************************************************************************
    //            
    //            
    // 3 +        
    //   |\         eta
    //   | \       ^        REFERENCE VOLUME ELEMENT
    //   |  \      | 
    // 1 +---+ 2   +--> xi
    
    ElemShape s;
    s.dim = 2;
    s.order = 1;
    s.numVertices = 3;
    s.numNodes = 3;
    s.numEdges = 3;
    s.numFaces = 1;
    Double midPoint[2] = {1.0/0.3, 1.0/0.3};
    Double nodeCoords[] = 
    { 
      0.0,  0.0, // #1
      1.0,  0.0, // #2
      0.0,  1.0, // #3
    };
    UInt edgeVertices[] = 
    { 
     1, 2, // #1
     2, 3, // #2
     1, 3, // #3
    };
    UInt numEdgeNodes[] = 
    {
     2, // #1
     2, // #2
     2, // #3
    };
    UInt * edgeNodes  = edgeVertices; // nodes = vertices
    UInt numEdgesPerDir[] =
    {
     2, // #edges in xi-dir
     2  // #edges in eta-dir
    };
    UInt locDirEdges[] =
    {
     1,2, // xi
     2,3  // eta
    };
    UInt numFaceVertices[] = 
    {
     3 // #1
    };
    UInt faceVertices[] = 
    {
     1, 2, 3 // #1
    };
    UInt * numFaceNodes = numFaceVertices;
    UInt * faceNodes = faceVertices;
    UInt numFacesPerDir[] = 
    {
     1, // #faces in xi-dir
     1  // #faces in eta-dir
    };
    UInt locDirFaces[][2] = 
    {
     {1,0}, // #faces in xi-dir & component
     {1,1}, // #faces in eta-dir & component
    };

    SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                 edgeNodes, numEdgesPerDir, locDirEdges,
                 numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                 numFacesPerDir, locDirFaces );
    Elem::shapes[Elem::ET_TRIA3] = s;
    
    
    // ************************************************************************
    // TRIA6
    // ************************************************************************
    
    // ************************************************************************
    // QUAD4
    // ************************************************************************
    {
      //      eta
      //       ^
      // 4 +---|---+ 3    
      //   |   |   |      
      //   |   0---|-> xi     REFERENCE VOLUME ELEMENT
      //   |       |
      // 1 +-------+ 2
      
      ElemShape s;
      s.dim = 2;
      s.order = 1;
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
      UInt edgeVertices[] = 
      { 
       1, 2, // #1
       2, 3, // #2
       4, 3, // #3
       1, 4  // #4 
      };
      UInt numEdgeNodes[] = 
      {
       2, // #1
       2, // #2
       2, // #3
       2  // #4
      };
      UInt * edgeNodes  = edgeVertices; // nodes = vertices
      UInt numEdgesPerDir[] =
      {
       2, // #edges in xi-dir
       2  // #edges in eta-dir
      };
      UInt locDirEdges[] =
      {
       1,3, // xi
       2,4  // eta
      };
      UInt numFaceVertices[] = 
      {
       4 // #1
      };
      UInt faceVertices[] = 
      {
       1, 2, 3, 4 // #1
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      UInt numFacesPerDir[] = 
      {
       1, // #faces in xi-dir
       1  // #faces in eta-dir
      };
      UInt locDirFaces[][2] = 
      {
       {1,0}, // #faces in xi-dir & component
       {1,1}, // #faces in eta-dir & component
      };

      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, numEdgesPerDir, locDirEdges,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   numFacesPerDir, locDirFaces );
      Elem::shapes[Elem::ET_QUAD4] = s;
    }

    // ************************************************************************
    // QUAD8
    // ************************************************************************
    {
      //      eta
      //       ^
      //       7
      // 4 +---+---+ 3    
      //   |   |   |      
      // 8 +   0---+6> xi     REFERENCE VOLUME ELEMENT
      //   |       |
      // 1 +---+---+ 2
      //       5
       ElemShape s;
       s.dim = 2;
       s.order = 1;
       s.numVertices = 4;
       s.numNodes = 8;
       s.numEdges = 4;
       s.numFaces = 1;
       Double midPoint[2] = {0.0, 0.0};
       Double nodeCoords[] = 
       { 
        -1.0, -1.0, // #1
         1.0, -1.0, // #2
         1.0,  1.0, // #3
        -1.0,  1.0, // #4
         0.0, -1.0, // #5
         1.0,  0.0, // #6 
         0.0,  1.0, // #7
        -1.0,  0.0  // #8
       };
       UInt edgeVertices[] = 
       { 
        1, 2, // #1
        2, 3, // #2
        4, 3, // #3
        1, 4  // #4 
       };
       UInt numEdgeNodes[] = 
       {
        3, // #1
        3, // #2
        3, // #3
        3  // #4
       };
       UInt edgeNodes[] = 
       { 
        1, 2, 5, // #1
        2, 3, 6, // #2
        4, 3, 7, // #3
        1, 4, 8  // #4 
       };
       UInt numEdgesPerDir[] =
       {
        2, // #edges in xi-dir
        2  // #edges in eta-dir
       };
       UInt locDirEdges[] =
       {
        1,3, // xi
        2,4  // eta
       };
       UInt numFaceVertices[] = 
       {
        4 // #1
       };
       UInt faceVertices[] = 
       {
        1, 2, 3, 4 // #1
       };
       UInt numFaceNodes[] = 
       {
        8 // #1
       };
       UInt faceNodes[] = 
       {
        1, 2, 3, 4, 5, 6, 7, 8 // #1
       };
       UInt numFacesPerDir[] = 
       {
        1, // #faces in xi-dir
        1  // #faces in eta-dir
       };
       UInt locDirFaces[][2] = 
       {
        {1,0}, // #faces in xi-dir & component
        {1,1}, // #faces in eta-dir & component
       };

       SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, numEdgesPerDir, locDirEdges,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    numFacesPerDir, locDirFaces );
       Elem::shapes[Elem::ET_QUAD8] = s;
     }
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
      //                    zeta 
      //                     ^ eta 
      //    8 +-------+ 7    |/
      //     /|      /|      0--> xi 
      //    / |     / |
      // 5 +--+----+6 |   
      //   |  +-- -|- + 3    
      //   | / 4   | /    REFERENCE VOLUME ELEMENT
      //   |/      |/
      // 1 +-------+ 2

      ElemShape s;
      s.dim = 3;
      s.order = 1;
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
      UInt edgeVertices[] = 
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
      UInt numEdgeNodes[] = 
      {
       2, // #1
       2, // #2
       2, // #3
       2, // #4
       2, // #5
       2, // #6
       2, // #7
       2, // #8
       2, // #9
       2, // #10
       2, // #11
       2  // #12
      };
      UInt * edgeNodes = edgeVertices;
      UInt numEdgesPerDir[] =
      {
       4, // #edges in xi-dir
       4, // #edges in eta-dir
       4, // #edges in zeta-dir
      };
      UInt locDirEdges[] =
      {
       1,3,9,11,  // xi
       2,4,10,12  // eta
      };
      UInt numFaceVertices[] = 
      {
       4, // #1
       4, // #2
       4, // #3
       4, // #4
       4, // #5
       4  // #6
      };
      UInt faceVertices[] = 
      {
       1, 2, 3, 4, // #1
       1, 2, 6, 5, // #2
       2, 3, 7, 6, // #3
       4, 3, 7, 8, // #4
       1, 4, 8, 5, // #5
       5, 6, 7, 8  // #6
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      UInt numFacesPerDir[] = 
      {
       4, // #faces in xi-dir
       4, // #faces in eta-dir
       4  // #faces in zeta-dir
      };
      UInt locDirFaces[][2] = 
      {
       {1,0},{2,0},{4,0},{6,0}, // faces in xi-dir & component
       {1,1},{3,0},{5,0},{6,1}, // faces in eta-dir & component
       {2,1},{3,1},{4,1},{5,1}
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, numEdgesPerDir, locDirEdges,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    numFacesPerDir, locDirFaces );
      Elem::shapes[Elem::ET_HEXA8] = s;
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
    
      
    Elem::ShapeType Elem::GetShapeType( Elem::FEType type )  {
      ShapeType ret = ST_UNDEF;
      switch(type) {
        case ET_UNDEF:
          ret = ST_UNDEF;
          break;
        case ET_POINT:
          ret = ST_POINT;
          break;
        case ET_LINE2:
        case ET_LINE3:
          ret = ST_LINE;
          break;
        case ET_TRIA3:
        case ET_TRIA6:
          ret = ST_TRIA;
          break;
        case ET_QUAD4:
        case ET_QUAD8:
        case ET_QUAD9:
          ret = ST_QUAD;
          break;
        case ET_TET4:
        case ET_TET10:
          ret = ST_TET;
          break;
        case ET_HEXA8:
        case ET_HEXA20:
        case ET_HEXA27:
          ret = ST_HEXA;
          break;
        case ET_PYRA5:
        case ET_PYRA13:
          ret = ST_PYRA;
          break;
        case ET_WEDGE6:
        case ET_WEDGE15:
          ret = ST_WEDGE;
          break;
      }
      return ret;
    }
    
    ElemShape& Elem::GetShape( Elem::ShapeType st) {

      Elem::FEType feType = Elem::ET_UNDEF;
      switch( st ) {
        case ST_POINT:
          feType = ET_POINT;
          break;
        case ST_LINE:
          feType = ET_LINE2;
          break;
        case ST_TRIA:
          feType = ET_TRIA3;
          break;
        case ST_QUAD:
          feType = ET_QUAD4;
          break;
        case ST_TET:
          feType = ET_TET4;
          break;
        case ST_HEXA:
          feType = ET_HEXA8;
          break;
        case ST_PYRA:
          feType = ET_PYRA5;
          break;
        case ST_WEDGE:
          feType = ET_WEDGE6;
          break;
        case ST_UNDEF:
          EXCEPTION("Unknown shape type");
      }
      return Elem::shapes[feType]; 
    }
    // ************************************************************************
    // ENUM INITIALIZATION
    // ************************************************************************
    

    // Definition of finite element type mappings
    static EnumTuple elemShapeTuples[] = {
       EnumTuple(Elem::ST_UNDEF,  "UNDEF"), 
       EnumTuple(Elem::ST_POINT,   "POINT"),
       EnumTuple(Elem::ST_LINE,   "LINE"),
       EnumTuple(Elem::ST_TRIA,   "TRIA"),
       EnumTuple(Elem::ST_QUAD,   "QUAD"),
       EnumTuple(Elem::ST_TET,    "TET"),
       EnumTuple(Elem::ST_HEXA,   "HEX"),
       EnumTuple(Elem::ST_PYRA,   "PYRA"),
       EnumTuple(Elem::ST_WEDGE,  "WEDGE")
    };
    Enum<Elem::ShapeType> Elem::shapeType = \
       Enum<Elem::ShapeType>("Finite Element Shape Types",
           sizeof(elemShapeTuples) / sizeof(EnumTuple),
           elemShapeTuples);
    
    // Definition of finite element type mappings
     static EnumTuple feTypeTuples[] = {
       EnumTuple(Elem::ET_UNDEF,   "UNDEF"),
       EnumTuple(Elem::ET_POINT,   "POINT"),
       EnumTuple(Elem::ET_LINE2,   "LINE2"),
       EnumTuple(Elem::ET_LINE3,   "LINE3"),
       EnumTuple(Elem::ET_TRIA3,   "TRIA3"),
       EnumTuple(Elem::ET_TRIA6,   "TRIA6"),
       EnumTuple(Elem::ET_QUAD4,   "QUAD4"),
       EnumTuple(Elem::ET_QUAD8,   "QUAD8"),
       EnumTuple(Elem::ET_QUAD9,   "QUAD9"),
       EnumTuple(Elem::ET_TET4,    "TET4"),
       EnumTuple(Elem::ET_TET10,   "TET10"),
       EnumTuple(Elem::ET_HEXA8,   "HEXA8"),
       EnumTuple(Elem::ET_HEXA20,  "HEXA20"),
       EnumTuple(Elem::ET_HEXA27,  "HEXA27"),
       EnumTuple(Elem::ET_PYRA5,   "PYRA5"),
       EnumTuple(Elem::ET_PYRA13,  "PYRA13"),
       EnumTuple(Elem::ET_WEDGE6,  "WEDGE6"),
       EnumTuple(Elem::ET_WEDGE15, "WEDGE15")      
     };

     Enum<Elem::FEType> Elem::feType = \
     Enum<Elem::FEType>("Finite Element Types",
         sizeof(feTypeTuples) / sizeof(EnumTuple),
         feTypeTuples); 

  
}
