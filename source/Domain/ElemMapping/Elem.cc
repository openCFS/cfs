#include "Elem.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "General/Enum.hh"

namespace CoupledField {

// static declaration
std::map<Elem::FEType,ElemShape> Elem::shapes;

  Elem::Elem() : 
    elemNum(0),
    type(ET_UNDEF),
    regionId(NO_REGION_ID),
    extended(NULL)
  {}

  Elem::~Elem(){
    delete extended;
  }

//

  std::string Elem::ToString() const {
    std::ostringstream os;
    os << "elemNum=" << elemNum << " region=" << regionId;
    return os.str();
  }
  
  std::string Elem::ToString(const StdVector<Elem*>& vec)
  {
    std::ostringstream os;
    for(unsigned int i = 0; i < vec.GetSize(); i++)
      os << "(" << vec[i]->elemNum << "/" << vec[i]->regionId << ") ";

    return os.str();
  }


 Elem & Elem::operator=(const Elem& t) {
   if (this!=&t) {
     elemNum = t.elemNum;
     type = t.type;
     regionId=t.regionId;
     connect = t.connect;
     if(t.extended){
       extended = new ExtendedElementInfo;
       extended->edges = t.extended->edges;
       extended->faces = t.extended->faces;
       extended->faceFlags = t.extended->faceFlags;
     }
   }
   return *this;
 }

 
 
 void Elem::GetFaceNodes( UInt faceNum, StdVector<UInt>& nodes ) const {
   // check, if face is defined at all
   Integer locFaceIndex = this->extended->faces.Find(faceNum);
   if( locFaceIndex < 0 ) {
     EXCEPTION("Could not find face " << faceNum << " for element #" 
               << this->elemNum );
   }
   // copy nodes, according to the element shape
    StdVector<UInt> faceNodes = shapes[type].faceNodes[locFaceIndex];
    UInt numNodes = faceNodes.GetSize();
    nodes.Resize(numNodes);
   for( UInt i = 0; i < numNodes; ++i ) {
     nodes[i] = connect[faceNodes[i]-1];
   }
 }
 
 void Elem::GetEdgeNodes( UInt edgeNum, StdVector<UInt>& nodes ) const {
   // check, if edge is defined at all
   UInt locEdgeIndex = 0;
   bool found = false;
   while(locEdgeIndex < extended->edges.GetSize() ) {
     if( std::abs(extended->edges[locEdgeIndex]) == static_cast<Integer>(edgeNum) ) {
       found = true;
       break;
     }
     locEdgeIndex++;
   }
   if( !found ) {
     EXCEPTION("Could not find edge " << edgeNum << " for element #" 
               << this->elemNum );
   }
   // copy nodes, according to the element shape
   StdVector<UInt> faceNodes = shapes[type].edgeVertices[locEdgeIndex];
   UInt numNodes = faceNodes.GetSize();
   nodes.Resize(numNodes);
   for( UInt i = 0; i < numNodes; ++i ) {
     nodes[i] = connect[faceNodes[i]-1];
   }
 }
 
 void Elem::CorrectConnectivity( const Grid& grid )
 {
   UInt dummy;

   // This funtion rearanges the connectivity of the element so
   // that the orientation is in such a way that the Jacobian determinant
   // will always be positive.
   switch(type)
   {
   case ET_TRIA6:
     dummy = connect[4];
     connect[4] = connect[5];
     connect[5] = dummy;
     // intentionally no break
   case ET_TRIA3:
     dummy = connect[0];
     connect[0] = connect[1];
     connect[1] = dummy;
     break;
   case ET_QUAD8:
   case ET_QUAD9:
     dummy = connect[4];
     connect[4] = connect[7];
     connect[7] = dummy;
     dummy = connect[5];
     connect[5] = connect[6];
     connect[6] = dummy;
     // intentionally no break
   case ET_QUAD4:
     dummy = connect[1];
     connect[1] = connect[3];
     connect[3] = dummy;
     break;
   case ET_TET10:
     dummy = connect[4];
     connect[4] = connect[7];
     connect[7] = connect[9];
     connect[9] = connect[5];
     connect[5] = dummy;
     dummy = connect[6];
     connect[6] = connect[8];
     connect[8] = dummy;
     // intentionally no break
   case ET_TET4:
     dummy = connect[3];
     connect[3] = connect[2];
     connect[2] = connect[1];
     connect[1] = connect[0];
     connect[0] = dummy;
     break;
   case ET_HEXA20:
     for (int n1=8; n1<=11; n1++) {
       dummy = connect[n1];
       connect[n1] = connect[n1+4];
       connect[n1+4] = dummy;
     }
     // intentionally no break
   case ET_HEXA8:
     for (int n1=0; n1<=3; n1++) {
       dummy = connect[n1];
       connect[n1] = connect[n1+4];
       connect[n1+4] = dummy;
     }
     break;
   case ET_PYRA13:
     dummy = connect[8];
     connect[8] = connect[6];
     connect[6] = dummy;
     dummy = connect[7];
     connect[7] = connect[5];
     connect[5] = dummy;
     dummy = connect[12];
     connect[12] = connect[10];
     connect[10] = dummy;
     // intentionally no break
   case ET_PYRA5:
     dummy = connect[1];
     connect[1] = connect[3];
     connect[3] = dummy;
     break;
   case ET_WEDGE15:
     for (int n1=6; n1<=8; n1++) {
       dummy = connect[n1];
       connect[n1] = connect[n1+3];
       connect[n1+3] = dummy;
     }
     // intentionally no break
   case ET_WEDGE6:
     for (int n1=0; n1<=2; n1++) {
       dummy = connect[n1];
       connect[n1] = connect[n1+3];
       connect[n1+3] = dummy;
     }
     break;
   default:
     EXCEPTION("Connectivity for " << feType.ToString(type) << " element "
               << elemNum << " in region "
               << grid.GetRegion().ToString(regionId)
               << " is not properly oriented!" );
     break;
   }
  }
  

 std::ostream& operator<< ( std::ostream& os , const Elem& elem)
 {
   os << elem.elemNum << " ";
   return os;
 }

  void SetElemInfo( ElemShape& shape,
                      Double midPoint[],
                      Double nodeCoords[], UInt edgeVertices[],
                      UInt numEdgeNodes[], UInt edgeNodes[],
                      Integer edgeLocDirs[],
                      //UInt numEdgesPerDir[],  UInt locDirEdges[],
                      Elem::FEType surfaceETypes[] = NULL,
                      UInt numFaceVertices[] = NULL, UInt faceVertices[] = NULL,
                      UInt numFaceNodes[] = NULL, UInt faceNodes[] = NULL,
                      Integer faceLocDirs[][2] = NULL
                      //UInt numFacesPerDir[] = NULL, UInt locDirFaces[][2] = NULL
                      
                      ) {

      //subelement information
      UInt surfEl = shape.numSurfElems;
      shape.surfElemTypes.Resize(surfEl);
      for(UInt eType = 0; eType < surfEl; ++eType){
        shape.surfElemTypes[eType] = surfaceETypes[eType];
      }

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
          assert(nodeCoords != nullptr);
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
      
      // edge local directions
      pos = 0;
      shape.edgeLocDirs.Resize( shape.numEdges );
      for( UInt iEdge = 0; iEdge < shape.numEdges; iEdge++ ) {
        shape.edgeLocDirs[iEdge] = edgeLocDirs[pos++];
      }
      
//      // locDirEdges
//      pos = 0;
//      shape.locDirEdges.Resize( shape.dim );
//      for( UInt iDim = 0; iDim < shape.dim; iDim++ ) {
//        shape.locDirEdges[iDim].Resize( numEdgesPerDir[iDim]);
//        for( UInt iEdge = 0; iEdge < numEdgesPerDir[iDim]; iEdge++ )
//          shape.locDirEdges[iDim][iEdge] = locDirEdges[pos++];
//      }

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
      
      // local face directions
      if( shape.numFaces > 0 ) {
        pos = 0;
        shape.faceLocDirs.Resize( shape.numFaces );
        for( UInt iFace = 0; iFace < shape.numFaces; iFace++ ) {
          shape.faceLocDirs[iFace][0] = faceLocDirs[pos][0];
          shape.faceLocDirs[iFace][1] = faceLocDirs[pos][1];
          pos++;
        }
      }
      
//      // locDirFaces
//      if( shape.numFaces > 0 ) {
//        pos = 0;
//        shape.locDirFaces.Resize( shape.dim );
//        for( UInt iDim = 0; iDim < shape.dim; iDim++ ) {
//          shape.locDirFaces[iDim].Resize( numFacesPerDir[iDim]);
//          for( UInt iFace = 0; iFace < numFacesPerDir[iDim]; iFace++ ) {
//            shape.locDirFaces[iDim][iFace] = 
//                std::pair<UInt,UInt>(locDirFaces[pos][0], 
//                                     locDirFaces[pos][1]);
//            pos++;
//          }
//        }
//      }
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
    {
      ElemShape s;
      s.dim = 0;
      s.order = 0;
      s.numVertices = 1;
      s.numNodes = 1;
      s.numEdges = 0;
      s.numFaces = 0;
      s.numSurfElems = 0;
      s.volume = 0.0;

      Double midPoint[1] = { 0.0 };
      Double nodeCoords[] = { 0.0 }; // #1
      // no dummy ada needed as the num* is 0
      SetElemInfo( s, midPoint, nodeCoords, nullptr, nullptr,
                   nullptr, nullptr, NULL,
                   NULL, NULL, NULL, NULL );
      Elem::shapes[Elem::ET_POINT] = s;
    }
    
    // ************************************************************************
    // LINE2
    // ************************************************************************
    {
      ElemShape s;
      s.dim = 1;
      s.order = 1;
      s.numVertices = 2;
      s.numNodes = 2;
      s.numEdges = 1;
      s.numFaces = 0;
      s.numSurfElems = 2;
      s.volume = 2.0;
      
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
      Integer edgeLocDirs[] =
      {
       0 // #1
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_POINT,Elem::ET_POINT
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeVertices, edgeLocDirs, surfElems,
                   NULL, NULL, NULL, NULL );
      Elem::shapes[Elem::ET_LINE2] = s;
    }
    
    // ************************************************************************
    // LINE3
    // ************************************************************************
    {
      //                   eta
      //                    ^
      //                    |
      // 1 +----[3]---+ 2   +--> xi

      ElemShape s;
      s.dim = 1;
      s.order = 2;
      s.numVertices = 2;
      s.numNodes = 3;
      s.numEdges = 1;
      s.numFaces = 0;
      s.numSurfElems = 2;
      s.volume = 2.0;

      Double midPoint[1] = {0.0};
      Double nodeCoords[] =
      {
       -1.0, // #1
        1.0,  // #2
        0.0  // #3
      };
      UInt edgeVertices[] =
      {
       1, 2 // #1
      };
      UInt numEdgeNodes[] =
      {
       2 /* former 3*/ // #1
      };
      Integer edgeLocDirs[] =
      {
       0 // #1
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_POINT,Elem::ET_POINT
      };

      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeVertices, edgeLocDirs, surfElems,
                   NULL, NULL, NULL, NULL );
      Elem::shapes[Elem::ET_LINE3] = s;
    }
    
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
    s.numSurfElems = 3;
    s.volume = 0.5;

    Double midPoint[2] = {1.0/3.0, 1.0/3.0};
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
     3, 1, // #3
    };
    UInt numEdgeNodes[] = 
    {
     2, // #1
     2, // #2
     2, // #3
    };
    UInt * edgeNodes  = edgeVertices; // nodes = vertices
    Integer edgeLocDirs[] =
    {
      0, // #1
     -1, // #2
      1  // #3
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
    Integer faceLocDirs[][2] =
    {
     {0,1} // #1
    };

    Elem::FEType surfElems[] =
    {
     Elem::ET_LINE2,Elem::ET_LINE2,Elem::ET_LINE2
    };
    SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                 edgeNodes, edgeLocDirs, surfElems,
                 numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                faceLocDirs );
    Elem::shapes[Elem::ET_TRIA3] = s;
    
    
    // ************************************************************************
    // TRIA6
    // ************************************************************************
    {
      // 3 +
      //   |`\
      //   |  `\
      //  [6]   [5]         eta
      //   |      `\        ^        REFERENCE VOLUME ELEMENT
      //   |        `\      |
      // 1 +----[4]---+ 2   +--> xi

      ElemShape s;
      s.dim = 2;
      s.order = 2;
      s.numVertices = 3;
      s.numNodes = 6;
      s.numEdges = 3;
      s.numFaces = 1;
      s.numSurfElems = 3;
      s.volume = 0.5;

      Double midPoint[2] = {1.0/3.0, 1.0/3.0};
      Double nodeCoords[] =
      { 0.0,  0.0, // #1
        1.0,  0.0, // #2
        0.0,  1.0, // #3
        0.5,  0.0, // #4
        0.5,  0.5, // #5
        0.0,  0.5, // #6
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        3, 1, // #3
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
      };
      UInt edgeNodes[] =
      { 1, 2, 4, // #1
    	2, 3, 5, // #2
    	3, 1, 6  // #3
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
      -1, // #2
       1  // #3
      };
      UInt numFaceVertices[] =
      {
       3 // #1
      };
      UInt faceVertices[] =
      {
       1, 2, 3 // #1
      };
      UInt numFaceNodes[] =
      {
       6, // #1
      };
      UInt faceNodes[] =
      {
       1, 2, 3, 4, 5, 6,   // #1
      };
      Integer faceLocDirs[][2] =
      {
       {0,1} // #1
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_LINE3,Elem::ET_LINE3,Elem::ET_LINE3
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_TRIA6] = s;
    }
    
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
      s.numSurfElems = 4;
      s.volume = 4.0;
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
      Integer edgeLocDirs[] =
      {
        0, // #1
        1, // #2
        0, // #3
        1  // #4
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
      Integer faceLocDirs[][2] =
      {
       {0,1} // #1
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_LINE2,Elem::ET_LINE2,Elem::ET_LINE2,Elem::ET_LINE2
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
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
       s.order = 2;
       s.numVertices = 4;
       s.numNodes = 8;
       s.numEdges = 4;
       s.numFaces = 1;
       s.numSurfElems = 4;
       s.volume = 4.0;
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
       Integer edgeLocDirs[] =
       {
        0, // #1
        1, // #2
        0, // #3
        1  // #4
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
       Integer faceLocDirs[][2] =
       {
        {0,1} // #1
       };
       Elem::FEType surfElems[] =
       {
        Elem::ET_LINE3,Elem::ET_LINE3,Elem::ET_LINE3,Elem::ET_LINE3
       };
       SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
       Elem::shapes[Elem::ET_QUAD8] = s;
     }
    // ************************************************************************
    // QUAD9
    // ************************************************************************
    {
      //      eta
      //       ^
      //       7
      // 4 +---+---+ 3
      //   |   |9  |
      // 8 +   0---+6> xi     REFERENCE VOLUME ELEMENT
      //   |       |
      // 1 +---+---+ 2
      //       5
       ElemShape s;
       s.dim = 2;
       s.order = 2;
       s.numVertices = 4;
       s.numNodes = 9;
       s.numEdges = 4;
       s.numFaces = 1;
       s.numSurfElems = 4;
       s.volume = 4.0;
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
        -1.0,  0.0, // #8
         0.0,  0.0  // #9
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
       Integer edgeLocDirs[] =
       {
        0, // #1
        1, // #2
        0, // #3
        1  // #4
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
        9 // #1
       };
       UInt faceNodes[] =
       {
        1, 2, 3, 4, 5, 6, 7, 8, 9 // #1
       };
       Integer faceLocDirs[][2] =
       {
        {0,1} // #1
       };
       Elem::FEType surfElems[] =
       {
        Elem::ET_LINE3,Elem::ET_LINE3,Elem::ET_LINE3,Elem::ET_LINE3
       };
       SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
       Elem::shapes[Elem::ET_QUAD9] = s;
     }
    
    // ************************************************************************
    // TET4
    // ************************************************************************
    {
      //                      zeta
      //                    .
      //                  ,/
      //                 /
      //              4
      //            ,/|`\
      //          ,/  |  `\
      //        ,/    '.   `\
      //      ,/       |     `\
      //    ,/         |       `\
      //   1-----------'.--------3 --> eta
      //    `\.         |      ,/
      //       `\.      |    ,/
      //          `\.   '. ,/
      //             `\. |/
      //                `2
      //                   `\.
      //                      ` xi
      //

      ElemShape s;
      s.dim = 3;
      s.order = 1;
      s.numVertices = 4;
      s.numNodes = 4;
      s.numEdges = 6;
      s.numFaces = 4;
      s.numSurfElems = 4;
      s.volume = 1./6.0;
      Double midPoint[3] = {1.0/4.0, 1.0/4.0, 1.0/4.0};
      Double nodeCoords[] =
      { 0.0,  0.0,  0.0, // #1
        1.0,  0.0,  0.0, // #2
        0.0,  1.0,  0.0, // #3
        0.0,  0.0,  1.0, // #4
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        1, 3, // #3
        1, 4, // #4
        2, 4, // #5
        3, 4, // #6
      };
      UInt numEdgeNodes[] =
      {
       2, // #1
       2, // #2
       2, // #3
       2, // #4
       2, // #5
       2, // #6
      };
      UInt * edgeNodes = edgeVertices;
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
      -1, // #3
      -1, // #4
      -1, // #5
      -1  // #6
      };
      UInt numFaceVertices[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
      };
      UInt faceVertices[] =
      {
       1, 2, 3,    // #1
       1, 2, 4,    // #2
       1, 3, 4,    // #3
       2, 3, 4,    // #4
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      Integer faceLocDirs[][2] =
      {
       {-1,-1}, // #1
       {-1,-1}, // #2
       {-1,-1}, // #3
       {-1,-1}  // #4
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_TRIA3,Elem::ET_TRIA3,Elem::ET_TRIA3,Elem::ET_TRIA3
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_TET4] = s;
    }
    
    // ************************************************************************
    // TET10
    // ************************************************************************
    {
      //                      zeta
      //                    .
      //                  ,/
      //                 /
      //              4
      //            ,/|`\
      //          ,/  |  `\
      //        [8]   '.  [10]
      //      ,/       |     `\
      //    ,/         |       `\
      //   1--------[7][9]-------3 --> eta
      //    `\.         |      ,/
      //       `\.      |    [6]
      //         [5].   '. ,/
      //             `\. |/
      //                `2
      //                   `\.
      //                      ` xi
      //

      ElemShape s;
      s.dim = 3;
      s.order = 2;
      s.numVertices = 4;
      s.numNodes = 10;
      s.numEdges = 6;
      s.numFaces = 4;
      s.numSurfElems = 4;
      s.volume = 1./6.0;

      Double midPoint[3] = {1.0/4.0, 1.0/4.0, 1.0/4.0};
      Double nodeCoords[] =
      { 0.0,  0.0,  0.0, // #1
        1.0,  0.0,  0.0, // #2
        0.0,  1.0,  0.0, // #3
        0.0,  0.0,  1.0, // #4
        0.5,  0.0,  0.0, // #5
        0.5,  0.5,  0.0, // #6
        0.0,  0.5,  0.0, // #7
        0.0,  0.0,  0.5, // #8
        0.5,  0.0,  0.5, // #9
        0.0,  0.5,  0.5, // #10
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        1, 3, // #3
        1, 4, // #4 
        2, 4, // #5
        3, 4, // #6
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
      };
      UInt edgeNodes[] =
      { 1, 2, 5,  // #1
        2, 3, 6,  // #4
        1, 3, 7,  // #2
        1, 4, 8,  // #3
        2, 4, 9,  // #5
        3, 4, 10  // #6
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
      -1, // #3
      -1, // #4
      -1, // #5
      -1  // #6
      };
      UInt numFaceVertices[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
      };
      UInt faceVertices[] =
      {
       1, 2, 3,    // #1
       1, 2, 4,    // #2
       1, 3, 4,    // #3
       2, 3, 4,    // #4
      };
      UInt numFaceNodes[] =
      {
       6, // #1
       6, // #2
       6, // #3
       6, // #4
      };
      UInt faceNodes[] =
      {
       1, 2, 3, 5, 6, 7,   // #1
       1, 2, 4, 5, 9, 8,   // #2
       1, 3, 4, 7, 10, 8,  // #3
       2, 3, 4, 6, 10, 9   // #4
      };
      Integer faceLocDirs[][2] =
      {
       {-1,-1}, // #1
       {-1,-1}, // #2
       {-1,-1}, // #3
       {-1,-1}  // #4
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_TET10] = s;
    }
    
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
      s.numSurfElems = 6;
      s.volume = 8.0;
      
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
        5, 6, // #5
        6, 7, // #6
        8, 7, // #7
        5, 8, // #8
        1, 5, // #9
        2, 6, // #10
        3, 7, // #11
        4, 8  // #12
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
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
       0, // #5
       1, // #6
       0, // #7
       1, // #8
       2, // #9
       2, // #10
       2, // #11
       2, // #12
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
       1, 2, 6, 5, // #1
       2, 3, 7, 6, // #2
       4, 3, 7, 8, // #3
       1, 4, 8, 5, // #4
       1, 2, 3, 4, // #5
       5, 6, 7, 8  // #6
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      Integer faceLocDirs[][2] =
      {
       {0, 2}, // #1
       {1, 2}, // #2
       {0, 2}, // #3
       {1, 2}, // #4
       {0, 1}, // #5
       {0 ,1}  // #6
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD4,Elem::ET_QUAD4,Elem::ET_QUAD4,
       Elem::ET_QUAD4,Elem::ET_QUAD4,Elem::ET_QUAD4
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
      Elem::shapes[Elem::ET_HEXA8] = s;
    }
    // ************************************************************************
    // HEXA20
    // ************************************************************************
    // ************************************************************************
    {

      //          8+----*----- +7
      //          /|   15    /|           zeta
      //      16 * |      14* |           ^ eta
      //        /20*  13   /  * 19        |/
      //      5+---|-*----+6  |           0-> xi 
      //       |   |   11 |   |
      //       |  4+----*-|---+ 3
      //     17*  /       *18/            REFERENCE VOLUME ELEMENT
      //       | *12      | * 10
      //       |/         |/
      //      1+----*-----+ 2
      //            9
      ElemShape s;
      s.dim = 3;
      s.order = 2;
      s.numVertices = 8;
      s.numNodes = 20;
      s.numEdges = 12;
      s.numFaces = 6;
      s.numSurfElems = 6;
      s.volume = 8.0;

      Double midPoint[3] = {0.0, 0.0, 0.0};
      Double nodeCoords[] = 
      { -1.0, -1.0, -1.0, //  #1
         1.0, -1.0, -1.0, //  #2
         1.0,  1.0, -1.0, //  #3
        -1.0,  1.0, -1.0, //  #4
        -1.0, -1.0,  1.0, //  #5
         1.0, -1.0,  1.0, //  #6
         1.0,  1.0,  1.0, //  #7
        -1.0,  1.0,  1.0, //  #8
         0.0, -1.0, -1.0, //  #9
         1.0,  0.0, -1.0, // #10
         0.0,  1.0, -1.0, // #11
        -1.0,  0.0, -1.0, // #12
         0.0, -1.0,  1.0, // #13
         1.0,  0.0,  1.0, // #14
         0.0,  1.0,  1.0, // #15
        -1.0,  0.0,  1.0, // #16
        -1.0, -1.0,  0.0, // #17
         1.0, -1.0,  0.0, // #18
         1.0,  1.0,  0.0, // #19
        -1.0,  1.0,  0.0  // #20

      };
      UInt edgeVertices[] = 
      { 1, 2, // #1
        2, 3, // #2
        4, 3, // #3
        1, 4, // #4
        5, 6, // #5
        6, 7, // #6
        8, 7, // #7
        5, 8, // #8
        1, 5, // #9
        2, 6, // #10
        3, 7, // #11
        4, 8  // #12
      };
      UInt numEdgeNodes[] = 
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
       3, // #9
       3, // #10
       3, // #11
       3  // #12
      };
      UInt edgeNodes[] =
      { 
       1, 2,  9, // #1
       2, 3, 10, // #2
       4, 3, 11, // #3
       1, 4, 12, // #4
       5, 6, 13, // #5
       6, 7, 14, // #6
       8, 7, 15, // #7
       5, 8, 16, // #8
       1, 5, 17, // #9
       2, 6, 18, // #10
       3, 7, 19, // #11
       4, 8, 20  // #12
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
       0, // #5
       1, // #6
       0, // #7
       1, // #8
       2, // #9
       2, // #10
       2, // #11
       2  // #12
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
       1, 2, 6, 5, // #1
       2, 3, 7, 6, // #2
       4, 3, 7, 8, // #3
       1, 4, 8, 5, // #4
       1, 2, 3, 4, // #5
       5, 6, 7, 8  // #6
      };
      UInt numFaceNodes[] = 
      {
       8, // #1
       8, // #2
       8, // #3
       8, // #4
       8, // #5
       8  // #6
      };
      UInt faceNodes[] = 
      {
       1, 2, 6, 5,  9, 18, 13, 17, // #2
       2, 3, 7, 6, 10, 19, 14, 18, // #3
       4, 3, 7, 8, 11, 19, 15, 20, // #4
       1, 4, 8, 5, 12, 20, 16, 17, // #5
       1, 2, 3, 4,  9, 10, 11, 12, // #1
       5, 6, 7, 8, 13, 14, 15, 16  // #6
      };
      Integer faceLocDirs[][2] =
      {
       {0, 2}, // #2
       {1, 2}, // #3
       {0, 2}, // #4
       {1, 2}, // #5
       {0, 1}, // #1
       {0 ,1}  // #6
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD8,Elem::ET_QUAD8,Elem::ET_QUAD8,Elem::ET_QUAD8,Elem::ET_QUAD8,Elem::ET_QUAD8
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_HEXA20] = s;
    }
    // ************************************************************************
    // HEXA27
    // ************************************************************************
    {
      ElemShape s;
      s.dim = 3;
      s.order = 2;
      s.numVertices = 8;
      s.numNodes = 27;
      s.numEdges = 12;
      s.numFaces = 6;
      s.numSurfElems = 6;
      s.volume = 8.0;
      Double midPoint[3] = {0.0, 0.0, 0.0};
      Double nodeCoords[] =
      { -1.0, -1.0, -1.0, //  #1
         1.0, -1.0, -1.0, //  #2
         1.0,  1.0, -1.0, //  #3
        -1.0,  1.0, -1.0, //  #4
        -1.0, -1.0,  1.0, //  #5
         1.0, -1.0,  1.0, //  #6
         1.0,  1.0,  1.0, //  #7
        -1.0,  1.0,  1.0, //  #8
         0.0, -1.0, -1.0, //  #9
         1.0,  0.0, -1.0, // #10
         0.0,  1.0, -1.0, // #11
        -1.0,  0.0, -1.0, // #12
         0.0, -1.0,  1.0, // #13
         1.0,  0.0,  1.0, // #14
         0.0,  1.0,  1.0, // #15
        -1.0,  0.0,  1.0, // #16
        -1.0, -1.0,  0.0, // #17
         1.0, -1.0,  0.0, // #18
         1.0,  1.0,  0.0, // #19
        -1.0,  1.0,  0.0, // #20

         0.0, -1.0,  0.0, // #21
         1.0,  0.0,  0.0, // #22
         0.0,  1.0,  0.0, // #23
        -1.0,  0.0,  0.0, // #24

         0.0,  0.0, -1.0, // #25
         0.0,  0.0,  1.0, // #26
         0.0,  0.0,  0.0, // #27
         
      };
      UInt edgeVertices[] = 
      { 1, 2, // #1
        2, 3, // #2
        4, 3, // #3
        1, 4, // #4
        5, 6, // #5
        6, 7, // #6
        8, 7, // #7
        5, 8, // #8
        1, 5, // #9
        2, 6, // #10
        3, 7, // #11
        4, 8  // #12
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
       3, // #9
       3, // #10
       3, // #11
       3  // #12
      };
      UInt edgeNodes[] =
      { 
       1, 2,  9, // #1
       2, 3, 10, // #2
       4, 3, 11, // #3
       1, 4, 12, // #4
       5, 6, 13, // #5
       6, 7, 14, // #6
       8, 7, 15, // #7
       5, 8, 16, // #8
       1, 5, 17, // #9
       2, 6, 18, // #10
       3, 7, 19, // #11
       4, 8, 20  // #12
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
       0, // #5
       1, // #6
       0, // #7
       1, // #8
       2, // #9
       2, // #10
       2, // #11
       2  // #12
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
       1, 2, 6, 5, // #1
       2, 3, 7, 6, // #2
       4, 3, 7, 8, // #3
       1, 4, 8, 5, // #4
       1, 2, 3, 4, // #5
       5, 6, 7, 8  // #6
      };
      UInt numFaceNodes[] =
      {
       9, // #1
       9, // #2
       9, // #3
       9, // #4
       9, // #5
       9  // #6
      };
      UInt faceNodes[] =
      {
       1, 2, 6, 5,  9, 18, 13, 17, 21, // #2
       2, 3, 7, 6, 10, 19, 14, 18, 22, // #3
       4, 3, 7, 8, 11, 19, 15, 20, 23, // #4
       1, 4, 8, 5, 12, 20, 16, 17, 24, // #5
       1, 2, 3, 4,  9, 10, 11, 12, 25, // #1
       5, 6, 7, 8, 13, 14, 15, 16, 26  // #6
      };
      Integer faceLocDirs[][2] =
      {
       {0, 2}, // #2
       {1, 2}, // #3
       {0, 2}, // #4
       {1, 2}, // #5
       {0, 1}, // #1
       {0 ,1}  // #6
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD9,Elem::ET_QUAD9,Elem::ET_QUAD9,Elem::ET_QUAD9,Elem::ET_QUAD9,Elem::ET_QUAD9
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_HEXA27] = s;
    }
    // ************************************************************************
    // PYRA5
    // ************************************************************************
    {
      //
      //                 5
      //               ,/|\
      //             ,/ .'|\
      //           ,/   | | \
      //         ,/    .' | `.
      //       ,/      |  '.  \
      //     ,/       .'   |   \            zeta
      //   ,/         |    |    \            ^
      //  3----------.'----2    `.           |
      //   `\        |      `\    \          |
      //     `\     .'        `\   \         +------> eta
      //       `\   |           `\  \        `\
      //         `\.'             `\`          `\
      //            4----------------1           xi
      //
      //
      //
      ElemShape s;
      s.dim = 3;
      s.order = 1;
      s.numVertices = 5;
      s.numNodes = 5;
      s.numEdges = 8;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 4.0/3.0;

      Double midPoint[3] = {0.0, 0.0, 1./4};
      Double nodeCoords[] =
      {  1.0,  1.0,  0.0, // #1
        -1.0,  1.0,  0.0, // #2
        -1.0, -1.0,  0.0, // #3
         1.0, -1.0,  0.0, // #4
         0.0,  0.0,  1.0, // #5
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        3, 4, // #3
        4, 1, // #4
        1, 5, // #5
        2, 5, // #6
        3, 5, // #7
        4, 5, // #8
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
      };
      UInt * edgeNodes = edgeVertices;
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
      -1, // #5
      -1, // #6
      -1, // #7
      -1  // #8
      };
      UInt numFaceVertices[] =
      {
       4, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
      };
      UInt faceVertices[] =
      {
       1, 2, 3, 4, // #1
       1, 2, 5,    // #2
       2, 3, 5,    // #3
       3, 4, 5,    // #4
       4, 1, 5,    // #5
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      Integer faceLocDirs[][2] =
      {
       {0, 1},   // #1
       {-1, -1}, // #2
       {-1, -1}, // #3
       {-1, -1}, // #4
       {-1, -1}  // #5
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD4,Elem::ET_TRIA3,Elem::ET_TRIA3,Elem::ET_TRIA3,Elem::ET_TRIA3
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
      Elem::shapes[Elem::ET_PYRA5] = s;
    }

    // ************************************************************************
    // PYRA13
    // ************************************************************************
    {
      //
      //                 5
      //               ,/|\
      //             ,/ .'|\
      //           ,/   | | \
      //        [12]    .' | `.
      //       ,/      | [11] \
      //     ,/      [13]  | [10]            zeta
      //   ,/         |    |    \            ^
      //  3-------[7]-.'---2    `.           |
      //   `\        |      `\    \          |
      //     [8]    .'        [6]  \         +------> eta
      //       `\   |           `\  \        `\
      //         `\.'             `\`          `\
      //            4-------[9]------1           xi
      //
      //
      //
      ElemShape s;
      s.dim = 3;
      s.order = 2;
      s.numVertices = 5;
      s.numNodes = 13;
      s.numEdges = 8;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 4.0/3.0;
      Double midPoint[3] = {0.0, 0.0, 1./4};
      Double nodeCoords[] =
      {  1.0,  1.0,  0.0, // #1
        -1.0,  1.0,  0.0, // #2
        -1.0, -1.0,  0.0, // #3
         1.0, -1.0,  0.0, // #4
         0.0,  0.0,  1.0, // #5
         0.0,  1.0,  0.0, // #6
        -1.0,  0.0,  0.0, // #7
         0.0, -1.0,  0.0, // #8
         1.0,  0.0,  0.0, // #9
         0.5,  0.5,  0.5, // #10
        -0.5,  0.5,  0.5, // #11
        -0.5, -0.5,  0.5, // #12
         0.5, -0.5,  0.5, // #13
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        4, 3, // #3
        1, 4, // #4
        1, 5, // #5
        2, 5, // #6
        3, 5, // #7
        4, 5, // #8
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
      };
      UInt edgeNodes[] =
      { 1, 2, 6, // #1
        2, 3, 7, // #2
        4, 3, 8, // #3
        1, 4, 9, // #4
        1, 5, 10, // #5
        2, 5, 11, // #6
        3, 5, 12, // #7
        4, 5, 13, // #8
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
      -1, // #5
      -1, // #6
      -1, // #7
      -1  // #8
      };
      UInt numFaceVertices[] =
      {
       4, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
      };
      UInt faceVertices[] =
      {
       1, 2, 3, 4, // #1
       1, 2, 5,    // #2
       2, 3, 5,    // #3
       3, 4, 5,    // #4
       4, 1, 5,    // #5
      };
      UInt numFaceNodes[] =
      {
       8, // #1
       6, // #2
       6, // #3
       6, // #4
       6, // #5
      };
      UInt faceNodes[] =
      {
       1, 2, 3, 4, 6, 7, 8, 9, // #1
       1, 2, 5, 6, 11, 10,     // #2
       2, 3, 5, 7, 12, 11,     // #3
       3, 4, 5, 8, 13, 12,     // #4
       4, 1, 5, 9, 10, 13      // #5
      };
      Integer faceLocDirs[][2] =
      {
       { 0,  1},   // #1
       {-1, -1}, // #2
       {-1, -1}, // #3
       {-1, -1}, // #4
       {-1, -1}  // #5
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD8,Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
      Elem::shapes[Elem::ET_PYRA13] = s;
    }

    // ************************************************************************
    // PYRA14
    // ************************************************************************
    {
      //
      //                 5
      //               ,/|\
      //             ,/ .'|\
      //           ,/   | | \
      //        [12]    .' | `.
      //       ,/      | [11] \
      //     ,/      [13]  | [10]            zeta
      //   ,/         |    |    \            ^
      //  3-------[7]-.'---2    `.           |
      //   `\        |      `\    \          |
      //     [8]    .'[14]    [6]  \         +------> eta
      //       `\   |           `\  \        `\
      //         `\.'             `\`          `\
      //            4-------[9]------1           xi
      //
      //
      //
      ElemShape s;
      s.dim = 3;
      s.order = 2;
      s.numVertices = 5;
      s.numNodes = 14;
      s.numEdges = 8;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 4.0/3.0;
      Double midPoint[3] = {0.0, 0.0, 1.0/4.0};
      Double nodeCoords[] =
      {  1.0,  1.0,  0.0, // #1
        -1.0,  1.0,  0.0, // #2
        -1.0, -1.0,  0.0, // #3
         1.0, -1.0,  0.0, // #4
         0.0,  0.0,  1.0, // #5
         0.0,  1.0,  0.0, // #6
        -1.0,  0.0,  0.0, // #7
         0.0, -1.0,  0.0, // #8
         1.0,  0.0,  0.0, // #9
         0.5,  0.5,  0.5, // #10
        -0.5,  0.5,  0.5, // #11
        -0.5, -0.5,  0.5, // #12
         0.5, -0.5,  0.5, // #13
         0.0,  0.0,  0.0, // #14
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        4, 3, // #3
        1, 4, // #4
        1, 5, // #5
        2, 5, // #6
        3, 5, // #7
        4, 5, // #8
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
      };
      UInt edgeNodes[] =
      { 1, 2, 6, // #1
        2, 3, 7, // #2
        4, 3, 8, // #3
        1, 4, 9, // #4
        1, 5, 10, // #5
        2, 5, 11, // #6
        3, 5, 12, // #7
        4, 5, 13, // #8
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
       1, // #2
       0, // #3
       1, // #4
      -1, // #5
      -1, // #6
      -1, // #7
      -1  // #8
      };
      UInt numFaceVertices[] =
      {
       4, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
      };
      UInt faceVertices[] =
      {
       1, 2, 3, 4, // #1
       1, 2, 5,    // #2
       2, 3, 5,    // #3
       3, 4, 5,    // #4
       4, 1, 5,    // #5
      };
      UInt numFaceNodes[] =
      {
       9, // #1
       6, // #2
       6, // #3
       6, // #4
       6, // #5
      };
      UInt faceNodes[] =
      {
       1, 2, 3, 4, 6, 7, 8, 9, 14, // #1
       1, 2, 5, 6, 11, 10,     // #2
       2, 3, 5, 7, 12, 11,     // #3
       3, 4, 5, 8, 13, 12,     // #4
       4, 1, 5, 9, 10, 13      // #5
      };
      Integer faceLocDirs[][2] =
      {
       {0, 1},   // #1
       {-1, -1}, // #2
       {-1, -1}, // #3
       {-1, -1}, // #4
       {-1, -1}  // #5
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_QUAD9,Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_TRIA6
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                    edgeNodes, edgeLocDirs, surfElems,
                    numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                    faceLocDirs );
      Elem::shapes[Elem::ET_PYRA14] = s;
    }

    // ************************************************************************
    // WEDGE6
    // ************************************************************************
    {
      //      + 6    
      //     /|\
      //    / |  \           zeta
      // 4 +----- + 5         ^  eta
      //   |  + 3 |           |/ 
      //   | / \  |           0--> xi
      //   |/    \|   
      // 1 +------+ 2

      ElemShape s;
      s.dim = 3;
      s.order = 1;
      s.numVertices = 6;
      s.numNodes = 6;
      s.numEdges = 9;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 1.0;
      Double midPoint[3] = {1.0/3.0, 1.0/3.0, 0.0};
      Double nodeCoords[] = 
      { 0.0,  0.0, -1.0, // #1
        1.0,  0.0, -1.0, // #2
        0.0,  1.0, -1.0, // #3
        0.0,  0.0,  1.0, // #4
        1.0,  0.0,  1.0, // #5
        0.0,  1.0,  1.0  // #6
      };
      UInt edgeVertices[] = 
      { 1, 2, // #1
        2, 3, // #2
        3, 1, // #3
        4, 5, // #4
        5, 6, // #5
        6, 4, // #6
        1, 4, // #7
        2, 5, // #8
        3, 6  // #9
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
       2  // #9
      };
      UInt * edgeNodes = edgeVertices;
      Integer edgeLocDirs[] =
       {
        0, // #1
       -1, // #2
        1, // #3
        0, // #4
       -1, // #5
        1, // #6
        2, // #7
        2, // #8
        2  // #9
       };
      UInt numFaceVertices[] = 
      {
       3, // #1
       3, // #2
       4, // #3
       4, // #4
       4, // #5
      };
      UInt faceVertices[] = 
      {
       1, 2, 3,    // #1
       4, 5, 6,    // #2
       1, 2, 5, 4, // #3
       2, 3, 6, 5, // #4
       3, 1, 4, 6, // #5
      };
      UInt * numFaceNodes = numFaceVertices;
      UInt * faceNodes = faceVertices;
      Integer faceLocDirs[][2] =
      {
       { 0, 1}, // #1
       { 0, 1}, // #2
       { 0, 2}, // #3
       {-1, 2}, // #4
       { 1, 2} // #5
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_TRIA3,Elem::ET_TRIA3,Elem::ET_QUAD4,Elem::ET_QUAD4,Elem::ET_QUAD4
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_WEDGE6] = s;
    }
    
    // ************************************************************************
    // WEDGE15
    // ************************************************************************
    {
      //        + 6
      //       /|\
      //      / | \
      //  12 *  |  * 11
      //    /   *15 \        zeta
      // 4 +----*----+ 5      ^  eta
      //   |  10|    |        |/ 
      //   |    + 3  |        0--> xi
      //   |   / \   |
      // 13*  /   \  * 14
      //   | *9   8* |
      //   |/       \|
      // 1 +----*----+ 2
      //        7
      ElemShape s;
      s.dim = 3;
      s.order = 1;
      s.numVertices = 6;
      s.numNodes = 15;
      s.numEdges = 9;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 1.0;

      Double midPoint[3] = {1.0/3.0, 1.0/3.0, 0.0};
      Double nodeCoords[] = 
      { 0.0,  0.0, -1.0, //  #1
        1.0,  0.0, -1.0, //  #2
        0.0,  1.0, -1.0, //  #3
        0.0,  0.0,  1.0, //  #4
        1.0,  0.0,  1.0, //  #5
        0.0,  1.0,  1.0, //  #6
        0.5,  0.0, -1.0, //  #7
        0.5,  0.5, -1.0, //  #8
        0.0,  0.5, -1.0, //  #9
        0.5,  0.0,  1.0, // #10
        0.5,  0.5,  1.0, // #11
        0.0,  0.5,  1.0, // #12
        0.0,  0.0,  0.0, // #13
        1.0,  0.0,  0.0, // #14
        0.0,  1.0,  0.0, // #15

      };
      UInt edgeVertices[] = 
      { 1, 2, // #1
        2, 3, // #2
        3, 1, // #3
        4, 5, // #4
        5, 6, // #5
        6, 4, // #6
        1, 4, // #7
        2, 5, // #8
        3, 6  // #9
      };
      UInt numEdgeNodes[] = 
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
       3  // #9
      };
      UInt edgeNodes[] = 
      { 1, 2,  7, // #1
        2, 3,  8, // #2
        3, 1,  9, // #3
        4, 5, 10, // #4
        5, 6, 11, // #5
        6, 4, 12, // #6
        1, 4, 13, // #7
        2, 5, 14, // #8
        3, 6, 15  // #9
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
      -1, // #2
       1, // #3
       0, // #4
      -1, // #5
       1, // #6
       2, // #7
       2, // #8
       2  // #9
      };
      UInt numFaceVertices[] = 
      {
       3, // #1
       3, // #2
       4, // #3
       4, // #4
       4, // #5
      };
      UInt faceVertices[] = 
      {
       1, 2, 3,    // #1
       4, 5, 6,    // #2
       1, 2, 5, 4, // #3
       2, 3, 6, 5, // #4
       3, 1, 4, 6, // #5
      };
      UInt numFaceNodes[] = 
      {
       6, // #1
       6, // #2
       8, // #3
       8, // #4
       8, // #5
      };
      UInt faceNodes[] = 
      {
       1, 2, 3,     7,  8,  9,     // #1
       4, 5, 6,    10, 11, 12,     // #2
       1, 2, 5, 4,  7, 14, 10, 13, // #3
       2, 3, 6, 5,  8, 15, 11, 14, // #4
       3, 1, 4, 6,  9, 13, 12, 15  // #5
      };
      Integer faceLocDirs[][2] =
       {
        { 0, 1}, // #1
        { 0, 1}, // #2
        { 0, 2}, // #3
        {-1, 2}, // #4
        { 1, 2} // #5
       };
      Elem::FEType surfElems[] =
      {
       Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_QUAD8,Elem::ET_QUAD8,Elem::ET_QUAD8
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_WEDGE15] = s;
    }

    // ************************************************************************
    // WEDGE18
    // ************************************************************************
    {
      //        + 6
      //       /|\
      //      / | \
      //  12 *  |  * 11
      //    /   *15 \        zeta
      // 4 +----*----+ 5      ^  eta
      //   |  10| 17 |        |/
      //   |18  + 3  |        0--> xi
      //   |   / \   |
      // 13*  /16 \  * 14
      //   | *9   8* |
      //   |/       \|
      // 1 +----*----+ 2
      //        7
      ElemShape s;
      s.dim = 3;
      s.order = 1;
      s.numVertices = 6;
      s.numNodes = 18;
      s.numEdges = 9;
      s.numFaces = 5;
      s.numSurfElems = 5;
      s.volume = 1.0;
      Double midPoint[3] = {1.0/3.0, 1.0/3.0, 0.0};
      Double nodeCoords[] =
      { 0.0,  0.0, -1.0, //  #1
        1.0,  0.0, -1.0, //  #2
        0.0,  1.0, -1.0, //  #3
        0.0,  0.0,  1.0, //  #4
        1.0,  0.0,  1.0, //  #5
        0.0,  1.0,  1.0, //  #6
        0.5,  0.0, -1.0, //  #7
        0.5,  0.5, -1.0, //  #8
        0.0,  0.5, -1.0, //  #9
        0.5,  0.0,  1.0, // #10
        0.5,  0.5,  1.0, // #11
        0.0,  0.5,  1.0, // #12
        0.0,  0.0,  0.0, // #13
        1.0,  0.0,  0.0, // #14
        0.0,  1.0,  0.0, // #15
        0.5,  0.0,  0.0, // #16
        0.5,  0.5,  0.0, // #17
        0.0,  0.5,  0.0, // #18
      };
      UInt edgeVertices[] =
      { 1, 2, // #1
        2, 3, // #2
        3, 1, // #3
        4, 5, // #4
        5, 6, // #5
        6, 4, // #6
        1, 4, // #7
        2, 5, // #8
        3, 6  // #9
      };
      UInt numEdgeNodes[] =
      {
       3, // #1
       3, // #2
       3, // #3
       3, // #4
       3, // #5
       3, // #6
       3, // #7
       3, // #8
       3  // #9
      };
      UInt edgeNodes[] =
      { 1, 2,  7, // #1
        2, 3,  8, // #2
        3, 1,  9, // #3
        4, 5, 10, // #4
        5, 6, 11, // #5
        6, 4, 12, // #6
        1, 4, 13, // #7
        2, 5, 14, // #8
        3, 6, 15  // #9
      };
      Integer edgeLocDirs[] =
      {
       0, // #1
      -1, // #2
       1, // #3
       0, // #4
      -1, // #5
       1, // #6
       2, // #7
       2, // #8
       2  // #9
      };
      UInt numFaceVertices[] =
      {
       3, // #1
       3, // #2
       4, // #3
       4, // #4
       4, // #5
      };
      UInt faceVertices[] =
      {
       1, 2, 3,    // #1
       4, 5, 6,    // #2
       1, 2, 5, 4, // #3
       2, 3, 6, 5, // #4
       3, 1, 4, 6, // #5
      };
      UInt numFaceNodes[] =
      {
       6, // #1
       6, // #2
       9, // #3
       9, // #4
       9, // #5
      };
      UInt faceNodes[] =
      {
       1, 2, 3,     7,  8,  9,     // #1
       4, 5, 6,    10, 11, 12,     // #2
       1, 2, 5, 4,  7, 14, 10, 13, 16, // #3
       2, 3, 6, 5,  8, 15, 11, 14, 17, // #4
       3, 1, 4, 6,  9, 13, 12, 15, 18  // #5
      };
      Integer faceLocDirs[][2] =
      {
       { 0, 1}, // #1
       { 0, 1}, // #2
       { 0, 2}, // #3
       {-1, 2}, // #4
       { 1, 2} // #5
      };
      Elem::FEType surfElems[] =
      {
       Elem::ET_TRIA6,Elem::ET_TRIA6,Elem::ET_QUAD9,Elem::ET_QUAD9,Elem::ET_QUAD9
      };
      SetElemInfo( s, midPoint, nodeCoords, edgeVertices, numEdgeNodes,
                   edgeNodes, edgeLocDirs, surfElems,
                   numFaceVertices, faceVertices, numFaceNodes, faceNodes,
                   faceLocDirs );
      Elem::shapes[Elem::ET_WEDGE18] = s;
    }
    } 
    
      
    Elem::ShapeType Elem::GetShapeType( Elem::FEType type ) {
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
        case ET_PYRA14:
          ret = ST_PYRA;
          break;
        case ET_WEDGE6:
        case ET_WEDGE15:
        case ET_WEDGE18:
          ret = ST_WEDGE;
          break;
        case ET_POLYGON:
          ret = ST_POLYGON;
          break;
        case ET_POLYHEDRON:
          ret = ST_POLYHEDRON;
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
        case ST_POLYGON:
          feType = ET_POLYGON;
          break;
        case ST_POLYHEDRON:
          feType = ET_POLYHEDRON;
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
       EnumTuple(Elem::ST_POINT,  "POINT"),
       EnumTuple(Elem::ST_LINE,   "LINE"),
       EnumTuple(Elem::ST_TRIA,   "TRIA"),
       EnumTuple(Elem::ST_QUAD,   "QUAD"),
       EnumTuple(Elem::ST_TET,    "TET"),
       EnumTuple(Elem::ST_HEXA,   "HEX"),
       EnumTuple(Elem::ST_PYRA,   "PYRA"),
       EnumTuple(Elem::ST_WEDGE,  "WEDGE"),
       EnumTuple(Elem::ST_POLYGON,  "POLYGON"),
       EnumTuple(Elem::ST_POLYHEDRON,  "POLYHEDRON")
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
       EnumTuple(Elem::ET_PYRA14,  "PYRA14"),
       EnumTuple(Elem::ET_WEDGE6,  "WEDGE6"),
       EnumTuple(Elem::ET_WEDGE15, "WEDGE15"),
       EnumTuple(Elem::ET_WEDGE18, "WEDGE18"),
       EnumTuple(Elem::ET_POLYGON, "POLYGON"),
       EnumTuple(Elem::ET_POLYHEDRON, "POLYHEDRON")
     };

     Enum<Elem::FEType> Elem::feType = \
     Enum<Elem::FEType>("Finite Element Types",
         sizeof(feTypeTuples) / sizeof(EnumTuple),
         feTypeTuples); 

  
}
