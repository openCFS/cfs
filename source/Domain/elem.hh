// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include <bitset>
#include <map>

#include "General/Enum.hh"
#include "Utils/StdVector.hh"
#include "MatVec/vector.hh"
//#include "Utils/tools.hh"

namespace CoupledField
{

  // forward definition
class ElemShape;

  //! Class for description of a volume finite element

  //! This class describes a volume finite element, where volume means the
  //! highest dimensional element entities in the current mesh.
  //! It has to be very lightweight, since this object is created many times.
  //! The finite element is described by:
  //! - corner node numbers
  //! - element number
  //! - element subdomain identifier
  //! - refinement flag / number

  struct Elem {
  
  public:

    //! Dummy constructor
    Elem();

    //! Dummy destructor
    virtual ~Elem() {;}

  public:
    // ========================================================================
    //  Public Enumeration Types
    // ========================================================================
    
    //@{ \name Enumeration types
    
    //! Definition of geometric shapes of elements
    typedef enum 
    {
      ST_UNDEF  = 0, 
      ST_LINE   = 1,
      ST_TRIA   = 2,
      ST_QUAD   = 3,
      ST_TET    = 4,
      ST_HEXA   = 5,
      ST_PYRA   = 6,
      ST_WEDGE  = 7
    } ShapeType;

    //! Static Enum for conversion of ElemShapeType
    static Enum<ShapeType> shapeType;

    //! Definition of supported geometric elements (i.e. Lagrangian elements)
    typedef enum { 
      ET_UNDEF   =  0,
      ET_POINT   =  1,
      ET_LINE2   =  2,
      ET_LINE3   =  3,
      ET_TRIA3   =  4,
      ET_TRIA6   =  5,
      ET_QUAD4   =  6,
      ET_QUAD8   =  7,
      ET_QUAD9   =  8,
      ET_TET4    =  9,
      ET_TET10   = 10,
      ET_HEXA8   = 11,
      ET_HEXA20  = 12,
      ET_HEXA27  = 13,
      ET_PYRA5   = 14,
      ET_PYRA13  = 15,
      ET_WEDGE6  = 16,
      ET_WEDGE15 = 17
    } FEType;

    //! Static Enum for FEType
    static Enum<FEType> feType;

    //@}
    
    // ======================================================
    // GEOMETRICAL INFORMATION
    // ======================================================

    //@{ \name Geometrical Information
    
    //! Global element number
    UInt elemNum;
    
    //! Type of element
    Elem::FEType type;

    //! Identifier for region
    RegionIdType regionId;

    //! Array with node numbers
    StdVector<UInt> connect;

    //! Array with edge numbers
    StdVector<Integer> edges;

    //! Array with face numbers
    StdVector<Integer> faces;

    //! Bitset describing the orientation of the faces (3 for each)
    StdVector<std::bitset<3> > faceFlags;

#ifdef ADAPTGRID
    //! Flag for refinement
    bool refinementFlag;

    //! Number of refinement for the element
    UInt refinementNumber;

#endif
    //@}

    // ======================================================
    // HELPER METHODS
    // ======================================================
    //@{ \name Helper Methods

    //! Overloading operator =
    Elem & operator=(const Elem& t);

    // Fix problems due to negative Jacobian determinants
    void CorrectConnectivity( );
    
    //! Obtain string representation
    std::string ToString() const;
    
    //@}
   
  public:
    
    //! Global collection of reference element shape
    static std::map<Elem::FEType,ElemShape> shapes;
  };

  
  
  //! Description of geometry of reference elements
  
  //! This struct contains the geometric information about a reference element.
  //! i.e. the nodal positions, edge/face indices etc.
  //! A collection of all element shapes can be found in the variable 
  //! Elem::shapes.
  struct ElemShape {

    //! Constructor
    ElemShape();
    
    // ========================================================================
    //  Public Data Members
    // ========================================================================

    //! Dimension of element
    UInt dim;

    //! Order of geometric element (linear, quadratic, cubic)
    UInt order;
    
    //! Number of vertices (corners)
    UInt numVertices;
    
    //! Number of nodes
    UInt numNodes;

    //! Number of edges
    UInt numEdges;

    //! Number of faces
    UInt numFaces;

    //! Coordinate of element midpoint
    Vector<Double> midPointCoord;

    //! Coordinates of nodes
    StdVector<StdVector<Double> > nodeCoords; 

    //! Contains for each edge the vertex node numbers
    StdVector<StdVector<UInt> > edgeVertices;
    
    //! Contains for each edge all node numbers
    StdVector<StdVector<UInt> > edgeNodes;
    
    //! Contains for each face the corner node numbers
    StdVector<StdVector<UInt> > faceVertices;

    //! Contains for each face all node numbers
    StdVector<StdVector<UInt> > faceNodes;

    // ========================================================================
    //  PUBLIC METHODS
    // ========================================================================

    //! Return for given FEtype the corresponding ShapeType 
    Elem::ShapeType GetShapeType( Elem::FEType type ) const;

    //! Initialize struct (only required once) 
    static void Initialize();
  };
   

} // end of namespace
#endif

