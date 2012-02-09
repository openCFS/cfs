#ifndef FILE_CFS_ELEM_2002
#define FILE_CFS_ELEM_2002

#include <bitset>
#include <map>

#include "General/Enum.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"

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
      ST_POINT  = 1,
      ST_LINE   = 2,
      ST_TRIA   = 3,
      ST_QUAD   = 4,
      ST_TET    = 5,
      ST_HEXA   = 6,
      ST_PYRA   = 7,
      ST_WEDGE  = 8
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
      ET_PYRA14  = 19,
      ET_WEDGE6  = 16,
      ET_WEDGE15 = 17,
      ET_WEDGE18 = 18
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
    
    //! This array contains the global edge numbers of the element 
    //! (absolute value of the entry!), 
    //! In addition, the sign of the edge is encoded in the sign.
    //! In the case of higher order elements, we have to guarantee a globally
    //! unique orientation of edges. By convention, we assume a positive
    //! orientation of an edge, if the nodes have ascending order:
    //! 
    //!   +----+    +-> I = A->B   Global orientation: index(A) < index(B)
    //!   A    B
    //!
    //! To get the correct mapping of the local xi-direction, we can use the
    //! following peace of code:
    //! \code
    //! double xi = dirI;
    //! xi *= edges[i] < 0 ? -1.0 : 1.0;
    //! \code
    //! 
    StdVector<Integer> edges;

    //! Array with face numbers
    StdVector<Integer> faces;

    //! Array with mapping of local(xi/eta) to global(I/II)  orientation

    //! In the case of higher order elements, we have to guarantee a global
    //! orientation of faces, described by the direction I and II. 
    //! By convention, we assume that the global orientation of a face is defined 
    //! by an ascending ordering of the nodal connectivity in such a way:
    //!
    //!  C +----+ C   ^ II    Global orientation:  index(A) < index(B) < index(C)
    //!    |    |     |       
    //!    |    |     +--> I         I-direction: A -> B
    //!  A +----+ B                 II-direction: A -> C                    
    //!         
    //! Based on: 
    //!   Solin, Segeth: "Higher-Order Finite-Element Methods", 2004, p. 164-170
    //! 
    //!
    //! For quadrilateral faces, we thus have 8 possible orientations of the face,
    //! i.e. 8 different possibilities how the local element directions (xi/eta)
    //! relate to the global directions (I/II).  
    //!
    //! To describe them, we introduce for every face a 3-bit array, where the 
    //! flags have the following meaning:
    //!
    //!  2   1   0     indices( 2 = most significant bit)
    //! +-----------+
    //! | 1 | 1 | 1 |  faceFlags(can be also set in integer representation, 
    //! +-----------+            i.e. 111=7)
    //!   |   |   |
    //!   |   |   \-- true  (=1), if  xi = +dirI/dirII
    //!   |   |       false (=0), if  xi = -dirI/dirII
    //!   |   \------ true  (=1), if eta = +dirI/dirII
    //!   |           false (=0), if eta = -dirI/dirII
    //!   \---------- true  (=1), if local and global directions match in same  
    //!                           order (|xi| = |I|,  |eta| = |II|)
    //!               false (=0), if local and global directions match in reverse 
    //!                           order (|xi| = |II|, |eta| = |I|)
    //!
    //! Here are two examples for the orientation (Note: numbers in the interior
    //! denote the orientation of the reference elements, numbers on the outside 
    //! denote the nodal connectivity):
    //!
    //! Case 1
    //! ======
    //!
    //! connect:   1,2,3,4
    //! faceFlags: [111] 
    //!
    //!   4+---------+3             
    //!    |4  ^eta 3|     ^ I      xi  =   I
    //!    |   |     |     |        eta =  II
    //!    |   +->xi |     +-> II   
    //!    |         |              faceFlags[2]=1/true // order: |xi|=|I|, |eta|=|II|
    //!    |1       2|              faceFlags[1]=1/true // eta = eta
    //!   1+---------+2             faceFlags[0]=1/true //  xi = xi
    //!
    //! Case 2
    //! ======
    //!
    //! connect:   4,1,2,3
    //! faceFlags: [0,1,0] 
    //!
    //!   3+---------+2
    //!    |4  ^eta 3|     ^ I      xi  = -II
    //!    |   |     |  II |        eta =   I
    //!    |   +->xi |   <-+   
    //!    |         |              faceFlags[2]=0/false // order: |xi|=|II|, |eta|=|I|
    //!    |1       2|              faceFlags[1]=1/true  // xi = xi 
    //!   4+---------+1             faceFlags[0]=0/false // eta = -eta
    //!
    //!
    //! To get a correct mapping of directions, one can typically use the following
    //! peace of code:
    //! 
    //! \code
    //! double xi = dirI, eta = dirII;
    //! 
    //! // check, if xi and eta have to get interchanged
    //! if( faceFlags[2]) 
    //!     std::swap(xi, eta);
    //! 
    //! // check for sign
    //! xi  = faceFlags[0] ? xi  : -xi;
    //! eta = faceFlags[1] ? eta : -eta;
    //! \endcode  
    //! 
    //! For a triangular face, the concept remains the same. 
    StdVector<std::bitset<3> > faceFlags;

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
    
    //! Return for given FEtype the corresponding ShapeType 
    static Elem::ShapeType GetShapeType( Elem::FEType type );
    
    //! Return for given ShapeType corresponding element shape
    
    //! This method returns for a given shape (ST_QUAD, ST_LINE)
    //! the corresponding element shape. In our case, we always
    //! return the shape to the 1st order elements.
    static ElemShape& GetShape( Elem::ShapeType );
    
    //@}
   
  public:
    
    //! Global collection of reference element shape
    static std::map<Elem::FEType,ElemShape> shapes;
  };

  
  
  //! Description of geometry of reference elements
  
  //! This struct contains the geometric information about a reference element,
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

    //! Coordinates of nodes (outer vector: number of nodes, inner: dim)
    StdVector<StdVector<Double> > nodeCoords; 

    //! Contains for each edge the vertex node numbers
    StdVector<StdVector<UInt> > edgeVertices;
    
    //! Contains for each edge all node numbers
    StdVector<StdVector<UInt> > edgeNodes;
    
    //! Contains for each local direction the edges pointing in this direction
    StdVector<StdVector<UInt> > locDirEdges;
    
    //! Contains for each face the corner node numbers
    StdVector<StdVector<UInt> > faceVertices;

    //! Contains for each face all node numbers
    StdVector<StdVector<UInt> > faceNodes;
    
    //! Contains for each local direction faces pointing in this direction

    //! This array contains for every local direction (fist index) all
    //! faces (first component of the pair) where at least one edge points
    //! in this direction (second component of the pair denotes the face-local
    //! direction component)
    StdVector<StdVector<std::pair<UInt,UInt> > > locDirFaces;

    // ========================================================================
    //  PUBLIC METHODS
    // ========================================================================

    //! Initialize struct (only required once) 
    static void Initialize();
  };
   

} // end of namespace
#endif

