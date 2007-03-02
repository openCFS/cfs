#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001

#include <Domain/grid.hh>
#include "Domain/edgeFace.hh"

#ifdef ADAPTGRID
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#include "Tetrahedron.h"
#include "Octahedron.h"
#include "MultilevelGrid.h"
#include "MeshReader.h"
#include "TetrahedronMeasure.h"
#include "MeshWriter.h"
#endif

namespace CoupledField
{

  // Forward class declarations
  struct Elem;
  struct SurfElem;

  //! Implementation of a simple, one level grid.

  //! This class implements the base class Grid. Is is a simple
  //! grid class, which is able to handle one level of geomatry without
  //! mesh refinement and multilevel elements.
  class GridCFS : public Grid
  {
  public:


    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
  
    //! Constructor 

    //! Standard Constructor 
    GridCFS(); 
  
    //! Destructor
    virtual ~GridCFS();

    //! Trigger mapping of elements' faces

    //! This method calculates global surface numbers and 
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    void MapFaces();
    
    //! Trigger mapping of edges

    //! This method calculates global edge and surface numbers and 
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    void MapEdges();


      virtual void FinishInit();
    
    //@}


    // =======================================================================
    // GENERAL GRID INFORMATION
    // =======================================================================
    //@{ \name General Grid Information

    //! Search and add directivity nodes, if defined in input file,
    //! to list of saved nodes     
    void GetNodesForDirectivity();

    //! Return if grid uses quadratic elements
    bool IsQuadratic() {return isQuadratic_; }

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    UInt GetNumElemOfType( FEType type );

    //! Return dimension of mesh

    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    UInt GetDim();  

    //! Return maximum number of nodes
  
    //! Returns the maximum node number in the finite element grid.
    UInt GetNumNodes();

      virtual void AddNodes(const UInt numNodes);
      

      virtual void SetNodeCoordinate(const UInt inode, const Point & rfPoint);
      
      virtual void SetNodeCoordinate(const UInt numNode, const Vector<Double> & rfPoint);
      

    //! Returns the number of nodes contained in given region
    UInt GetNumNodes( const StdVector<RegionIdType> & regions );

    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName );



      virtual void AddElems(UInt nElems);
      
      virtual void SetElemData(UInt ielem,
                               FEType type,
                               RegionIdType region,
                               const UInt* connect);

      virtual void GetElemData(const UInt ielem,
                               FEType & type,
                               RegionIdType & region,
                               UInt* connect) const;

    //! Returns the total number of elements in the grid
    UInt GetNumElems();
    
    //! Return maximum number of elements 
  
    //! Returns the total number of volume elements in the  grid
    UInt GetNumVolElems();

    //! Returns the total number of volume elements in the  grid
    UInt GetNumSurfElems();
  
    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of 
    UInt GetNumElems( const StdVector<RegionIdType> & regions );
  
    //! Get list with names of all named nodes
    
    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    virtual void GetListNodeNames( StdVector<std::string> & nodeNames);
      virtual void GetListNodeNames( std::vector<std::string> & nodeNames);

    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    virtual void GetListElemNames( StdVector<std::string> & elemNames);
      virtual void GetListElemNames( std::vector<std::string> & elemNames);

    virtual void GetElemsByName( std::vector<Elem*> & elems,
                                 const std::string & elemsName );
    virtual void GetNodesByName( std::vector<UInt> & nodes,
                                 const std::string & nodesName );
      
      virtual void GetElemNumsByName( std::vector<UInt> & elemsNums,
                                      const std::string & elemsName );
      
    //@}


    // =======================================================================
    // NODE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Node Access Functions

    //! Get list of nodes by their name

    //! Returns a list of nodes, which are defined by a name.
    //! These names can specify e.g. boundary nodes or nodes for
    //! saving results.
    //! \param nodeList (out) list with node numbers
    //! \param name (in) name of nodes
    void GetNodesByName( StdVector<UInt> & nodeList,
                         const std::string & name );

    //! Get list of nodes contained in a region

    //! Returns a list of all nodes, which are contained in a 
    //! volume- or surface-region.
    //! \param nodeList (out) list with node numbers
    //! \param regionId (in) region identifier
    void GetNodesByRegion( StdVector<UInt> & nodeList,
                           const RegionIdType regionId );
    
    //! Get coordinates of node with global number inode
    //! \param rfPoint (out) coordinates of point 2D
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    void GetNodeCoordinate( Point & rfPoint,
                            const UInt inode,
                            bool updated);
  

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode,
                            bool updated );
    //@}

    // =======================================================================
    // ELEMENT ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Element Access Functions

    //! Get element with given element number
    
    //! Returns a single element with the given element number
    //! \param elemNr element number
    const Elem * GetElem( UInt elemNr );

      virtual const UInt GetMaxNumNodesPerElem()
      {
          return maxNumElemNodes_;
      }
      
      virtual void SetElemRegion(UInt ielem, RegionIdType region);


    //! Get list of elements (surface / volumes)
    
    //! Returns all elems for a given surface / volume region. If the desired 
    //! region consists of surface elements, they are up-casted into Elem*.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetElems( StdVector<Elem*> & elems,
                   const RegionIdType regionId );

    
    //! Get list of volume elements

    //! Returns all elements for a given volume region.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetVolElems( StdVector<Elem*> & elems, 
                      const RegionIdType regionId );

    //! Get list of surface elements
  
    //! Returns all elements for a given surface region 
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetSurfElems( StdVector<SurfElem*> & elems, 
                       const RegionIdType regionId );

    //! Get list
    void GetElemsByName( StdVector<Elem*> & elems,
                         const std::string & elemsName );

    //! Get node numbers of given element
  
    //! Returns the node numbers of a  given element.
    //! \param connect (out) contains global node numbers
    //! \param iElem (in) element number
    void GetElemNodes( StdVector<UInt> & connect, 
                       const UInt iElem );


      virtual void AddNamedNodes( std::string name, StdVector<UInt> & nodeNums);
      virtual void AddNamedNodes( std::string name, std::vector<UInt> & nodeNums);

      virtual void AddNamedElems( std::string name, StdVector<UInt> & elemNums);
      virtual void AddNamedElems( std::string name, std::vector<UInt> & elemNums);


    //! Get coordinates of element nodes

    //! Returns the coordinates of all nodes of one element.
    //! \param coordMat (out) coordinates of the element nodes 
    //!                         (spaceDim \f$\times\f$ nrNodes);
    //! \param connect (in) global node numbers of element
    //! \param updated (in) flag indicating if updated geometry should be used
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect,
                            bool updated );
  
    //! Get elements associated with given nodes

    //! Returns a list of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes of 
    //!                          nodeList
    //! \param nodeList (in) list of nodes for which neighbouring elements 
    //!                         are needed
    //! \param regionIds (in) identifiers for the regions, where the 
    //!                          neihgbouring elements are searched in
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds );

    //! Get volume elements lying next to given surface elements
  
    //! Returns for a given list of surface elements those (volume) elements, 
    //! which are neighbouring / have a face in common and are within a given 
    //! listof regions. This can be used to determine interfaces.
    //!
    //! \note A surface element is considered to be a neighbour of a volume
    //! element only, if all nodes of the surface element are common with 
    //! a volume element
    //! 
    //! \param neighbours (out) vector of neighbouring elements
    //! \param surfElems (in) vector of surface elements
    //! \param neighRegions (in) region ids, where the volume elems 
    //!                             must lie in
    //!
    //! \note If not each surface element was assigned to EXACTLY ONE volume
    //! element, an error is thrown. If the search was successfull, the
    //! i-th entry in the surfElems-vector corresponds to the i-th
    //! entry in the volElems-vector
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                &neighRegions );
    
    //@}

    // =======================================================================
    // ELEMENT FACE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Surface Access Functions

    //! Get total number of faces in the grid
    UInt GetNumFaces();

    //! Return face with given face
    const Face& GetFace( UInt faceNr);

    //@}

    // =======================================================================
    // ELEMENT EDGE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Edge Access Functions
    
    //! Get number ofe edges
    UInt GetNumEdges();
    
    //! Return edge with given number
    const Edge& GetEdge( UInt edgeNr );

    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    //! Calculates area of a element
    
    //! Calculates the are of an element
    //! \param elem (in) element object
    Double CalcElemArea( const Elem* elem );

    //! Returns surface element normal without defined orientation

    //! This method calculates the normal of an surface element. The direction
    //! is unspecified and dependend on the order or nodes of the element.
    //! However, each surface element stores a flag \a normalSign, which 
    //! indicates the direction of the resulting normal w.r.t. its first
    //! volume element pointer.
    //! \param n (out) vector containing surface normal
    //! \param surfElem (in) reference to surface element
    //! \param updated (in) flag indicating if updated geometry should be used
    void CalcSurfNormal( Vector<Double> & n, 
                         const Elem & surfElem,
                         bool updated );

    //! Returns surface element normal with defined orientation

    //! Calculates the surface normal pointing OUT OF the neighbouring
    //! volume element
    //! \param n (out) normal vector
    //! \param surfElem (in) surface element
    //! \param volElem (in) volume element
    //! \param updated (in) flag indicating if updated geometry should be used

    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem,
                                 bool updated );

    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume. 
    //! 'Volume' here means, that for 2D elements the third dimension is 
    //! assumed to be 1m.
    //! \param regionId (in) region identifier 
    //! \param isaxi (in) flag indicating axial symmetry
    //! \param updated (in) flag indicating if updated geometry should be used
    Double CalcVolumeOfRegion( const RegionIdType regionId ,
                               bool isaxi = false,
                               bool updated = false);
    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous  
 
    //! Returns node numbers of a list of Elements

    //! This method returns the unique node numbers of
    //! a list of given elements. Ther are no duplicate entries.
    //! \param nodeList (out) list of unique node numbers in elemList
    //! \param elemList (in) list of elements
    //! \param onlyLinNodes (in) if true, only the corner nodes are retrieved
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     bool onlyLinNodes = false);
    

    //! Returns the names of all regions

    //! This method return the names of all (surface and volume) regions.
    //! \param regionNames (out) vector containing all region Names
      virtual void GetRegionNames( StdVector<std::string> & regionNames );
      virtual void GetRegionNames( std::vector<std::string> 
                                   & regionNames );

    
    //! Set offset for coordinates due to updated Lagrangian formulation
    void SetNodeOffset( const StdVector<UInt>& nodes, 
                        const Vector<Double>& offsets );

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    bool HasNodalOffset();
    //@}

#ifdef ADAPTGRID
    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Adaptivity Section

    void putNodesFromGrid_RG(grd::MultilevelGrid * grid, const UInt level);
  
    void putElemsFromGrid_RG(grd::MultilevelGrid * grid, const UInt level);
  
    void Refine(grd::MultilevelGrid& grid);
  
    void ReRefine(grd::MultilevelGrid& grid);

    void RefineUniform(grd::MultilevelGrid& grid);

    //@}

#endif


  private:

    // =======================================================================
    // Helper Methods
    // =======================================================================
    //@{ \name Helper Methods

    //! Creates the surface elements

    //! This method creates the surface elements, by assigning each surface 
    //! element one or two volume neighbours. Also the flag for indicating
    //! the direction of the surface normal is calculated.
    //! \param elems (input) vector containing surface elements which are not
    //!                      yet converted to \a SurfElem*
    void CreateSurfaceElements( StdVector<StdVector<Elem*> > & elems);

    
    

    //! Prints information about the grid into the .info file
    void PrintGridInfo() const;
    //@}

    // =======================================================================
    // General attributes
    // =======================================================================
    //@{ \name General Attributes

    //! Flag for initialization status
    bool isInitialized_;

    //! Dimension of grid
    UInt dim_;

    //! Total number of nodes
    UInt numNodes_;

    //! Total number of elements
    UInt numElems_;

    //! Total number of faces
    UInt numFaces_;

    //! Total number of edges
    UInt numEdges_;

    //! Flag indicating if edges are already mapped
    bool edgesMapped_;

    //! Flag indicating if faces are already mapped
    bool facesMapped_;
    
    //! Flag indicating use of quadratic elements
    bool isQuadratic_;

    //@}

    // =======================================================================
    // Mesh attributes
    // =======================================================================
    //@{ \name Mesh Attributes
  
    //! Vector with nodal coordinates
    StdVector<Point> coords_;

    //! Vector with nodal coordinate offsets
    StdVector<Point> deltCoords_;
  
    //! Vector with elements (surface and volume), ordered by element number
    StdVector<Elem*> orderedElems_;
  
    //! Map containing number elements of each type
    std::map<FEType, UInt> numElemTypes_;

      UInt maxNumElemNodes_;
    //@}
  
    //! Map containing face number for each face
    std::map<Face,UInt> faceNums_;

    //! Map containing edge number for each edge
    std::map<Edge,UInt> edgeNums_;

    //! Vector containing all edges
    StdVector<Face> faces_;

    //! Vector containing all edges
    StdVector<Edge> edges_;
    
    // =======================================================================
    // Named Entities
    // =======================================================================
    //@{ \name Named Entities

    //! Vector with named Nodes
    StdVector<StdVector<UInt> > namedNodes_;
  
    //! Vector with names of named nodes
    StdVector<std::string> namedNodeNames_;

    //! Vector with named elements
    StdVector<StdVector<UInt> > namedElems_;
  
    //! Vector with names of named elements
    StdVector<std::string> namedElemNames_;

    //@}
 
    
#ifdef ADAPTGRID
    
    //@{ \name To We Need This ?
    //! procedure for forming list with neighbors
    void FormNeighborsLists();

    //! list with neighbors for element
    StdVector<StdVector<Elem*> >** elNeighbors_;

    //! list with neighbors for nodes
    StdVector<StdVector<Elem*> > vtNeighbors_;

    //! class for mapping Grid_RG and GridCFS
    struct ElementMap {
      int sd;
      StdVector<UInt> map;
    };
  
    StdVector<ElementMap*> elemMap_; //!< mapping between GridRG and GridCFS

    //! only for test
    void SetRefinementFlag();
    //@}
#endif

  };

} // end of namespace
#endif // FILE_GRID_CFS
