#ifndef FILE_CFS_GRID_STRUCT_2004
#define FILE_CFS_GRID_STRUCT_2004

#include "DataInOut/filetype.hh"
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

#include "DataInOut/WriteInfo.hh"
#endif

namespace CoupledField
{


  // Forward class declarations
  struct Elem;
  struct SurfElem;

  /// Class for working with grid
  template<UInt DIM> 
  class GridStruct
  {
  public:
    //! Constructor with parameter - pointer to FileType for reading initial grid
    GridStruct(UInt DIM);

    //! Deconstructor
    virtual ~GridStruct() {}; 

    //! Read Grid Information
    //for normal Grid
    void Read();
    
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

   // =======================================================================
    // SPECIAL METHODS for STRUCTGRID
    // =======================================================================
    //@{ \name General Grid Information

    //! generate the structured grid
    void GenGridStruct(const UInt elemx,const UInt elemy,const UInt elemz);

    //! Transorm Grid
    //void TransformGridStruct(UInt & nodeShift, UInt & shiftFactor, const UInt flag);
      void TransformGridStruct(UInt & shiftFactor, UInt & nodeShift,
			UInt & elemgrid, Double &  meshsize, const UInt flag);
    //! Get maximum elements in "dir"-direction
    Integer GetMaxElem(std::string dir);

    //! Set Dirichlet BCs nodes
    void SetNamedNodes();

    //@}

    //====================================================================================

    // =======================================================================
    // GENERAL GRID INFORMATION
    // =======================================================================
    //@{ \name General Grid Information

    //! Return if grid uses quadratic elements
    bool IsQuadratic() {return false; }

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    UInt GetNumElemOfType( FEType type );

    //! Return dimension of mesh

    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    UInt GetDim() 
    { return dim_; };  

    //! Return maximum number of nodes
  
    //! Returns the maximum node number in the finite element grid.
    UInt GetNumNodes()
    { return numNodes_; };  

    //! Returns the number of nodes contained in given region
    UInt GetNumNodes( const StdVector<RegionIdType> & regions )
    {Error("Method not supported by GridStruct-Class",__FILE__,__LINE__); return 0;};

    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName );

    //! Returns the total number of elements in the grid
    UInt GetNumElems()
    { return numVolElems_ + numSurfElems_;}
    

    //! Return maximum number of elements 
  
    //! Returns the total number of volume elements in the  grid
    UInt GetNumVolElems()
    { return numVolElems_; };  

    //! Returns the total number of volume elements in the  grid
    UInt GetNumSurfElems()
    { return numSurfElems_; };  

  
    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of 
    UInt GetNumElems( const StdVector<RegionIdType> & regions )
    //  { Error("Makes no sense for Slicing technique"); return 0;};
    { return numVolElems_; }; 

    //! Get vector with all (surface and volume) region identifiers
    
    //! Return a vector with names of all region identifiers in the
    //! current mesh (surface and volume)
    //! \param volRegions (out) vector with region identifiers
    void GetRegionIds( StdVector<RegionIdType> & regions );

    
    //! Get vector with all volume region identifiers

    //! Return a vector with names of all volume region identifiers in the
    //! current mesh.
    //! \param volRegions (out) vector with volume region identifiers
    void GetVolRegionIds( StdVector<RegionIdType> & volRegions )
    {volRegions = volRegionIds_;} ;
  
    //! Get vector with all surface region identifiers

    //! Return a vector with names of all surface region identifiers in the
    //! current mesh.
    //! \param surfRegions (out) vector with volume region identifiers
    void GetSurfRegionIds( StdVector<RegionIdType> & surfRegions )
    {surfRegions = surfRegionIds_;} ;

    //! Get list with names of all named nodes
    
    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    void GetListNodeNames( StdVector<std::string> & nodeNames)
    { nodeNames =  namedNodeNames_;}
          
    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    void GetListElemNames( StdVector<std::string> & elemNames) 
    { elemNames =  namedElemNames_;}
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
    void GetNodeCoordinate( Point<DIM> & rfPoint,
                            const UInt inode,
                            bool updated );
  

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
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds )
    {Error("Method not supported by GridStruct-Class",__FILE__,__LINE__); };

    //! Get volume elements lying next to given surface elements
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                &neighRegions )
    {Error("Method not supported by GridStruct-Class",__FILE__,__LINE__); };
   
    //@}


    // =======================================================================
    // ELEMENT FACE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Surface Access Functions

    //! Get total number of faces in the grid
    UInt GetNumFaces() {
      Error( "Not implemented!", __FILE__, __LINE__ );
    }

    //! Return face with given face
    Face& GetFace( UInt faceNr) {
      Error( "Not implemented!", __FILE__, __LINE__ );
    }

    //@}

    // =======================================================================
    // ELEMENT EDGE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Edge Access Functions
          

    //! Get number ofe edges
    UInt GetNumEdges() {
      Error( "Not implemented!", __FILE__, __LINE__ );
    }
    
    Edge& GetEdge( UInt edgeNr ) {
      Error( "Not implemented!", __FILE__, __LINE__ );
    }

    //@}

    
    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    //! Returns surface element normal without defined orientation
    void CalcSurfNormal( Vector<Double> & n, 
                         const Elem & surfElem, bool updated  );

    //! Returns surface element normal with defined orientation
    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem,
                                 bool updated )
    {Error("Method not supported by GridStruct-Class",__FILE__,__LINE__); };

    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume. 
    //! 'Volume' here means, that for 2D elements the third dimension is 
    //! assumed to be 1m.
    //! \param regionId (in) region identifier 
    //! \param isaxi (in) flag indicating axial symmetry
    Double CalcVolumeOfRegion( const RegionIdType regionId,
                               bool isaxi = false,
                               bool updated = false ) {
      Error("Method not supported by GridStruct-Class",__FILE__,__LINE__); 
      return -1.0;
    }
    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous  
 
    //! Returns node numbers of a list of Elements
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     			     bool onlyLinNodes = false)
    {Error("Method not supported by GridStruct-Class"__FILE__,__LINE__); };


    //! Returns the names of all regions

    //! This method return the names of all (surface and volume) regions.
    //! \param regionNames (out) vector containing all region Names
    void GetAllRegionNames( StdVector<std::string> & regionNames );
    //@}

    void SetNodeOffset( const StdVector<UInt>& nodes, 
                        const Vector<Double>& offsets ) {
      Error( "Not implemented for structured grid!",
             __FILE__, __LINE__ );
    }

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    bool HasNodalOffset() {
      return false;
    }


  protected:
 
  private:
    // =======================================================================
    // General attributes
    // =======================================================================
    //@{ \name General Attributes

    bool isInitialized_;

    //! Dimension of grid
    UInt dim_;

    //! Total number of nodes
    UInt numNodes_;

    //! Total number of volume elements
    UInt numVolElems_;

    //! Total number of surface elements
    UInt numSurfElems_;

    //! Total number of elements
    UInt numElems_;

    //@}

    // =======================================================================
    // Mesh attributes
    // =======================================================================
    //@{ \name Mesh Attributes
  
    //! Names of regions
    StdVector<std::string> regionNames_;

    //<! pointer to array with coordinates of all nodes
    Point<DIM> * ptCoordinate_;

    //! Vector with elements for each volume region
    StdVector<StdVector<Elem*> > volElems_;  

    //! Vector with node numbers for each volume region
    StdVector<StdVector<UInt> > volElemNodes_;

    //! Vector with region ids for each volume region
    StdVector<RegionIdType> volRegionIds_;
  
    //! Vector with elements for each surface region
    StdVector<StdVector<SurfElem*> > surfElems_;  

    //! Vector with node numbers for each surface region
    StdVector<StdVector<UInt> > surfElemNodes_;

    //! Vector with region ids for each surface region
    StdVector<RegionIdType> surfRegionIds_;

    //! Vector with elements (surface and volume), ordered by element number
    StdVector<Elem*> orderedElems_;

    //@}
  
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
    

    // =======================================================================
    // Struct Grid specifics
    // =======================================================================
    //@{ \name Named Entities

    //<! dimension of grid
    std::string pdename_;

    //<! maximum number of elements in x direction
    Integer maxnumelemx_;

    //<! maximum number of elements in y direction
    Integer maxnumelemy_;

    //<! maximum number of elements in z direction
    Integer maxnumelemz_;

    Double dimX_;
    Double dimY_;
    Double dimZ_;
    Double waveLength_;
    Double pulseTime_;
    Double soundSpeed_;
    Double safetyRegion_;
    bool ABC;

    Integer elementsPerWavelength_;
    std::string subname_;
    Integer tol_;
    //@}

  };

} // end of namespace
#endif // FILE_GRID_Struct
