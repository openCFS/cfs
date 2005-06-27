#ifndef FILE_CFS_GRID_STRUCT_2004
#define FILE_CFS_GRID_STRUCT_2004

#include "DataInOut/filetype.hh"

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
    ~GridStruct() {}; 

    //! Read Grid Information
    //for normal Grid
    void Read();

   // =======================================================================
    // SPECIAL METHODS for STRUCTGRID
    // =======================================================================
    //@{ \name General Grid Information

    //! generate the structured grid
    void GenGridStruct(const UInt elemx,const UInt elemy,const UInt elemz);

    //! Transorm Grid
    void TransformGridStruct(UInt & nodeShift, UInt & shiftFactor, const UInt flag);

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
    {Error("Method not supported by GridStruct-Class"); return 0;};

    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName );

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
    void GetNodeCoordinate( Point<DIM> & rfPoint,
                            const UInt inode );
  

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode );
    //@}

    // =======================================================================
    // ELEMENT ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Element Access Functions

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
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect );
  
    //! Get elements associated with given nodes
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds )
    {Error("Method not supported by GridStruct-Class"); };

    //! Get volume elements lying next to given surface elements
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                &neighRegions )
    {Error("Method not supported by GridStruct-Class"); };
   
    //@}
    
    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    //! Calculates area of a element
    
    //! Calculates the are of an element
    //! \param elem (in) element object
    Double CalcElemArea( const Elem* elem )
    {Error("Method not supported by GridStruct-Class"); return 0;};

    //! Returns surface element normal without defined orientation
    void CalcSurfNormal( Vector<Double> & n, 
                         const Elem & surfElem )
    {Error("Method not supported by GridStruct-Class"); };

    //! Returns surface element normal with defined orientation
    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem )
    {Error("Method not supported by GridStruct-Class"); };
    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous  
 
    //! Returns node numbers of a list of Elements
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     			     Boolean onlyLinNodes = FALSE)
    {Error("Method not supported by GridStruct-Class"); };


    //! Returns the names of all regions

    //! This method return the names of all (surface and volume) regions.
    //! \param regionNames (out) vector containing all region Names
    void GetAllRegionNames( StdVector<std::string> & regionNames );
    //@}

  protected:
 
  private:
    // =======================================================================
    // General attributes
    // =======================================================================
    //@{ \name General Attributes

    Boolean isInitialized_;

    //! Dimension of grid
    UInt dim_;

    //! Total number of nodes
    UInt numNodes_;

    //! Total number of volume elements
    UInt numVolElems_;

    //! Total number of surface elements
    UInt numSurfElems_;

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
    Integer elementsPerWavelength_;
    std::string subname_;
    Integer tol_;
    Integer sik_;

    //@}

  };

#ifdef __GNUC__
  template class GridStruct<3>;
  template class GridStruct<2>;
#endif

} // end of namespace
#endif // FILE_GRID_Struct
