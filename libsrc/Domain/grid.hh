#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include <list>

#include "Domain/elem.hh"
#include "Domain/surfElem.hh"

namespace CoupledField
{

  // Forward class declarations
  class FileType;

  //! Class representing geometrical entities (elements, nodes, ...) of a
  //! FE simulation.

  //! This class represents an abstract interface to the geometrical entities 
  //! (elements, nodes, ...) in a FE simulation.
  //! It allows access to elements, nodes and coordinates.
  //! The underlying structure is hidden by the according interface classes.
  //!
  //! \note Some General Remarks:
  //! - Node and element numbers are always counted one-based
  //! - The Type \c regionIdType is guraranteed to be of a countable type,
  //! e.g. \c Integer or \c short \c int. Therefore it is always zero-based
  //! and can be directly used as a index for vectors / arrays.
  //! - The region identifiers are used for surface AND volume elements

  class Grid {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
  
    //! Constructor 

    //! Standard Constructor 
    //! \param aptFileType pointer to FileType for reading initial grid
    Grid(FileType * aptFileType); 
  
    //! Destructor
    virtual ~Grid();

    //! Reads the grid from input file
    virtual void Read()=0;

    //@}

    // =======================================================================
    // GENERAL GRID INFORMATION
    // =======================================================================
    //@{ \name General Grid Information

    //! Return if grid uses quadratic elements
    virtual Boolean IsQuadratic() = 0;

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    virtual UInt GetNumElemOfType( FEType type ) = 0;

    //! Return dimension of mesh

    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    virtual UInt GetDim() = 0;  

    //! Return maximum number of nodes
  
    //! Returns the maximum node number in the finite element grid.
    virtual UInt GetNumNodes() = 0;

    //! Returns the number of nodes contained in given region
    virtual UInt GetNumNodes( const StdVector<RegionIdType> 
                                 & regions ) = 0;

    //! Returns the number of nodes in the given nodelist
    virtual UInt GetNumNodes( const std::string & nodesName ) = 0;

    //! Returns the total number of elements in the grid
    virtual UInt GetNumElems() = 0;

    //! Return maximum number of volume elements 
  
    //! Returns the total number of volume elements in the finite element grid
    virtual UInt GetNumVolElems() = 0;
  
    //! Returns the total number of surface elements in the grid 
    virtual UInt GetNumSurfElems() = 0;

    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of 
    virtual UInt GetNumElems( const StdVector<RegionIdType> 
                                 & regions ) = 0;
  
    //! Get vector with all (surface and volume) region identifiers
    
    //! Return a vector with names of all region identifiers in the
    //! current mesh (surface and volume)
    //! \param volRegions (out) vector with region identifiers
    virtual void GetRegionIds( StdVector<RegionIdType> & regions ) = 0;

    //! Get vector with all volume region identifiers

    //! Return a vector with names of all volume region identifiers in the
    //! current mesh.
    //! \param volRegions (out) vector with volume region identifiers
    virtual void GetVolRegionIds( StdVector<RegionIdType> & volRegions ) = 0;
  
    //! Get vector with all surface region identifiers

    //! Return a vector with names of all surface region identifiers in the
    //! current mesh.
    //! \param surfRegions (out) vector with volume region identifiers
    virtual void GetSurfRegionIds( StdVector<RegionIdType> 
                                   & surfRegions ) = 0;

    //! Get list with names of all named nodes

    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    virtual void GetListNodeNames( StdVector<std::string> & nodeNames) = 0;

    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    virtual void GetListElemNames( StdVector<std::string> & elemNames) = 0;

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
    virtual void GetNodesByName( StdVector<UInt> & nodeList,
                                 const std::string & name ) = 0;

    //! Get list of nodes contained in a region

    //! Returns a list of all nodes, which are contained in a 
    //! volume- or surface-region.
    //! \param nodeList (out) list with node numbers
    //! \param regionId (in) region identifier 
    virtual void GetNodesByRegion( StdVector<UInt> & nodeList,
                                   const RegionIdType regionId ) = 0;
    
    //! Get coordinates of node with global number inode
    //! \param rfPoint (out) coordinates of point 2D
    //! \param inode (in) node number
    virtual void GetNodeCoordinate( Point<2> & rfPoint,
                                    const UInt inode )
    { Error( "Not implemented", __FILE__, __LINE__ ); }  

    //! Get coordinates of node with global number inode
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    virtual void GetNodeCoordinate( Point<3> & rfPoint,
                                    const UInt inode )
    { Error( "Not implemented", __FILE__, __LINE__ ); }  

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    virtual void GetNodeCoordinate( Vector<Double> & rfPoint,
                                    const UInt inode ) = 0;
  
    //@}

    // =======================================================================
    // ELEMENT ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Element Access Functions

    //! Get element with given element number

    //! Returns a single element with the given element number
    //! \param elemNr element number
    virtual const Elem * GetElem( UInt elemNr ) = 0;


    //! Get list of elements (surface / volumes)
    
    //! Returns all elems for a given surface / volume region. If the desired 
    //! region consists of surface elements, they are up-casted into Elem*.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetElems( StdVector<Elem*> & elems,
                           const RegionIdType regionId ) = 0;

    //! Get list of volume elements

    //! Returns all elements for a given volume region.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetVolElems( StdVector<Elem*> & elems, 
                              const RegionIdType regionId ) = 0;

    //! Get list of surface elements
  
    //! Returns all elements for a given surface region 
    //! \param surfElems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetSurfElems( StdVector<SurfElem*> & surfElems, 
                               const RegionIdType regionId ) = 0;
    
    //! Get list of elements by their name
    virtual void GetElemsByName( StdVector<Elem*> & elems,
                                 const std::string & elemsName ) = 0;

    //! Get node numbers of given element
  
    //! Returns the node numbers of a  given element.
    //! \param connect (out) contains global node numbers
    //! \param iElem (in) element number
    virtual void GetElemNodes( StdVector<UInt> & connect, 
                               const UInt iElem ) = 0;


    //! Get coordinates of element nodes

    //! Returns the coordinates of all nodes of one element.
    //! \param coordMat (out) coordinates of the element nodes 
    //!                         (spaceDim \f$\times\f$ nrNodes);
    //! \param connect (in) global node numbers of element
    virtual void GetElemNodesCoord( Matrix<Double> & coordMat,  
                                    const StdVector<UInt> & connect ) = 0;
  
    //! Get elements associated with given nodes

    //! Returns a list of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes 
    //!                          of nodeList
    //! \param nodeList (in) list of nodes for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the 
    //!                       neihgbouring elements are searched in
    virtual void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                                      const StdVector<UInt> & nodeList,
                                      const StdVector<RegionIdType> 
                                      & regionIds ) = 0;

    //! Get volume elements lying next to given surface elements
  
    //! Returns for a given list of surface elements those (volume) elements, 
    //! which are neighbouring / have a face in common and are within a given 
    //! listof regions. This can be used to determine interfaces.
    //! \note A surface element is considered to be a neighbour of a volume
    //! element only, if all nodes of the surface element are common with 
    //! a volume element
    //! \param neighbours (out) vector of neighbouring elements
    //! \param surfElems (in) vector of surface elements
    //! \param neighRegions (in) region ids, where the volume elems must lie
    //! \note If not each surface element was assigned to EXACTLY ONE volume
    //! element, an error is thrown. If the search was successfull, the
    //! i-th entry in the surfElems-vector corresponds to the i-th
    //! entry in the volElems-vector
    virtual void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                        const StdVector<Elem*> & surfElems,
                                        const StdVector<RegionIdType> 
                                        &neighRegions ) = 0;
    
    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation

    //! Returns surface element normal without defined orientation

    //! This method calculates the normal of an surface element. The direction
    //! is unspecified and dependend on the order or nodes of the element.
    //! However, each surface element stores a flag \a normalSign, which 
    //! indicates the direction of the resulting normal w.r.t. its first
    //! volume element pointer.
    //! \param n (out) vector containing surface normal
    //! \param surfElem (in) reference to surface element
    virtual void CalcSurfNormal( Vector<Double> & n, 
                                 const Elem & surfElem ) = 0;

    //! Returns surface element normal with defined orientation

    //! Calculates the surface normal pointing OUT OF the neighbouring
    //! volume element
    //! \param n (out) normal vector
    //! \param surfElem (in) surface element
    //! \param volElem (in) volume element
    virtual void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                         const Elem & surfElem,
                                         const Elem & volElem ) = 0;

    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume. 
    //! 'Volume' here means, that for 2D elements the third dimension is 
    //! assumed to be 1m.
    //! \param regionId (in) region identifier 
    //! \param isaxi (in) flag indicating axial symmetry
    virtual Double CalcVolumeOfRegion( const RegionIdType regionId,
                                       Boolean isaxi = FALSE ) = 0;
    //@}


    // =======================================================================
    // ADAPTIVITY SECTION
    // =======================================================================
    //@{ \name Mesh Adaptivity

    //! Prolongation of solution
    //! \param sol_coarse (in) solution on coarse grid
    //! \param sol (out) contains the solution on the new grid (fine grid)
    //! \param alevel (in) index in multilevel hierarchy
    virtual void ProlongSol( const Vector<Double> sol_coarse, 
                             Vector<Double> & sol,
                             const UInt alevel )
    { Error( "Not implemented", __FILE__, __LINE__); }

    //! Do refinement of elements, which we mark through function 
    //! SetRefinementFlag
    virtual void Refine(const UInt numLoops = 1)
    { Error( "This fnc is valid only for adaptgrid. Change your config-file",
             __FILE__ ,__LINE__); }

    //! Do uniform refinement of elements, which we mark through function 
    //! 'SetRefinementFlag'
    virtual void RefineUniform()
    { Error( "Not implemented", __FILE__, __LINE__); }

    //! Update nodes for boundary conditions
    //! \param bcs list of boundary nodes
    virtual void UpdateBCs(std::list<UInt> * bcs)
    { Error( "Not implemented", __FILE__, __LINE__); }
  
    //! Return vector of element-neighbors for the element with number noOfElem
    //! \param noOfElem (in) element level
    //! \param color (in) subdomain
    virtual  StdVector<Elem*> *GetNeighboursOfElem(const UInt noOfElem, 
                                                   std::string color)
    { 
      Error(" Not implemented",__FILE__,__LINE__);
      StdVector<Elem*> *dummy = NULL;
      return dummy;
    }
 

    //! Return vector of element-neighbors for the node with number noOfNode
    //! \param noOfNode (in) global number of node
    //! \param neighbours (out) list with neighbors
    virtual void GetNeighboursOfNode(const UInt noOfNode,
                                     StdVector<Elem*> * neighbours)
    { Error(" Not implemented",__FILE__,__LINE__);}
    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous  
 

    //! Returns for a given region identifier the associated name
    std::string RegionIdToName( const RegionIdType regionId );

    //! Returns for a given region name the associated identifier
    RegionIdType RegionNameToId( const std::string & regionName );

    //! Returns for a list of region identifiers the associated names
    void RegionIdToName( StdVector<std::string> & regionNames,
                         const StdVector<RegionIdType> & regionId );
  
    //! Returns for a list of region names  the associated identifiers
    void RegionNameToId( StdVector<RegionIdType> & regionIds,
                         const StdVector<std::string> 
                         & regionName );
  
    //! Returns node numbers of a list of Elements

    //! This method returns the unique node numbers of
    //! a list of given elements. Ther are no duplicate entries.
    //! \param nodeList (out) list of unique node numbers in elemList
    //! \param elemList (in) list of elements
    //! \param onlyLinNodes (in) if true, only the corner nodes are retrieved
    virtual void GetNodesOfElemList( StdVector<UInt> & nodeList,
                                     const StdVector<Elem*> & elemList,
				     Boolean onlyLinNodes = FALSE ) = 0;
  
    //! Resets the integration type of all known elements
    //! \deprecated Does anyone need this function?
    void SetIntTypeAllElems( IntegrationType aIntType );

    //@}
  

    // =======================================================================
    // Methods just needed for StructGrid-Class
    // =======================================================================
    //@{ \name StructGrid  

    //! reads the grid
    virtual void GenGridStruct(const UInt elemx,const UInt elemy, 
			       const UInt elemz) {;};

    //! Transforms the grid
    virtual void TransformGridStruct(UInt & shiftFactor, UInt & nodeShift,
				     UInt & elemgrid, Double &  meshsize, const UInt flag){;};

    //! Returns maximum number of elements in x,y,z-direction
    virtual Integer GetMaxElem(std::string dir)
    { Error("GetMaxElem in Grid-Class not implemented",__FILE__, __LINE__); 
      return 1;};

     //@}

  
  protected:

    //! Get vector containing all region names

    //! Get a vecor which contains all region nodes. The order is in that way,
    //! that one can access directly the elements of the vector by using a
    //! RegionId and get the according entry of the vector.
    virtual void GetAllRegionNames( StdVector<std::string> 
                                    & regionNames ) = 0;

    //! Pointer to mesh input file
    FileType * ptFileType; 
  
    //! List of region identifiers
    StdVector<std::string> regionNames_;
  

  private:
    ///
  };


} // end of namespace
#endif // FILE_GRID
