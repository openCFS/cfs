// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include <list>
#include <set>
#include <def_use_interpolation.hh>

#ifdef USE_INTERPOLATION
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#endif

#include "Domain/elem.hh"
#include "Elements/basefe.hh"
#include "Domain/surfElem.hh"
#include "Domain/ncElem.hh"
#include "Domain/edgeFace.hh"
#include "Domain/entityList.hh"
#include "DataInOut/Scripting/scriptable.hh"
#include "MatVec/vector.hh"


namespace CoupledField
{

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

  class Grid : public Scriptable {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization

    //! Constructor

    //! Standard Constructor
    //! \param aptFileType pointer to FileType for reading initial grid
    Grid();

    //! Destructor
    virtual ~Grid();

    //! Trigger mapping of elements' faces

    //! This method calculates global surface numbers and
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    virtual void MapFaces() = 0;

    //! Trigger mapping of edges

    //! This method calculates global edge and surface numbers and
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    virtual void MapEdges() = 0;

    //@}

    // =======================================================================
    // GENERAL GRID INFORMATION
    // =======================================================================
    //@{ \name General Grid Information

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ NODE INFORMATION +++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    virtual void AddNodes(const UInt numNodes)
    { EXCEPTION( "Not implemented" ); }

    virtual void SetNodeCoordinate(const UInt numNode, const Point & rfPoint)
    { EXCEPTION( "Not implemented" ); }

    virtual void SetNodeCoordinate(const UInt numNode, const Vector<Double> & rfPoint)
    { EXCEPTION( "Not implemented" ); }

    //! Return dimension of mesh

    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    virtual UInt GetDim() = 0;

    //! Returns the maximum node number in the finite element grid.
    virtual UInt GetNumNodes() = 0;

    //! Return if grid uses quadratic elements
    virtual bool IsQuadratic() = 0;

    //! Get coordinates of node with global number inode
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetNodeCoordinate( Point & rfPoint,
                                    const UInt inode,
                                    bool updated = false )
    { EXCEPTION( "Not implemented" ); }

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (out) coordinates of point 3D
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetNodeCoordinate( Vector<Double> & rfPoint,
                                    const UInt inode,
                                    bool updated = false ) = 0;

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

    //! Get coordinates of element nodes

    //! Returns the coordinates of all nodes of one element.
    //! \param coordMat (out) coordinates of the element nodes
    //!                         (spaceDim \f$\times\f$ nrNodes);
    //! \param connect (in) global node numbers of element
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetElemNodesCoord( Matrix<Double> & coordMat,
                                    const StdVector<UInt> & connect,
                                    bool updated = false ) = 0;


    //! Set offset for a single node, called by shape optimization
    virtual void SetNodeOffset( const UInt node, const Point& offset ) = 0;

    //! Get offset for a single node, called by shape optimization
    virtual void GetNodeOffset(const UInt node, Point& offset) = 0;

    //! Set offset for coordinates due to updated Lagrangian formulation
    virtual void SetNodeOffset( const StdVector<UInt>& nodes,
                                const Vector<Double>& offsets ) = 0;

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    virtual bool HasNodalOffset() = 0;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ ELEM INFORMATION +++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    virtual void AddElems(UInt nElems) = 0;

    virtual void SetElemData(UInt ielem,
                             Elem::Elem::FEType type,
                             RegionIdType region,
                             const UInt* connect) = 0;

    virtual void GetElemData(const UInt ielem,
                             Elem::FEType & type,
                             RegionIdType & region,
                             UInt* connect) const
    { EXCEPTION( "Not implemented" ); }


    //! Returns the total number of elements in the grid
    virtual UInt GetNumElems() = 0;

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    virtual UInt GetNumElemOfType( Elem::FEType type ) = 0;

    //! Returns the total number of volume elements in the finite element grid
    virtual UInt GetNumVolElems() = 0;

    //! Returns the total number of surface elements in the grid
    virtual UInt GetNumSurfElems() = 0;

    //! Get element with given element number

    //! Returns a single element with the given element number
    //! \param elemNr element number
    virtual const Elem * GetElem( UInt elemNr ) = 0;


    //! Get node numbers of given element

    //! Returns the node numbers of a  given element.
    //! \param connect (out) contains global node numbers
    //! \param iElem (in) element number
    virtual void GetElemNodes( StdVector<UInt> & connect,
                               const UInt iElem ) = 0;

    //! Returns surface element normal without defined orientation

    //! This method calculates the normal of an surface element. The direction
    //! is unspecified and dependend on the order or nodes of the element.
    //! However, each surface element stores a flag \a normalSign, which
    //! indicates the direction of the resulting normal w.r.t. its first
    //! volume element pointer.
    //! \param n (out) vector containing surface normal
    //! \param surfElem (in) reference to surface element
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void CalcSurfNormal( Vector<Double> & n,
                                 const Elem & surfElem,
                                 bool updated = false ) = 0;

    //! Returns surface element normal with defined orientation

    //! Calculates the surface normal pointing OUT OF the neighbouring
    //! volume element
    //! \param n (out) normal vector
    //! \param surfElem (in) surface element
    //! \param volElem (in) volume element
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                         const Elem & surfElem,
                                         const Elem & volElem,
                                         bool updated = false ) = 0;

    //! Returns node numbers of a list of Elements

    //! This method returns the unique node numbers of
    //! a list of given elements. Ther are no duplicate entries.
    //! \param nodeList (out) list of unique node numbers in elemList
    //! \param elemList (in) list of elements
    //! \param onlyLinNodes (in) if true, only the corner nodes are retrieved
    virtual void GetNodesOfElemList( StdVector<UInt> & nodeList,
                                     const StdVector<Elem*> & elemList,
				     bool onlyLinNodes = false ) = 0;

    //! Returns the global midpoint of an element

    //! This method return the global midpoint of an element
    virtual void GetGlobalElemMidPoint( UInt elemNum, Vector<Double>& coord ) = 0;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ REGION INFORMATION +++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    virtual void AddRegion(const std::string regionName, RegionIdType & rId);
    virtual void AddRegions(const StdVector<std::string> & regionNames,
                            StdVector<RegionIdType> & rIds);

    virtual UInt GetNumRegions();
    virtual UInt GetNumVolRegions();
    virtual UInt GetNumSurfRegions();

    //! Get vector with all (surface and volume) region identifiers

    //! Return a vector with names of all region identifiers in the
    //! current mesh (surface and volume)
    //! \param volRegions (out) vector with region identifiers
    virtual void GetRegionIds( StdVector<RegionIdType> & regions );

    //! Get vector with all volume region identifiers

    //! Return a vector with names of all volume region identifiers in the
    //! current mesh.
    //! \param volRegions (out) vector with volume region identifiers
    virtual void GetVolRegionIds( StdVector<RegionIdType> & volRegions );

    //! Get vector with all surface region identifiers

    //! Return a vector with names of all surface region identifiers in the
    //! current mesh.
    //! \param surfRegions (out) vector with volume region identifiers
    virtual void GetSurfRegionIds( StdVector<RegionIdType>
                                   & surfRegions );

    //! Returns the number of nodes contained in given region
    virtual UInt GetNumNodes( const StdVector<RegionIdType>
                              & regions ) = 0;

    //! Get list of nodes contained in a region

    //! Returns a list of all nodes, which are contained in a
    //! volume- or surface-region.
    //! \param nodeList (out) list with node numbers
    //! \param regionId (in) region identifier
    virtual void GetNodesByRegion( StdVector<UInt> & nodeList,
                                   const RegionIdType regionId ) = 0;

    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of
    virtual UInt GetNumElems( const StdVector<RegionIdType>
                              & regions ) = 0;


    virtual UInt GetMaxNumNodesPerElem()
    { EXCEPTION( "Not implemented" );  return 0;}

    virtual void SetElemRegion(UInt ielem, RegionIdType region)
    { EXCEPTION( "Not implemented" ); }


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


    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume.
    //! 'Volume' here means, that for 2D elements the third dimension is
    //! assumed to be 1m.
    //! \param regionId (in) region identifier
    //! \param isaxi (in) flag indicating axial symmetry
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual Double CalcVolumeOfRegion( const RegionIdType regionId,
                                       bool isaxi = false,
                                       bool updated = false ) = 0;

    /** Calculate the volume spanned by all named nodes of all regions.
     * For most proper convex grids this will be faster than CalcVolumeOfRegion() but one needs to
     * be sure to have named nodes at the exterior points. Will hold for periodic
     * boundary problems */
    Double CalcVolumeSpannedByNamedNodes();

    //! Returns for a given region identifier the associated name
    virtual std::string RegionIdToName( const RegionIdType regionId );

    //! Returns for a list of region identifiers the associated names
    virtual void RegionIdToName( StdVector<std::string> & regionNames,
                                 const StdVector<RegionIdType> & regionId );

    //! Returns for a list of region names  the associated identifiers
    virtual void RegionNameToId( StdVector<RegionIdType> & regionIds,
                                 const StdVector<std::string> & regionName );

    //! Returns for a given region name the associated identifier
    virtual RegionIdType RegionNameToId( const std::string & regionName );



    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++ NAMED NODES INFORMATION ++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    virtual void AddNamedNodes( std::string name, StdVector<UInt> & nodeNums) = 0;

    //! Get list with names of all named nodes

    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    virtual void GetListNodeNames( StdVector<std::string> & nodeNames) = 0;

    virtual void GetElemNumsByName( StdVector<UInt> & elemsNums,
                                    const std::string & elemsName )
    { EXCEPTION( "Not implemented" ); }

    //! Returns the number of nodes in the given nodelist
    virtual UInt GetNumNodes( const std::string & nodesName ) = 0;

    //! Get list of nodes by their name

    //! Returns a list of nodes, which are defined by a name.
    //! These names can specify e.g. boundary nodes or nodes for
    //! saving results.
    //! \param nodeList (out) list with node numbers
    //! \param name (in) name of nodes
    virtual void GetNodesByName( StdVector<UInt> & nodeList,
                                 const std::string & name ) = 0;


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++ NAMED ELEMS INFORMATION ++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    virtual void AddNamedElems( std::string name, StdVector<UInt> & elemNums) = 0;

    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    virtual void GetListElemNames( StdVector<std::string> & elemNames) = 0;

    //! Get list of elements by their name
    virtual void GetElemsByName( StdVector<Elem*> & elems,
                                 const std::string & elemsName ) = 0;





    virtual void FinishInit() = 0;



    //@}

    // =======================================================================
    // ENTITIYLIST ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Entity Access Functions

    //! Get element list of given type with given entities

    //! Returns an EntityList (elements, surface-elements, nodes ) of  the
    //! specified type and with the given name
    //! \param type Type of EntityList to be created
    //! \param name Name of the elements (region) or nodes
    shared_ptr<EntityList> GetEntityList( EntityList::ListType type,
                                          const std::string& name ,
                                          EntityList::DefineType defineType);


    //@}


    //! Get vector containing all region names

    //! Get a vecor which contains all region nodes. The order is in that way,
    //! that one can access directly the elements of the vector by using a
    //! RegionId and get the according entry of the vector.
    virtual void GetRegionNames( StdVector<std::string>
                                 & regionNames ) = 0;

    // =======================================================================
    // NODE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Node Access Functions

    //@}

    // =======================================================================
    // ELEMENT ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Element Access Functions

    //@}


    // =======================================================================
    // SURFACE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Surface Access Functions

    //! Get total number of faces in the grid
    virtual UInt GetNumFaces() = 0;

    //! Return face with given face
    virtual const Face& GetFace( UInt faceNr) = 0;

    //@}

    // =======================================================================
    // ELEMENT EDGE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Edge Access Functions


    //! Get number of edges in the grid
    virtual UInt GetNumEdges() = 0;

    //! Return edge with given number
    virtual const Edge& GetEdge( UInt edgeNr ) = 0;

    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation


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
    { EXCEPTION( "Not implemented" ); }

    //! Do refinement of elements, which we mark through function
    //! SetRefinementFlag
    virtual void Refine(const UInt numLoops = 1)
    { EXCEPTION( "This fnc is valid only for adaptgrid. Change your config-file" ); }

    //! Do uniform refinement of elements, which we mark through function
    //! 'SetRefinementFlag'
    virtual void RefineUniform()
    { EXCEPTION( "Not implemented" ); }

    //! Update nodes for boundary conditions
    //! \param bcs list of boundary nodes
    virtual void UpdateBCs(std::list<UInt> * bcs)
    { EXCEPTION( "Not implemented" ); }

    //! Return vector of element-neighbors for the element with number noOfElem
    //! \param noOfElem (in) element level
    //! \param color (in) subdomain
    virtual  StdVector<Elem*> *GetNeighboursOfElem(const UInt noOfElem,
                                                   std::string color)
    {
      EXCEPTION( "Not implemented" );
      StdVector<Elem*> *dummy = NULL;
      return dummy;
    }


    //! Return vector of element-neighbors for the node with number noOfNode
    //! \param noOfNode (in) global number of node
    //! \param neighbours (out) list with neighbors
    virtual void GetNeighboursOfNode(const UInt noOfNode,
                                     StdVector<Elem*> * neighbours)
    { EXCEPTION( "Not implemented" ); }
    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous

    /** Convenience method for developers, dumps summary information to stdout */
    void Dump();


    /** Probes for the existence of a region name
     * @param name the name to check - case sensitive */
    bool HasRegion(const std::string & regionName)
    {
        return regionNames_.Find(regionName) != -1;
    }

    //! Computes the surface region between two volume regions (2D only)
    virtual void SurfRegionFromVolRegions(
        const std::string& surfRegionName,
        const std::string& region1,
        const std::string& region2);

    //! Computes the surface region around one volume region (outline, 2D only)
    virtual void SurfRegionFromSingleVolRegion(
        const std::string& surfRegionName,
        const std::string& region);
    //@}


    // =======================================================================
    // NONCONFORMING GRID SECTION
    // =======================================================================
    //@{ \name Methods for nonconforming grids

    enum ISecType
    {
      LINE_LINE,
      TRI_TRI,
      TRI_QUAD,
      QUAD_QUAD
    };

    ISecType itype;

    //! NC_SIMON: Main nonmatching grid intersection method
    bool InitNonmatchingInterfaces();

    //! NC_SIMON: check if NC interface is coplanar
    bool IsNcInterfaceCoplanar(const std::string &ncIfaceName);
    bool IsNcInterfaceCoplanar(RegionIdType regionId);

    //! NC_SIMON: intersect two line elements
    bool SideOnSide (SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                     bool coplanar, bool coordUpdate,
                     StdVector<NCElem*>& elemList);

    //! NC_SIMON: intersect two axiparallel quads
    bool RectangleOnRectangle(SurfElem *ifaceElem1, SurfElem *ifaceElem2,
                              bool coplanar, StdVector<NCElem*>& elemList);

    //! intersect two elements of arbitrary type
    bool PolygonOnPolygon(SurfElem *ifElem1, SurfElem *ifElem2,
                          bool coplanar, bool coordUpdate,
                          StdVector<NCElem*> &elemList,
                          Double absTol = 1e-12, Double relTol = 1e-4);

    //! NC_SIMON: add a node to the grid
    //! \param coord (in) coordinates of point
    //! \param inode (out) node number
    virtual void AddNode( const Point & coord, UInt & inode)
    { EXCEPTION( "Not implemented" ); }

    virtual void AddNode( const Vector<Double> & coord, UInt & inode )
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: add multiple nodes to the grid
    //! \param coords (in) coordinates of points
    //! \param inode (out) node numbers
    virtual void AddNodes( const StdVector< Point > & coords,
                           StdVector< UInt > & inodes)
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: Add surface elements
    //! \param regionId (in) elements will be added to region with this id
    //! \param surfelems (in) surface elements to be added
    //! \param elemids (out) element id numbers returned
    virtual void AddSurfaceElems( const RegionIdType regionid,
                                  const StdVector< SurfElem* > & surfelems,
                                  StdVector< UInt > & elemids)
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: Add volume elements
    //! \param regionId (in) elements will be added to region with this id
    //! \param volelems (in) volume elements to be added
    //! \param elemids (out) element id numbers returned
    virtual void AddVolumeElems(  const RegionIdType regionid,
                                  const StdVector< Elem* > & volelems,
                                  StdVector< UInt > & elemids)
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: Add a new Surface region to the grid
    //! \param name (in) name of the new region
    //! \param regionid (out) id of the new region
    virtual void AddSurfaceRegion( const std::string name,
                                   RegionIdType& regionid)
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: Add a new volume region to the grid
    //! \param name (in) name of the new region
    //! \param regionid (out) id of the new region
    virtual void AddVolumeRegion( const std::string name,
                                  RegionIdType& regionid)
    { EXCEPTION( "Not implemented" ); }

    //! NC_SIMON: Remove all elements from the given region
    //! \param regionid (in) id of the region
    virtual void ClearRegion( const RegionIdType regionid)
    { EXCEPTION( "Not implemented" ); }

    //@}


  protected:

    // =======================================================================
    // Method wrappers for function tracing
    // =======================================================================
    //@{ \name Scripting wrapper functions

    //! Trigger function registering
    void RegisterFunctions();

    void Wrap_GetNodeCoordinate();
    void Wrap_GetNodesByName();
    void Wrap_GetNodesByRegion();
    void Wrap_GetListNodeNames();
    void Wrap_GetListElemNames();
    void Wrap_GetRegionNames();
    void Wrap_GetNumNodes();
    void Wrap_GetNumElems();
    void Wrap_GetNumSurfElems();
    void Wrap_GetNumVolElems();
    void Wrap_GetNumNodesOfRegion();
    void Wrap_GetNumElemsOfRegion();
    //@}


    //! List of region identifiers
    StdVector<std::string> regionNames_;

    //! Vector with elements for each volume region
    StdVector<StdVector<Elem*> > volElems_;

    //! Vector with node numbers for each volume region
    StdVector< std::set<UInt> > volElemNodes_;

    //! Vector with region ids for each volume region
    StdVector<RegionIdType> volRegionIds_;

    //! Vector with elements for each surface region
    StdVector<StdVector<Elem*> > surfElems_;

    //! Vector with node numbers for each surface region
    StdVector< std::set<UInt> > surfElemNodes_;

    //! Vector with region ids for each surface region
    StdVector<RegionIdType> surfRegionIds_;

    //! Map from name to type of entity
    std::map<std::string, EntityList::DefineType> nameTypeMap_;


    // =======================================================================
    // Non-matching grid interface calculation
    // =======================================================================

    //@{ \name Non-matching grid helper functions

    //! type of intersection calculation for ncInterfaces
    enum IntersectType { LINE_INTERSECT, RECT_INTERSECT, POLYGON_INTERSECT };

    //! return codes of function CutLines
    enum LineIntersectType {
      INTERSECT_NONE,
      INTERSECT_OUTSIDE,
      INTERSECT_ON_LINE2,
      INTERSECT_CROSS,
      INTERSECT_IN_A,
      INTERSECT_IN_B,
      INTERSECT_IN_C,
      INTERSECT_IN_D,
      INTERSECT_A_EQ_C
    };

    //! struct for keeping all information of an ncInterface
    struct ncInterface {
      std::string name;
      RegionIdType region;
      RegionIdType masterRegion;
      RegionIdType slaveRegion;
      IntersectType intersectAlgo;
      Double tolAbs;
      Double tolRel;
      Double rotationAngle;
      std::string coordSysId;
      bool coplanar;
    };

    //! returns if a list of surface elements are all coplanar
    bool IsSurfacePlanar(const StdVector<SurfElem*>& surfElems);

    //! updates the intersection mesh of a ncInterface
    void UpdateNcIntersection(const ncInterface& ncIf);

    //! Calculates the intersection of two lines [a,b] and [c,d] and stores
    //! the intersection point (if any) in e.
    //! return codes:
    //! INTERSECT_NONE: no intersection
    //! INTERSECT_OUTSIDE: lines [a,inf) and [c,inf) intersect outside [c,d]
    //! INTERSECT_ON_LINE2: lines [a,inf) and [c,d] intersect on [c,d]
    //! INTERSECT_CROSS: lines (a,b) and (c,d) intersect in e
    //! INTERSECT_IN_A: intersection in a
    //! INTERSECT_IN_B: intersection in b
    //! INTERSECT_IN_C: intersection in c
    //! INTERSECT_IN_D: intersection in d
    //! INTERSECT_A_EQ_C: a equals c
    //! \param a (in) starting point of line 1
    //! \param b (in) end point of line 1
    //! \param c (in) starting point of line 2
    //! \param d (in) end point of line 2
    //! \param e (out) intersection point
    LineIntersectType CutLines(const Vector<Double> &a,
        const Vector<Double> &b, const Vector<Double> &c,
        const Vector<Double> &d, Vector<Double> &e) const;

    //! Calculates the intersection between two polygons
    //! \param p1 (in) first polygon to be intersected
    //! \param p2 (in) second polygon to be intersected
    //! \param r (out) resulting polygon
    bool CutPolys(StdVector< Vector<Double> > &p1,
	    StdVector< Vector<Double> > &p2, const bool coPlanarIface,
	    StdVector< Vector<Double> > &r) const;

    //! Returns if a point lies inside a convex polygon. Optionally the
    //! centroid of the polygon may be given in order to save computational
    //! time.
    //! \param p (in) coordinates of the point of interest
    //! \param poly (in) polygon to be tested
    //! \param c (in) pointer to the centroid of the polygon (may be NULL)
    bool PointInsidePoly(const Vector<Double> &p,
	    const StdVector< Vector<Double> > &poly,
	    const Vector<Double> *const c) const;

    //! Computes the centroid of a polygon and returns the radius of the
    //! surrounding circle.
    //! \param p (in) polygon for which the centroid is to be computed.
    //! \param c (out) coordinates of the centroid
    Double PolyCentroid(const StdVector< Vector<Double> > &p,
	    Vector<Double> &c) const;

    //! Computes a triangulation for a polygon.
    //! \param p (in) polygon to be triangulized
    //! \param tri (out) list of triangles
    UInt TriangulatePoly(const StdVector< Vector<Double> > &p,
	    StdVector<NCElem*> &tri);

    //! vector of all ncInterfaces defined
    StdVector<ncInterface> ncInterfaces_;

    //! variable for tolerance of absolute values (used in polygonal
    //! intersection algorithm)
    double polysectAbsTol_;

    //! variable for tolerance of relative values (used in polygonal
    //! intersection algorithm)
    double polysectRelTol_;

    //@}
    
    
    // this is interpolation stuff, but it must be defined,
    // even if you switch off interpolation
  public:
    
    struct ciTolerance {
      Double start;
      Double end;
      Double inc;
      
      ciTolerance(Double s, Double e, Double i) :
        start(s), end(e), inc(i) {};
    };

  private:

#ifdef USE_INTERPOLATION
    typedef CGAL::Box_intersection_d::Box_d<int,2> Box;

    //  typedef CGAL::Box_intersection_d::Box_d<double,2> Box2D;
    //  typedef CGAL::Box_intersection_d::Box_d<double,3> Box3D;
    //  typedef CGAL::Bbox_2                              BBox2D;
    typedef CGAL::Bbox_3                              BBox3D;
    typedef CGAL::Box_intersection_d
                ::Box_with_handle_d<double,3,const UInt*> HandleBox;

    typedef CGAL::Cartesian<double> K;
    typedef K::Point_2 Point2D;


    std::vector<HandleBox> elemBoxes_;

  // callback function object writing results to an output iterator
  template <class OutputIterator, class Grid>
  struct Report2 {
    OutputIterator it;
    Grid& grid;
    Report2(OutputIterator i, Grid& g) : it(i), grid(g) {} // store iterator in object
    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const HandleBox& a, const HandleBox& b) {
      UInt elemNum = *a.handle();
      UInt dim = grid.GetDim();
      UInt localDim;
      UInt numElemNodes;
      Elem::FEType type;
      RegionIdType region;
      std::vector<UInt> connect;
      Matrix<Double> coordMat;
      Matrix<Double> globCoordMat;
      Matrix<Double> localCoords;
      Vector<Double> point;
      StdVector<bool> coordsInside;
      Vector<Double> locCoords;

      connect.resize(64);

      std::cout << "Elem Number " << elemNum << std::endl;

      grid.GetElemData(elemNum, type, region, &connect[0]);
      numElemNodes = Elem::GetNumElemNodes(type);
      coordMat.Resize(dim, numElemNodes);
      globCoordMat.Resize(dim, 1);

      for(UInt i=0; i<numElemNodes; i++)
      {
        grid.GetNodeCoordinate(point, connect[i], true);
        // coordMat auffÃ¼llen!
        for(UInt j=0; j<dim; j++)
        {
          coordMat[j][i] = point[j];

          //          std::cout << "Corner " << (i+1) << " Coord "
          //                    << (j+1) << " " << point[j] << std::endl;
        }
      }

      for(UInt j=0; j<dim; j++)
      {
        globCoordMat[j][0] = b.min_coord(j);

        //        std::cout << "Glob Coord " << (j+1)
        //                  << " Coord " << b.min_coord(j) << std::endl;

      }

      grid.GetElem(elemNum)->ptElem->Global2LocalCoords(localCoords, globCoordMat, coordMat);

      grid.GetElem(elemNum)->ptElem->CoordsInsideElem(localCoords, 0, coordsInside);

      localDim = localCoords.GetNumRows();
      locCoords.Resize(localDim);

      for(UInt j=0; j<localDim; j++)
      {
        //        std::cout << "Local Coord " << (j+1)
        //                  << " Coord " << localCoords[j][0] << std::endl;

        locCoords[j] = localCoords[j][0];
      }

      if(coordsInside[0])
      {
        std::pair<const Elem*, Vector<Double> > pair;
        pair.first = grid.GetElem(elemNum);
        pair.second = locCoords;
        *it++ = pair;
      }

    }
  };
  template <class Iter, class Grid> // helper function to create the function object
  Report2<Iter, Grid> report2(Iter it, Grid& g) { return Report2<Iter, Grid>(it, g); }




  // callback function object writing results to an output iterator
  struct ConsInterpReportFunctor {
      //    const ElemList& destElemList_;
    //    const NodeList& sourceNodeList_;
    Grid* sourceGrid_;
    Grid* destGrid_;
    const std::vector< Vector<Double> >& nodeCoords_;
    Double localEpsilon_;
    std::vector< std::map<UInt, Double> >& consInterpWeights_;
    std::vector<UInt> connect_;
    UInt nodeCounter_;
    UInt numSourceNodes_;
    UInt percentage_;
    UInt oldPercentage_;
    std::map<UInt, UInt> destNodeNumToPosMap_;

    ConsInterpReportFunctor(const ElemList& destElemList,
                            const NodeList& sourceNodeList,
                            const std::vector< Vector<Double> >& nodeCoords,
                            Double localEpsilon,
                            std::vector< std::map<UInt, Double> >& consInterpWeights);

    // We write the id-number of box a to the output iterator assuming
    // that box b (the query box) is not interesting in the result.
    void operator()( const HandleBox& a, const HandleBox& b);
  };
    
    


    // helper function to create the function object
    ConsInterpReportFunctor
    GenConsInterpReportFunctor(const ElemList& destElemList,
                               const NodeList& sourceNodeList,
                               const std::vector< Vector<Double> >& nodeCoords,
                               Double localEpsilon,
                               std::vector< std::map<UInt, Double> >& consInterpWeights)
    {
      return ConsInterpReportFunctor(destElemList,
                                     sourceNodeList,
                                     nodeCoords,
                                     localEpsilon,
                                     consInterpWeights);
    }



    void intersection();

  public:
    const Elem* GetElemAtGlobalCoord(const Vector<double>& globCoord,
                                     Vector<double>& localCoords);

    void ComputeConservativeInterpolationWeights(const ElemList& destElemList,
            const NodeList& sourceNodeList,
            const std::string& coordSysId,
            ciTolerance& globalEpsilon,
            ciTolerance& localEpsilon,
            Double z,
            Double zEpsilon,
            std::vector< std::map<UInt, Double> >& consInterpWeights,
            StdVector<UInt> &unmapped_nodes);

#endif // USE_INTERPOLATION

    ///
  };


} // end of namespace
#endif // FILE_GRID
