#ifndef FILE_SCFE_GRID_2001
#define FILE_SCFE_GRID_2001

#include <def_use_cgal.hh>
#include <def_use_eigen.hh>
#include <def_use_libfbi.hh>

#include <set>
#include <map>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#if defined(USE_CGAL) && defined(USE_LIBFBI)
#error "Either USE_CGAL or USE_LIBFBI can be active, but not both!"
#endif

#ifdef USE_CGAL
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Simple_cartesian.h>
#ifdef USE_EIGEN
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Eigen_svd.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#endif // USE_EIGEN
#endif // USE_CGAL

#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/SurfElem.hh"
#include "Domain/ElemMapping/EdgeFace.hh"
#include "Domain/ElemMapping/EntityLists.hh"


#include "MatVec/Vector.hh"
#include "Forms/IntScheme.hh"

#include "Utils/ThreadLocalStorage.hh"

namespace CoupledField

{

  //! Forward class declaration
  class ResultHandler;
  class CoordSystem;
  class BaseNcInterface;

  //! Class representing geometrical entities (elements, nodes, ...) of a
  //! FE simulation.

  //! This class represents an abstract interface to the geometrical entities
  //! (elements, nodes, ...) in a FE simulation.
  //! It allows access to elements, nodes and coordinates.
  //! The underlying structure is hidden by the according interface classes.
  //!
  //! \note Some General Remarks:
  //! - Node and element numbers are always counted one-based
  //! - The Type \c regionIdType is guaranteed to be of a countable type,
  //! e.g. \c Integer or \c short \c int. Therefore it is always zero-based
  //! and can be directly used as a index for vectors / arrays.
  //! - The region identifiers are used for surface AND volume elements

  class Grid {
      
  public:

    // friend class declaration
    friend class BaseNcInterface;
    friend class MortarInterface;


    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization

    //! Constructor

    //! Standard Constructor
    //! \param param Pointer to <domain> parameter node
    Grid(PtrParamNode param, PtrParamNode infoNode );

    //! Destructor
    virtual ~Grid();

    enum RegionType { NOT_SET, VOLUME_REGION, SURFACE_REGION };


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

    //! Return dimension of mesh

    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    UInt GetDim() const {return dim_;}

    //! Return if grid uses quadratic elements
    virtual bool IsQuadratic() const = 0;

    /** exports the grid to a param node. For command line option --export-grid or for streaming with mesh. */
    virtual void ExportGrid(PtrParamNode out)
    { EXCEPTION( "Not implemented" ); }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ NODE INFORMATION +++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    virtual void AddNodes(const UInt numNodes)
    { EXCEPTION( "Not implemented" ); }

    //! Set coordinate of single node
    
    //! This method sets the coordinate of a single node. 
    //! \param nodenNum (in) node number (1-based)
    //! \param coord (in) coordinate (size: >= dimension of grid)
    //!
    //! \note The dimension of the coordinate must be 3 for 3D and can be 2 or 3
    //!       for 2D. If a 3-dimensional coordinate is provided for a point in 
    //!       2D, we assume the working plane is x-y and check, if the 3rd 
    //!       coordinate component is zero, otherwise an exception is thrown.
    virtual void SetNodeCoordinate( const UInt nodeNum, 
                                    const Vector<Double>& coord ) = 0;

    
    //! Returns the maximum node number in the finite element grid or of a special region
    virtual UInt GetNumNodes(RegionIdType = ALL_REGIONS) const = 0;

    //! Return if grid is axi-symmetric
    bool IsAxi() {return isAxi_; }
    
    //! Set if grid is axisymmetric
    void SetAxi(bool isAxi ) { isAxi_ = isAxi;}

    //! Depth for 2d plane
    Double GetDepth2dPlane() {return depth2dPlane_; }
    
    //! Set if grid is axisymmetric
    void SetDepth2dPlane(Double depth) { depth2dPlane_ = depth;}
    
    
    //! Get coordinates of node (dimension: grid dependent)
    
    //! This method returns the nodal coordinate of the point, where the dimension
    //! of the point corresponds to the grid dimension.
    //! \param rfPoint (out) coordinates of point (size: dim)
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetNodeCoordinate( Vector<Double>& rfPoint,
                                    const UInt inode,
                                    bool updated = false ) const = 0;

    //! Get coordinates of nodeList (dimension: grid dependent)

    //! Basically the same as GetNodeCoordinate but only with
    //! multiple nodes in nodeList
    //! \param nodeCoords (out) coordinates of points
    //! \param inode (in) node numbers
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetNodeCoordinates(StdVector<Vector<Double> >& nodeCoords,
                                    const StdVector<UInt>& nodeList,
                                    const bool updated ) const = 0;
    
    //! Get coordinates of node (dimension: 3D)

    //! This method returns the nodal coordinate of the point always with 3 components.
    //! \param rfPoint (out) coordinates of point (size: 3)
    //! \param inode (in) node number
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetNodeCoordinate3D( Vector<Double>& rfPoint,
                                      const UInt inode,
                                      bool updated = false ) const = 0;

    //! Get elements associated with given node

    //! Returns a list of elements, which have one or more of the given in
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes
    //!                          of nodeList
    //! \param node  (in) node for which neighbouring elements
    //!                      are needed
    virtual void GetElemsNextToNode( StdVector<const Elem*>& elemList,
                                      const UInt& node) = 0;

    //! Get elements associated with given node

    //! Returns a list of elements, which have one or more of the given in
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes
    //!                          of nodeList
    //! \param node  (in) node for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the
    //!                       neihgbouring elements are searched in
    virtual void GetElemsNextToNode( StdVector<const Elem*>& elemList,
                                      const UInt& node,
                                      const StdVector<RegionIdType>& regionIds) = 0;

    //! Get elements associated with given nodes

    //! Returns a list of elements, which have one or more of the given in
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes
    //!                          of nodeList
    //! \param nodeList (in) list of nodes for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the
    //!                       neihgbouring elements are searched in
    virtual void GetElemsNextToNodes( StdVector<const Elem*>& elemList,
                                      const StdVector<UInt>& nodeList,
                                      const StdVector<RegionIdType>& regionIds) = 0;

    //! Get number of elements associated with given nodes

    //! Returns the number of elements, which have one or more of the given in
    //! common. The elements are taken out of a given list of regions.
    //! IMPORTANT: Before using this method, SetNodeNeighbourMap() has to be
    //! called first.
    //! \param num (out) number of elements which have one or more nodes
    //!                          of nodeList
    //! \param node  (in) node for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the
    //!                       neihgbouring elements are searched in
    virtual void GetNumOfElemsNextToNodes( UInt& num,
        const UInt & node,
        const StdVector<RegionIdType>& regionIds) = 0;
        
    virtual void ClearNodeToElemConnectivity() = 0;

    //! Find for every node the number of neighbouring elements

    //! Methods fills the mapNodeToElems_ or mapNodeToElemsNew_ vector with a NodeNeighbourElems-
    //! entry for every volume-region
    //! \param *useNew optional parameter if volume regions are taken into account
    //!                by default, the old version is used
    virtual void SetNodesToElemsMap() = 0;


    //! Get coordinates of element nodes

    //! Returns the coordinates of all nodes of one element.
    //! \param coordMat (out) coordinates of the element nodes
    //!                         (spaceDim \f$\times\f$ nrNodes);
    //! \param connect (in) global node numbers of element
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual void GetElemNodesCoord( Matrix<Double>& coordMat,
                                    const StdVector<UInt>& connect,
                                    bool updated = false ) = 0;

    /** The slower convenient method */
    Matrix<double> GetElemNodesCoord(const Elem* elem, bool updated=true) {
      Matrix<double> ret; // trust return value optimization
      GetElemNodesCoord(ret, elem->connect, updated);
      return ret;
    }

    //! Set offset for coordinates due to updated Lagrangian formulation
    virtual void SetNodeOffset( const StdVector<UInt>& nodes,
                                const Vector<Double>& offsets ) = 0;

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    virtual bool HasNodalOffset() = 0;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ ELEM INFORMATION +++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //! Return the shape representation for a given element. Note that it is
    //! dangerous to call ElemShapeMap::SetElem()!

    //! This method returns the element shape map for a given geometrical
    //! element. In the case of Lagrangian-mapped elements, we keep internally
    //! fixed-instances of element shape maps and just updated their state.
    //!
    //! \param ptElem Pointer to the geometrical element
    //! \param updated Flag for updated Lagrangian geometry
    virtual shared_ptr<ElemShapeMap> GetElemShapeMap( const Elem* ptElem, bool updated = false);
    
    virtual void AddElems(UInt nElems) = 0;

    //! Reserve memory for a number of elements without adding them
    virtual void ReserveElems(UInt nElems) = 0;
    
    virtual void SetElemData(UInt ielem,
                             Elem::FEType type,
                             RegionIdType region,
                             const UInt* connect) = 0;

    virtual void GetElemData(const UInt ielem,
                             Elem::FEType& type,
                             RegionIdType& region,
                             UInt* connect) const
    { EXCEPTION( "Not implemented" ); }


    /* Returns the total number of elements in the grid or of a special region.
     @param reg_id ALL_REGIONS or a volume or surface region id */
    virtual UInt GetNumElems(RegionIdType = ALL_REGIONS) const = 0;

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    virtual UInt GetNumElemOfType( Elem::FEType type ) = 0;
    
    //! Return number of element of a given dimension
    virtual UInt GetNumElemOfDim( UInt dim ) = 0;

    /** Returns the number of volume elements for a special region or all regions
     * @param reg_id if ALL_REGIONS then all volume regions are summed up */
    virtual UInt GetNumVolElems(RegionIdType = ALL_REGIONS) const = 0;

    /** @see GetNumVolElems() */
    virtual UInt GetNumSurfElems(RegionIdType = ALL_REGIONS) const = 0;

    //! Get element with given element number

    //! Returns a single element with the given element number
    //! \param elemNr element number
    virtual const Elem* GetElem( UInt elemNr ) = 0;


    //! Get node numbers of given element

    //! Returns the node numbers of a  given element.
    //! \param connect (out) contains global node numbers
    //! \param iElem (in) element number
    virtual void GetElemNodes(StdVector<UInt>& connect, const UInt iElem) = 0;

    //! Returns node numbers of a list of Elements

    //! This method returns the unique node numbers of
    //! a list of given elements. There are no duplicate entries.
    //! \param nodeList (out) list of unique node numbers in elemList
    //! \param elemList (in) list of elements
    //! \param onlyLinNodes (in) if true, only the corner nodes are retrieved
    virtual void GetNodesOfElemList( StdVector<UInt>& nodeList,
                                     StdVector<const Elem*>& elemList,
				                             bool onlyLinNodes = false ) = 0;

    virtual void GetNodesOfElemList( StdVector<UInt>& nodeList,
                                     StdVector<Elem*>& elemList,
             bool onlyLinNodes = false ) = 0;

    //! Return element at global position and the locally projected coordinate
    
    //! This method returns a FE element, in which the global coordinate
    //! is located, and the related local coordinate. If a list
    //! of regionIds is passed, only elements within these regions are
    //! considered in the search.   
    //! \param globCoord (in) Global coordinate for which the element is requested
    //! \param locCoord (out) Local projection of the global coordinate in the
    //!                      reference element coordinate system.
    //! \param srcEntities(in) (Optional) List or elements, which are considered
    //!                       for the element search. If the set is empty,
    //!                       all (volume) regions are considered. 
    //! \return Element at global coordinate position.
    const Elem* GetElemAtGlobalCoord(const Vector<double>& globCoord, LocPoint& locCoord,
                                     const StdVector<shared_ptr<EntityList> >& srcEntities =StdVector<shared_ptr<EntityList> >(),
                                     bool printWarnings = true, bool updatedGeo = false);
    
    //! Return a list of elements and local coordinate for global coordinates
    void GetElemsAtGlobalCoords( const StdVector<Vector<double> >& globCoords,
                                StdVector< LocPoint >& localCoords,
                                StdVector< const Elem* >& elems,
                                const StdVector<shared_ptr<EntityList> >& srcEntities =
                                StdVector<shared_ptr<EntityList> >(),
                                Double globalTol = 1e-3, 
                                Double localTol = 1e-2,
                                bool printWarnings = true,
                                bool updatedGeo = false);
    
    //! Return for a given node number the element and a local coordinate
    
    //! This method is similar to \ref GetElemAtGlobalCoord, but it takes
    //! a node number rather than a global coordinate. Thus, it will definitely
    //! find an element and the corresponding local coordinate. 
    //! This method is especially useful for retrieving the coordinates of 
    //! mid-side nodes in a higher order FE-simulation, where only the
    //! vertex nodes have natural degree of freedom
    
    //! \note Upon the first call, the internal data structure gets built up,
    //! so that succeeding calls will perform faster.
    const Elem* GetElemAtNode( UInt nodeNum,
                               LocPoint& locCoord,
                               const std::set<RegionIdType>& srcRegions
                               = std::set<RegionIdType>() );


    //! Extract center of element
    virtual void GetElemCentroid(Vector<Double>& center, UInt eNUm, bool isupdated) = 0;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++ REGION INFORMATION +++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /** Adds a region name.
     * Is written to region and regionData but the rest of regionData is only set by FinishInit()
     * @return id the new index within regionData */
    RegionIdType AddRegion(const std::string& regionName, bool regular = false);

    //! Add a new Surface region to the grid
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED
    //! \param name (in) name of the new region
    //! \param regionid (out) id of the new region
    RegionIdType AddSurfaceRegion(const std::string& name) {
      return AddRegion(name, SURFACE_REGION);
    }


    //! Add a new volume region to the grid
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED
    //! \param name (in) name of the new region
    //! \param regionid (out) id of the new region
    RegionIdType AddVolumeRegion(const std::string& name) {
      return AddRegion(name, VOLUME_REGION);
    }


    /** multiple calls of AddRegion().
     * @see AddRegion() */
    void AddRegions(const StdVector<std::string>& names, StdVector<RegionIdType>& ids);

    UInt GetNumRegions() { return regionData.GetSize(); }
    UInt GetNumVolRegions();
    UInt GetNumSurfRegions();

    /** are all elements within this region of the same type and dimension */
    bool IsRegionRegular(RegionIdType reg_id) const { return regionData[reg_id].regular; }

    /** are all regions regular? */
    bool IsRegionRegular(StdVector<RegionIdType>& regions) const;

    /** does the grid consist of regular regions? */
    bool IsGridRegular() const;

    //! Get name of region using a region id
    std::string GetRegionName( RegionIdType id );
    
    //! Get vector containing all region names
    //! Get a vector which contains all region nodes. The order is in that way,
    //! that one can access directly the elements of the vector by using a
    //! RegionId and get the according entry of the vector.
    void GetRegionNames(StdVector<std::string>& out);

    //! Get the region id for the region with name
    RegionIdType GetRegionId(const std::string name);

    //! Get vector with all volume region identifiers

    //! Return a vector with names of all volume region identifiers in the
    //! current mesh.
    //! \param volRegions (out) vector with volume region identifiers
    void GetVolRegionIds( StdVector<RegionIdType>& volRegions );

    //! Get vector with all surface region identifiers

    //! Return a vector with names of all surface region identifiers in the
    //! current mesh.
    //! \param surfRegions (out) vector with volume region identifiers
    void GetSurfRegionIds(StdVector<RegionIdType>& surfRegions );

    //! Returns the number of nodes contained in given region
    UInt GetNumNodes(const StdVector<RegionIdType>& regions) const;

    //! Get list of nodes contained in a region

    //! Returns a list of all nodes, which are contained in a
    //! volume- or surface-region.
    //! \param nodeList (out) list with node numbers
    //! \param regionId (in) region identifier
    virtual void GetNodesByRegion( StdVector<UInt>& nodeList, const RegionIdType regionId ) = 0;

    //! Returns a list of all elements, which are contained in a
    //! volume- or surface-region.
    //! \param elementList (out) list with element numbers
    //! \param regionId (in) region identifier
    virtual void GetElementsByRegion( StdVector<UInt>& elementList, const RegionIdType regionId ) = 0;

    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of
    virtual UInt GetNumElems(const StdVector<RegionIdType>& regions) const = 0;


    virtual UInt GetMaxNumNodesPerElem()
    { EXCEPTION( "Not implemented" );  return 0;}

    virtual void SetElemRegion(UInt ielem, RegionIdType region)
    { EXCEPTION( "Not implemented" ); }

    virtual void GetElemRegion(UInt ielem, RegionIdType& region)
    { EXCEPTION( "Not implemented" ); }

    virtual RegionIdType GetElemRegion(UInt ielem)
    { EXCEPTION( "Not implemented" ); }

    /** Sets the element barycenters for the given region.
     * Checks the RegionData::barycenters and does nothing if already set.
     * @param updated handle updated coordinates?
     * @return the number of actually set barycenters. */
    unsigned int SetElementBarycenters(RegionIdType region, bool updated);

    /** set the element barycenters for all regions */
    unsigned int SetElementBarycenters(bool updated)
    {
      unsigned int total = 0;

      for(unsigned int i = 0; i < regionData.GetSize(); i++)
        total += SetElementBarycenters(regionData[i].id, updated);

      return total;
    };

    /** Determines the neighborhood of elements and store in within the elements.
     * It checks if already called and does nothing if called multiple times. */
    virtual void FindElementNeighorhood()
      { EXCEPTION( "Not implemented" ); }

    //! Get list of elements (surface / volumes)

    //! Returns all elems for a given surface / volume region. If the desired
    //! region consists of surface elements, they are up-casted into Elem*.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetElems(StdVector<Elem*>& elems, const RegionIdType regionId) = 0;


    //! Get list of volume elements

    //! Returns all elements for a given volume region.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetVolElems( StdVector<Elem*>& elems,
                              const RegionIdType regionId ) = 0;

    //! Get list of surface elements

    //! Returns all elements for a given surface region
    //! \param surfElems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    virtual void GetSurfElems( StdVector<SurfElem*>& surfElems,
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
    virtual void GetElemsNextToSurface( StdVector<Elem*>& neighbours,
                                        const StdVector<Elem*>& surfElems,
                                        const StdVector<RegionIdType>& neighRegions ) = 0;

    //! Get the surface element of which \param volElemNum is in ptVolElems[0,1]. Surf elem must be in region reg_id.
    //! If not found, Elem vector \parm surfEl (in) is empty.
    virtual void GetAdjacentSurfElem( const UInt volElemNum, StdVector<Elem *>& surfEl, const RegionIdType reg_id = ALL_REGIONS) = 0;

    //! Get list of volume regions attached to another region. e.g. get all the volume regions of a surface
    //! \param reg_id (in) query parameter
    //! \param volRegIds (out) list
    virtual void GetListOfVolumeRegions( const RegionIdType reg_id, StdVector<RegionIdType>& volRegIds ) = 0;

    /** total volume of a sparse or dense grid */
    virtual Double CalcHullVolume( bool updated = false ) = 0;

    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume.
    //! 'Volume' here means, that for 2D elements the third dimension is
    //! assumed to be 1m.
    //! \param regionId (in) region identifier
    //! \param updated (in) flag indicating if updated geometry should be used
    virtual Double CalcVolumeOfRegion( const RegionIdType regionId,bool updated = false ) = 0;
    
    /** calculates the bounding box of the entire grid. Is slow but cached
     * CalcVolumeSpannedByNamedNodes() in legacy cfs is faster!
     * Is guarded as critical block!
     * @param sys if NULL the the default one is used
     * @param force_3D third dimension is zero for 2D otherwise there is no third dimension in 2D
     * @return num rows == GetDim() or 3 for force_3D */
    Matrix<double>& CalcGridBoundingBox(CoordSystem* sys = NULL, bool force_3D = false);

    /** @see CalcBoundingBoxOfRegion() */
    Matrix<double> CalcRegionsBoundingBox(const StdVector<RegionIdType>& regs, CoordSystem* sys = nullptr);

    //! This method computes the x-y-z boundig box of a given region
    //! \param(in) reId the region to compute the box of
    //! \param(out) minMax Matrix of bounding box values each row
    //!             corresponds to a dimension and has 2 entries
    //!             first column for minimum and second entry for maximum value
    virtual void CalcBoundingBoxOfRegion (const RegionIdType regId,
                                          Matrix<Double>& minMax,
                                          CoordSystem* cSys = nullptr) = 0;

    /** if the grid is regular gives the element discretization, e.g. 20, 20, 1 for 2D.
     * Is reasonable fast but the result is not cached by itself.
     * @return if not regular returns an empty array otherwise always 3 entries of value at least 1*/
    virtual StdVector<unsigned int> CalcRegulardGridDiscretization() { EXCEPTION("not implemented"); }

    //! Returns the volume of a given entitylist (\see Grid::CalcVolumeOfRegion)
    virtual Double CalcVolumeOfEntityList( shared_ptr<EntityList> ent,
                                           bool updated = false ) = 0;


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++ NAMED NODES INFORMATION ++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    virtual void AddNamedNodes(const std::string& name, StdVector<UInt>& nodeNums) = 0;

    //! Get list with names of all named nodes

    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    virtual void GetListNodeNames(StdVector<std::string>& nodeNames) = 0;

    virtual void GetElemNumsByName(StdVector<UInt>& elemNums, const std::string& elemName ) = 0;

    //! Returns the number of nodes in the given nodelist
    virtual UInt GetNumNodes(const std::string& nodesName) const = 0;

    //! Get list of nodes by their name

    //! Returns a list of nodes, which are defined by a name.
    //! These names can specify e.g. boundary nodes or nodes for
    //! saving results.
    //! \param nodeList (out) list with node numbers
    //! \param name (in) name of nodes
    virtual void GetNodesByName(StdVector<UInt>& nodeList, const std::string& name) = 0;


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++ NAMED ELEMS INFORMATION ++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    virtual void AddNamedElems(const std::string& name, StdVector<UInt>& elemNums) = 0;

    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    virtual void GetListElemNames(StdVector<std::string>& elemNames) = 0;

    //! Get list of elements by their name
    virtual void GetElemsByName(StdVector<Elem*>& elems, const std::string& elemsName ) = 0;

    //! Get all elem neighbors for given node id
    virtual const StdVector<Elem*>& GetElemsByNode(UInt node) = 0;

    /** To be called when all regions are added.
     * Sets the internal element and region structures. */
    virtual void FinishInit() = 0;
    
    //! Create result for grid information (local directions etc., Jacobian
    //! determinant)
    virtual 
    void CreateGridInformation( ResultHandler* ptRes,
                                std::map<std::string, CoordSystem*>& coordSysMap ) = 0;
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
                                          const std::string& name );


    //@}

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

    //! Generate Surface elements used for DG calculation schemes.
    //! \param (in) regionList List of volume regions we want to create surfElems for
    //! \param (out) interiorSurfElems Set of surface elements which are inside of the
    //!              volume defined by the regionList. i.e. they have at least one neighbor
    //! \param (out) exteriorSurfElems Surface elements on the exterior boundary
    //!                                no neighbor element
    virtual void GenerateDGSurfaceElemes(std::set<RegionIdType> regionList,
                                         StdVector<shared_ptr<NcSurfElem> >& interiorSurfElems,
                                         StdVector<shared_ptr<NcSurfElem> >& exteriorSurfElems){
      EXCEPTION( "GenerateDGSurfaceElemes: Not implemented here!");
    }
    //@}

    // =======================================================================
    // Automatic Layer Generation and Geometry Computation
    // =======================================================================
  public:
    //! struct that collects StdVectors to store the node geometry 
    //! for a desired node set.
    //! resizes the vectors to a size numNodes when constructed.
    struct NodeGeometry{
      public:
        //! constructor
        NodeGeometry(UInt numNodes) {
          numNodes_ = numNodes;
          nodeIds_ = StdVector<UInt>(numNodes_);
          normalVectors_ = StdVector<Vector<Double>>(numNodes_);
          minPrincipalVectors_ = StdVector<Vector<Double>>(numNodes_);
          maxPrincipalVectors_ = StdVector<Vector<Double>>(numNodes_);
          minPrincipalCurvatures_ = StdVector<Double>(numNodes_);
          maxPrincipalCurvatures_ = StdVector<Double>(numNodes_);
        };

        //! add additional nodes from other NodeGeometry
        void AddNodes(shared_ptr<NodeGeometry> newNodes) {
          UInt numNewNodes = newNodes->numNodes_;
          // check if there are already nodes inserted
          if (numNodes_ > 0) {
            // append new nodes sorted by node Id, assuming that the new nodes and
            // old nodes are already sorted and that every node is contained only once
            if (newNodes->nodeIds_[0] > nodeIds_[numNodes_-1]) {
              nodeIds_.Append(newNodes->nodeIds_);
              normalVectors_.Append(newNodes->normalVectors_);
              minPrincipalVectors_.Append(newNodes->minPrincipalVectors_);
              maxPrincipalVectors_.Append(newNodes->maxPrincipalVectors_);
              minPrincipalCurvatures_.Append(newNodes->minPrincipalCurvatures_);
              maxPrincipalCurvatures_.Append(newNodes->maxPrincipalCurvatures_);
            }
            else if (newNodes->nodeIds_[numNewNodes-1] < nodeIds_[0]) {
              StdVector<UInt> oldNodeIds = nodeIds_;
              StdVector<Vector<Double>> oldNormalVectors = normalVectors_;
              StdVector<Vector<Double>> oldMinPrincipalVectors = minPrincipalVectors_;
              StdVector<Vector<Double>> oldMaxPrincipalVectors = maxPrincipalVectors_;
              StdVector<Double> oldMinPrincipalCurvatures = minPrincipalCurvatures_;
              StdVector<Double> oldMaxPrincipalCurvatures = maxPrincipalCurvatures_;

              nodeIds_ = newNodes->nodeIds_;
              normalVectors_ = newNodes->normalVectors_;
              minPrincipalVectors_ = newNodes->minPrincipalVectors_;
              maxPrincipalVectors_ = newNodes->maxPrincipalVectors_;
              minPrincipalCurvatures_ = newNodes->minPrincipalCurvatures_;
              maxPrincipalCurvatures_ = newNodes->maxPrincipalCurvatures_;

              nodeIds_.Append(oldNodeIds);
              normalVectors_.Append(oldNormalVectors);
              minPrincipalVectors_.Append(oldMinPrincipalVectors);
              maxPrincipalVectors_.Append(oldMaxPrincipalVectors);
              minPrincipalCurvatures_.Append(oldMinPrincipalCurvatures);
              maxPrincipalCurvatures_.Append(oldMaxPrincipalCurvatures);
            } else {
              EXCEPTION("NodeGeometry::AddNodes: Inserted nodes must be unique!")
            }
          } else {
            nodeIds_= newNodes->nodeIds_;
            normalVectors_ = newNodes->normalVectors_;
            minPrincipalVectors_ = newNodes->minPrincipalVectors_;
            maxPrincipalVectors_ = newNodes->maxPrincipalVectors_;
            minPrincipalCurvatures_ = newNodes->minPrincipalCurvatures_;
            maxPrincipalCurvatures_ = newNodes->maxPrincipalCurvatures_;
          }
          numNodes_ += numNewNodes;
        };

        // number of contained nodes
        UInt numNodes_;
        // vectors to store the data
        StdVector<UInt> nodeIds_;
        StdVector<Vector<Double>> normalVectors_;
        StdVector<Vector<Double>> minPrincipalVectors_;
        StdVector<Vector<Double>> maxPrincipalVectors_;
        StdVector<Double> minPrincipalCurvatures_;
        StdVector<Double> maxPrincipalCurvatures_;
    };

    //! Check if autoLayerGeneration parameters are specified for a region and call
    //! CreateExternalLayer if so. Return otherwise.
    virtual void TriggerAutoLayerGeneration() {
      EXCEPTION("Grid::TriggerAutoLayerGeneration not overwritten by child class");
    };

    //! Returns the geometry data for a given surface region if it is already computed, triggers the computation if not.
    //! \param geometry (out) pointer to the struct containing the geometry
    //! \param layerGenNode (in) the layer generation node that specifies the surface and the coefficients for geometry computation
    //! \param surfRegionId (in) the surface region on which the computation is performed
    void GetGeometryOnRegionNodes(shared_ptr<NodeGeometry>& geometry, const PtrParamNode layerGenNode, const RegionIdType& surfRegionId);

    //! checks if there are assigned surface regions to a passed volume region.
    //! Raises exception if the connection has not been set yet. 
    void GetConnectedSurfaceRegions(StdVector<RegionIdType>& connecedSurfRegionIds, const RegionIdType& volumeRegionId);

  protected:
    //! Computes an external grid layer that can be used as a PML region. 
    //! The actual function is implemented in GridCFS
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    virtual void GenerateExternalLayer(const RegionIdType surfRegionId, const PtrParamNode layerGenNode) {
      EXCEPTION("Grid::GenerateExternalLayer not overwritten by child class");
    };

    //! This function triggers the computation of the geometry (normal vectors, principal vectors, 
    //! and principal curvatures) on every node in a surface region.
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    virtual void ComputeGeometryOnSurfaceRegionNodes(const PtrParamNode layerGenNode) {
      EXCEPTION("Grid::ComputeGeometryOnSurfaceRegionNodes not overwritten by child class");
    };

    //! This function triggers the computation of the geometry (normal vectors, principal vectors, 
    //! and principal curvatures) on every node in a surface region.
    //! It exploits the knowledge of a simple geometry type to speed up computation and improve robustness.
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    virtual void ComputeGeometryOfSimpleType(const RegionIdType& surfRegionId, const PtrParamNode layerGenNode) {
      EXCEPTION("Grid::ComputeGeometryOfSimpleType not overwritten by child class");
    };

    //! This function triggers reading the geometry from a file (normal vectors, principal vectors, 
    //! and principal curvatures on every node in a surface region).
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    virtual void ReadGeometryFromFile(const PtrParamNode layerGenNode) {
      EXCEPTION("Grid::ReadGeometryFromFile not overwritten by child class");
    };

    //! This function triggers writing the geometry that is computed by ComputeGeometryOnSurfaceRegionNodes()
    //! to a file.
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    virtual void WriteGeometryToFile(const PtrParamNode layerGenNode) {
      EXCEPTION("Grid::WriteGeometryToFile not overwritten by child class");
    };

    //! map that stores the node geometry for a desired surface region
    //! the RegionIdType is intended to be the key holding the ID of the 
    //! surface region. The NodeGeometry is the struct that holds the data
    std::map<RegionIdType, shared_ptr<NodeGeometry>> geometryRegionMap_;

    //! map that allows to store connected surface and volume regions. 
    //! one volume (key) can hold multiple surface regions (value)
    //! the key is thus the RegionIdType of the volume
    //! the value is the StdVector<RegionIdType>> containing the connected
    //! surfaces.
    std::map<RegionIdType, StdVector<RegionIdType>> volumeSurfaceRegionMap_;


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous
  public:
    
    //! Get type of list denoted by string
    EntityList::DefineType GetEntityType( const std::string& name ) const;
    
    //! Get dimension of region / named elements / nodes with given name
    UInt GetEntityDim( const std::string& name ) const;

    /** Convenience method for developers, dumps summary information to stdout */
    void Dump();

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

    /** Computes for regular grid number of elements in each direction for specified region
        by using maximal and minimal values of barycenters
        @return result vector: [nx ny nz] returns 0 vector, if mesh is not regular */
    StdVector<UInt> GetRegularDiscretization(RegionIdType region);

    // =======================================================================
    // FINITE VOLUME REPRESENTATION SECTION
    // =======================================================================
    //@{ \name Methods for finite volume representations

  public:
  
    //! A struct represting a finite volume (FV) mesh (face based), which usually can not be represented by a CFS grid.
    //! Such a grid might be loaded from Ensight or CCM datasets, but also a FE grid could be interpreted as
    //! Finite volume grid
    struct FiniteVolumeRepresentation {
    
      FiniteVolumeRepresentation();
      
      // If the FV mesh is set
      bool isSet;
      
      // pointer to the FE grid
      Grid* grid;
      
      // A vector storing if elements of the FE grid are volume elements
      std::vector<bool> isVolumeElement;
      
      //! Master elements of the FE Elements which are the results of splitted FV cells into several FE elements.
      //! If the FV cell is a polyhedron splitted into 6 tetrahedron elements [1...6], all of these elements
      //! should have the same master element, e.g. 1. 
      StdVector<UInt> masterElement;
    
      //! Number of faces
      UInt faceCount;
      
      //! stores if a face is a boundary face
      std::vector<bool> isFaceBoundary;
      
      // master elements adjacent to the faces
      StdVector<UInt> faceMasterElement;
      
      // slace elements adjacent to a faces, in case of boundary elements equal to the master element
      StdVector<UInt> faceSlaveElement;
      
      //! index to search for in the facePoint vector to get te points/nodes of the faces
      //! e.g. the nodes for face 10 can be found in the range 
      //! from facePoint[facePointIndex[10]] to facePoint[facePointIndex[11]]
      //! Consequently this vector has to be of length faceCount + 1
      StdVector<UInt> facePointIndex;
      
      //! points/nodes of the faces
      StdVector<UInt> facePoint;
    };
    
    Grid::FiniteVolumeRepresentation& GetFiniteVolumeRepresentation();
    
  private:
    
    FiniteVolumeRepresentation fvr_;


    // =======================================================================
    // NONCONFORMING INTERFACES SECTION
    // =======================================================================
    //@{ \name Methods for nonconforming grids

  public:

    typedef UInt NcInterfaceId;
    
    //! Returns an NcInterface object identified by its ID
    shared_ptr<BaseNcInterface> GetNcInterface(NcInterfaceId ncId) const;
    
    //! Returns an NcInterface object identified by its name
    NcInterfaceId GetNcInterfaceId(const std::string& name) const;
    
    //! Initialize non-conforming interfaces from XML files
    //! Static interfaces are initialized first, followed by moving ones.
    virtual void InitNcInterfacesFromXML();
    
    //! Adds a new NcInterface to the grid and returns its ID
    NcInterfaceId AddNcInterface(shared_ptr<BaseNcInterface> ncIf);
    
    //! Computes if a list of surface elements are all coplanar
    bool IsSurfacePlanar(const StdVector<SurfElem*>& surfElems) const;

    //! Triggers calculation of node offsets for moving interfaces
    //! First, interfaces are triggered to be reset from last to first
    //! and then updated from first to last.
    void MoveNcInterfaces();

    bool HasNCI();


  protected:
    
    //! Add a node to the grid
    //! \param coord (in) coordinates of node
    //! \param inode (out) node number
    virtual void AddNode( const Vector<Double>& coord, UInt& inode )
    { EXCEPTION( "Not implemented" ); }

    //! add multiple nodes to the grid
    //! \param coords (in) coordinates of points
    //! \param inode (out) node numbers
    virtual void AddNodes( const StdVector< Vector<Double> >& coords,
                           StdVector< UInt >& inodes)
    { EXCEPTION( "Not implemented" ); }

    //! Add surface elements
    //! \param regionId (in) elements will be added to region with this id
    //! \param surfelems (in) surface elements to be added
    //! \param elemids (out) element id numbers returned
    virtual void AddSurfaceElems( const RegionIdType regionid,
                                  const StdVector< SurfElem* >& surfelems,
                                  StdVector< UInt >& elemids)
    { EXCEPTION( "Not implemented" ); }
  
    //! Add volume elements
    //! \param regionId (in) elements will be added to region with this id
    //! \param volelems (in) volume elements to be added
    //! \param elemids (out) element id numbers returned
    virtual void AddVolumeElems(  const RegionIdType regionid,
                                  const StdVector< Elem* >& volelems,
                                  StdVector< UInt >& elemids)
    { EXCEPTION( "Not implemented" ); }
    
    //! Remove all elements from the given region
    //! \param regionid (in) id of the region
    virtual void ClearRegion( const RegionIdType regionid)
    { EXCEPTION( "Not implemented" ); }
  
    //! Delete all nodes in a node list. 
    //! Checks if the nodelist is a block at the end, because in this
    //! case it will likely lead to no trouble. 
    //! If the nodes are not at the end, passes an exception because it 
    //! would change connectivity. 
    virtual void DeleteNamedNodes( const std::string& name ) {
      EXCEPTION("Not implemented here.");
    }
  
    //! map for storing ncInterfaces
    StdVector< shared_ptr<BaseNcInterface> > ncInterfaces_;

    //! mapping from ncInterface name to ID
    std::map< std::string, NcInterfaceId > nciNameMap_;


    // =======================================================================
    // Integration Scheme
    // =======================================================================
    
  public:
    
    shared_ptr<IntScheme> GetIntegrationScheme(){
      return integScheme_;
    }

    /** The enum is only for convenience and redundant to regionData.
     * With the addition of the mapping ALL_REGIONS <-> "all".
     * @return note that the size is not correct as ALL_REGIONS is added */
    const Enum<RegionIdType> GetRegion() const { return region_; }

    /** This stores the Region information.
     * The region_id is the index within the vector.
     * @see region */
    struct RegionData
    {
      RegionData();

      RegionIdType id;          // equals the position in the regionData vector
      std::string  name;
      RegionType   type;        // VOLUME_REGION or SURFACE_REGION
      int          type_idx;    // the index within surfRegionIds_ or volRegionIds_
      bool         regular;     // exclusively regular elements of the same size in the region
      bool         homogeneous; // everywhere the same material within the region. In topology optimization this is false
      bool         barycenters; // are the element barycenters calculated?
    };

    /** To be set from FinishInit() */
    StdVector<RegionData> regionData;

  protected:
    
    //! Flag for 2d-plane (default = 1m)
    Double depth2dPlane_;
    
    //! Flag for axi-symmetry
    bool isAxi_;

    /** service for AddSurfaceRegion() and AddVolumeRegion() */
    RegionIdType AddRegion(const std::string& name, RegionType type);

    
    //! Flag for initialization status
    bool isInitialized_;
    
    //! ParamNode of xml file
    PtrParamNode param_;
    
    //! Info node
    PtrParamNode info_;

    /** This is redundant to regionData but convenient for mapping regionId to region_name.
     * Note the the region_id is determined by AddRegion() and not an enum! */
    Enum<RegionIdType> region_;

    //! Vector with elements for each volume region
    StdVector<StdVector<Elem*> > volElems_;

    //! Vector with node numbers for each volume region
    StdVector<UInt> numVolElemNodes_;

    //! Vector with region ids for each volume region
    StdVector<RegionIdType> volRegionIds_;

    //! Vector with elements for each surface region
    StdVector<StdVector<Elem*> > surfElems_;

    //! Vector with node numbers for each surface region
    StdVector<UInt> numSurfElemNodes_;

    //! Vector with region ids for each surface region
    StdVector<RegionIdType> surfRegionIds_;

    //! Dimension of grid
    UInt dim_;

    /** Map from name to type of entity.
    * Therefore the names need to be unique across all element types
    * (e.g. a name must equal for a region and named nodes).
    * Note that the name is also a key for entityDim_ and GridCFS::namedNodeNames_, namedElemNames_, ...
    * The entity structure is split over nameTypeMap_, entityDim_ and GridCFS::namedNodeNames_, GridCFS::namedNodes_ */
    std::map<std::string, EntityList::DefineType> nameTypeMap_;
    
    //! Store dimension for each region / element / node group (name).
    
    //! This map stores for each region / element / node group the dimension.
    //! Node groups have dimension 0 by definition.
    std::map<std::string, UInt> entityDim_;

    // =======================================================================
    // Interation Scheme
    // =======================================================================
    shared_ptr<IntScheme> integScheme_;
    

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

  protected:

    //! Return for a list of entitylists all element numbers
    void GetElemNums( std::set<UInt>& elemNums, 
                      std::set<UInt>& dims,
                      const StdVector<shared_ptr<EntityList> >& entities );
    
    // =======================================================================
    //  ELEMENT / POINT MAPPING
    // =======================================================================

    //@{
    //! Store for mid-side nodes the elements and the local coordinates
    typedef StdVector<std::pair<const Elem*, LocPoint> > NodeElemMatch; 
    boost::unordered_map<UInt, NodeElemMatch> midNodeProjections_;
    //@}
    
    //! Trigger projection of mid-side nodes to element interior
    virtual void MapMidSideNodes() = 0;
    
    //! Structure for mapping global coordinates to a list of
    //! potentially interesting elements
    struct PointElemMatch {

      //! Global coordinate
      Vector<Double> globCoord;

      //! List of potential elements lying in the bounding box
      std::set<const Elem*> matches;

    };

    //! Map a list of global points to element local points
    
    //! This method maps each global coordinate (contained in the PointElemMatch
    //! struct) to a element-local coordinate. It searches all candidate elements
    //! for each point. If the point could not be mapped to any element,
    //! the corresponding Elem-pointer in the elems array is NULL.
    void MapGlobPointsToLoc( const StdVector<PointElemMatch>& matches,
                             StdVector<const Elem*>& elems,
                             StdVector<LocPoint>& lps,
                             Double tol = 1e-2,
                             bool printWarnings = true,
                             bool updatedGeo = false);

  public:
    //! Create a bounding box from a given element. Mapping LIBFBI 
    void CreateBBoxFromElement(const Elem* elem,
                               Double globToler,
                               Double* bbox,
                               double updated=false);

  protected:

#ifdef USE_CGAL
  public:
    //! Define 3-dimensional bounding box
    typedef CGAL::Bbox_3 BBox3D;

    //! Define box handler, which additionally stores an index
    typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,const UInt*> HandleBox;

    //! Define box handler just with an ID
    typedef CGAL::Box_intersection_d::Box_d<double,3> IdBox;

  protected:
    //! Return list of potential elements containing global points

    //! This method returns for every global coordinate a list
    //! of elements, where the point is contained in the bounding box.
    //! This algorithm makes use of the CGAL fast intersection algorithm.
    void MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                   const StdVector<shared_ptr<EntityList> >& srcEntities =
                                   StdVector<shared_ptr<EntityList> >(),
                                   Double tol = 0.0,
                                   bool updatedGeo = false );

    //! Map for each dimension (key) a list containing the "boxes" of elements (value)
    std::map<UInt, std::vector<HandleBox> > elemBoxes_;

    //! \param coords(in) The vector of global coordinates
    //! \param id(in) An identifier for this specific coordinate (e.g. index in a vector)
    //! \param tol(in) Tolerance in meters which determines the size of the bounding box
    HandleBox CreateBoxFromCoord( const Vector<double>& coords,
                                  UInt *id,
                                  Double tol = 0.0 );

#ifdef USE_EIGEN
  public:
    // typedefs to use CGAL for geometry computation (principal directions, curvatures, normal vectors)
    // implemented for use in the curvilinear PML formulation
    //! data type used in the data kernel
    typedef Double                                     DFT;
    //! cartesian data kernel to represent the geometric objects
    typedef CGAL::Simple_cartesian<DFT>                DataKernel;
    //! cartesian local data kernel to represent the geometric objects
    typedef CGAL::Simple_cartesian<DFT>                LocalKernel;
    //! define algebra type to solve linear system
    typedef CGAL::Eigen_svd                            SVD;
    //! representation of a point in 3D Euclidean space
    typedef DataKernel::Point_3                        DPoint;
    //! representation of a vector in 3D Euclidean space
    typedef DataKernel::Vector_3                       DVector;
    //! class to estimate local differential quantities at a given point
    typedef CGAL::Monge_via_jet_fitting<DataKernel, LocalKernel, SVD>    MongeViaJetFitting;
    //! class to store the Monge coordinate system
    typedef MongeViaJetFitting::Monge_form             MongeForm;
  protected:
    //! Function that instantiates a MongeForm, i.e., a surface description of the form z = F(x,y)
    //! The MongeForm can then be used to compute surface parameters, e.g. normal vectors, 
    //! principal directions, or principal curvatures
    //! \param mongeForm (out) created monge form
    //! \param degreePolynomFitting (in) degree of the fitted polynomial to approximate the surface
    //! \param degreeMongeCoeffs (in) degree of the monge coefficients
    //! \param nodeCoordinates (in) coordinates of the vertex (position 0 in vector) and surrounding nodes in DPoint format.
    void SetUpMongeForm(MongeForm& mongeForm, const UInt& degreePolynomFitting, const UInt& degreeMongeCoeffs, 
                        const std::vector<DPoint>& nodeCoordinates);

    //! function that converts a StdVector of Vectors into the Point_3 format that is used by CGAL
    //! !!currently only tested in 3D!!
    //! \param pointsAsPoint_3Format (out) representation as Point_3 
    //! \param pointsAsCfsVectors (in) representation as CFS Vector
    void ConvertVectorToPoint_3Format(std::vector<DPoint>& pointsAsPoint_3Format, StdVector<Vector<Double>>& pointsAsCfsVectors);

    //! function that converts the Point_3 format used by CGAL into a StdVector of Vectors
    //! always outputs 3 dimensional Vecors
    //! !! currently not in use and untested !!
    //! \param pointsAsCfsVectors (out) representation as Point_3 
    //! \param pointsAsPoint_3Format (in) representation as CFS Vector
    void ConvertVectorFromPoint_3Format(StdVector<Vector<Double>>& pointsAsCfsVectors, std::vector<DPoint>& pointsAsPoint_3Format);
#endif // USE_EIGEN
#elif USE_LIBFBI // USE_CGAL
    void MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                   const StdVector<shared_ptr<EntityList> >& srcEntities =
                                   StdVector<shared_ptr<EntityList> >(),
                                   Double tol = 1e-3,
                                   bool updatedGeo = false );
#else // USE_CGAL
    //! Return list of potential elements containing global points (slow version)

    //! This method returns for every global coordinate a list
    //! of elements, where the point is contained in the bounding box.
    //! This method uses the own (potentially slow) algorithm, to determine
    //! the bounding boxed
    void MapPointsToBoundingBoxes( StdVector<PointElemMatch>& matches,
                                   const StdVector<shared_ptr<EntityList> >& srcEntities =
                                   StdVector<shared_ptr<EntityList> >(),
                                   Double tol = 1e-3,
                                   bool updatedGeo = false );
    
    //! Define type for bounding boxes
    typedef std::pair<boost::array<Double,6>, UInt> BoxType;
    
    //! Define for each dimension type (key) bounding boxes (value)
    std::map<UInt, StdVector<BoxType> > elemBoxes_;
#endif // USE_CGAL
  private:
    /** Here the result from CalcGridBoundingBox() is stored and reused. */
    Matrix<double> grid_bounding_box_;
  };
} // end of namespace
#endif // FILE_GRID
