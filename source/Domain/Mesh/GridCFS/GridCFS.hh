#ifndef FILE_SCFE_GRID_CFS_2001
#define FILE_SCFE_GRID_CFS_2001


#include "Domain/Mesh/Grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include <unordered_map>

namespace CoupledField
{

  // Forward class declarations
  struct Elem;
  struct SurfElem;
  // either for nanoflann (kd-tree) or dummy
  struct NodeGridKDTree; 
  struct ElemGridKDTree; 

  //! Implementation of a simple, one level grid.

  //! This class implements the base class Grid. Is is a simple
  //! grid class, which is able to handle one level of geometry without
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
    GridCFS(UInt dim, PtrParamNode param, PtrParamNode infoNode, const std::string &id = "default");
  
    //! Destructor
    ~GridCFS();

    //! Trigger mapping of elements' faces

    //! This method calculates global surface numbers and 
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    void MapFaces() override;
    
    //! Trigger mapping of edges

    //! This method calculates global edge and surface numbers and 
    //! makes them available in the element definitions, so they can
    //! be used for higher order elements or edge functions.
    void MapEdges() override;

    //! Triggers CheckPatternRegion() to check if there is a region 
    //! pattern set in the regionList. Assigns the pattern to respective elements.
    //! Determines the maximum nodes occuring in any element of the grid and
    //! stores it into maxNumElemNodes_.
    //! Assigns a dimension to each region if not set already.
    //! Adds all surface elements to surfElems_, adds all volume elements 
    //! volElems_.
    //! Adds all volRegionIds_. Sets the numVolElemNodes_ for each region.
    //! Discovers which element is a surface element and sets it in orderedElems_.
    //! Sets the surfRegionIds_. Sets the numSurfElemNodes_.
    //! Checks if every region contains at least one volume or surface element.
    //! Calls CheckForRegularRegion() for each region.
    //! Calls CorrectElementConnectivities().
    //! Calls makeNameNodesFromLines().
    //! Calls CreateUserDefinedNodesElems().
    //! Calls CalcRegulardGridDiscretization() and assigns to the mathparser if the 
    //! grid is regular.
    void FinishInit() override;
    
    //! Create result for grid information (local directions etc., Jacobian
    //! determinant)
    void CreateGridInformation( ResultHandler* ptRes,
                                std::map<std::string, CoordSystem*>& coordSysMap )  override;

    //@}

    /**  @see Grid::ExportGrid() */
    void ExportGrid(PtrParamNode out) override;


    // =======================================================================
    // GENERAL GRID INFORMATION
    // =======================================================================
    //@{ \name General Grid Information

    //! Return if grid uses quadratic elements
    bool IsQuadratic() const override { return isQuadratic_; }

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    UInt GetNumElemOfType( Elem::FEType type ) override;
    
    //! Return number of element of a given dimension
    UInt GetNumElemOfDim( UInt dim ) override;

    //! Return dimension of mesh
    //! Returns the geometrical dimension of the mesh. Currently only
    //! two- and three-dimensional meshes are supported.
    using Grid::GetDim;

    //! Add a number of nodes to the grid;
    //! Resizes coords_ and increases numNodes_, but does not assign
    //! coordinates or Ids 
    void AddNodes(const UInt numNodes) override;

    /** @see Grid::GetNumNodes() */
    UInt GetNumNodes(RegionIdType reg_id = ALL_REGIONS) const  override;

    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName ) const override;

    //! increases the counter 'numElems_' by 'nElems';
    //! creates 'nElems' new 'Elem()' objects and adds them to
    //! orderedElems_.
    //! \param nElems (in) number of elements to add
    void AddElems(UInt nElems) override;
      
    //! Reserve memory for a number of elements without adding them
    void ReserveElems(UInt nElems) override;

    //! assign data to the elements with id 'ielem';
    //! \param ielem is the id of the element
    //! \param type is the element type (e.g. Elem::ET_HEXA8)
    //! \param region is the region ID where the element gets assigned to
    //! \param connect is a pointer to the connectivity list
    Elem* SetElemData(UInt ielem,
                             Elem::FEType type,
                             RegionIdType region,
                             const UInt* connect) override;

    void GetElemData(const UInt ielem,
                             Elem::FEType & type,
                             RegionIdType & region,
                             UInt* connect) const override;

    /** @see Grid::GetNumElems() */
    UInt GetNumElems(RegionIdType = ALL_REGIONS) const override;
    
    //! Return maximum number of elements 

    /** @see Grid::GetNumVolElems() */
    UInt GetNumVolElems(RegionIdType = ALL_REGIONS) const override;

    /** @see Grid::GetNumVolElems() */
    UInt GetNumSurfElems(RegionIdType = ALL_REGIONS) const override;
  
    //! Returns number of element contained in one region

    //! Returns the number of element, which belong to a list of given
    //! regions.
    //! \param regions (in) contains the regionIds of 
    UInt GetNumElems(const StdVector<RegionIdType> & regions) const override;
  
    //! Get list with names of all named nodes
    
    //! Get a list with names of all named nodes in the grid
    //! \param nodeNames list of names of nodes
    void GetListNodeNames( StdVector<std::string> & nodeNames) override;

    //! Get list with names of all named elements

    //! Get a list with names of all named elements in the grid
    //! \param elemNames list of names of elements
    void GetListElemNames( StdVector<std::string> & elemNames)  override;

    /** implements the 'elemList' and 'nodeList' handling within 'domain' in the .xml.
     * The user can specify nodes or elements (barycenter) coordinates and named nodes and elements are created.
     * 'coord' is handled by CreateUserDefinedByCoord(), 'bounds' and 'test' are handled by CreateUserDefinedTraverse() 
     * and 'expression' and 'list' are handled by CreateUserDefinedBySearch(). 
     * 'bounds' is much faster than 'list' but only since 03.2026, so almost all 'list' 
     * shall be 'bounds'.
     * @param param with param_->Get("domain") or from RegularMesh. Potential elemList and nodeList is processed. */
    void CreateUserDefinedNodesElems(const PtrParamNode& param);

    /** helper for CreateUserDefinedNodesElems() - handle single node/element defined by coordinate 
     * Only fills the found_entities vector, is in principle const. */
    void SelectByCoord(const std::string& name, bool isNode, const PtrParamNode& pn, StdVector<unsigned int>& found_entities);

    /** 'bounds' and 'test' in parallel: Traverse all elements and check if to add to nodes.
     * As the nodes are traversed in the original order, the entities are also created by sorted index 
     */
    void SelectByTraversal(const std::string& name, bool isNode, const PtrParamNode& pn, StdVector<unsigned int>& found_entities);

    /** Create virtual nodes/elements and searches the nearest in the mesh. 
     *  Implements 'expression' and 'list' */
    void SelectBySearch(const std::string& name, bool isNode, const PtrParamNode& pn, StdVector<unsigned int>& found_entities);

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
                         const std::string & name) override;

    //! Get list of nodes contained in a region

    //! Returns a list of all nodes, which are contained in a volume- or surface-region.
    //! Is quite fast
    //! \param nodeList (out) list with 1-based node numbers, shall be sorted
    //! \param regionId (in) region identifier
    void GetNodesByRegion( StdVector<UInt> & nodeList,
                           const RegionIdType regionId) override;

    //! Returns a list of all elements, which are contained in a
    //! volume- or surface-region.
    //! \param elementList (out) list with element numbers
    //! \param regionId (in) region identifier
    void GetElementsByRegion( StdVector<UInt>& elementList, const RegionIdType regionId ) override;
    
    //! \see Grid::GetNodeCoordinate
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode,
                            bool updated) const override;
    
    //! \see Grid::GetNodeCoordinates
    void GetNodeCoordinates(StdVector<Vector<Double>>& nodeCoords,
                            const StdVector<UInt> & nodeList,
                            const bool updated ) const override;

    /** simply gives access to the full 0-based node coordinates. Note that node numbers are 1-based */
    const StdVector<Vector<double> >& GetNodeCoordinates() const { return coords_;}

    //! \see Grid::GetNodeCoordinate3D
    void GetNodeCoordinate3D( Vector<Double> & rfPoint,
                              const UInt inode,
                              bool updated) const override;
    
    //! \see Grid::SetNodeCoordinate
    void SetNodeCoordinate( const UInt nodeNum, 
                            const Vector<Double> & coord) override;
    //@}

    // =======================================================================
    // ELEMENT ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Element Access Functions

    //! Get element with given element number
    
    //! Returns a single element with the given element number
    //! \param elemNr element number
    const Elem* GetElem(UInt elemNr) override
    {
      assert(orderedElems_[elemNr-1] != NULL);
      assert(orderedElems_[elemNr-1]->elemNum == elemNr);
      return orderedElems_[elemNr-1];
    }

    UInt GetMaxNumNodesPerElem() override
    {
      return maxNumElemNodes_;
    }
      
    void SetElemRegion(UInt ielem, RegionIdType region) override;

    /** @see Grid::FindElementNeighorhood() */
    void FindElementNeighorhood(bool only_volume = true) override;

    void GetElemRegion(UInt ielem, RegionIdType & region) override;

    RegionIdType GetElemRegion(UInt ielem) override;

    //! Get list of elements (surface / volumes)
    
    //! Returns all elems for a given surface / volume region. If the desired 
    //! region consists of surface elements, they are up-casted into Elem*.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetElems( StdVector<Elem*> & elems, const RegionIdType regionId ) override;

    
    //! Get list of volume elements

    //! Returns all elements for a given volume region.
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetVolElems( StdVector<Elem*> & elems, 
                      const RegionIdType regionId ) override;

    //! Get list of surface elements
  
    //! Returns all elements for a given surface region 
    //! \param elems (out) vector with elements for given regionId
    //! \param regionId (in) region identifier
    void GetSurfElems( StdVector<SurfElem*> & elems, 
                       const RegionIdType regionId ) override;

    //! Get list
    void GetElemsByName( StdVector<Elem*> & elems,
                         const std::string & elemsName ) override;

    //! Returns the element numbers of a region or element list.
    //! \param elemNums (out) vector with element numbers of given region
    //! \param elemName (in) name identifying a region or element list
    void GetElemNumsByName( StdVector<UInt> & elemNums,
                               const std::string & elemName ) override;
    
    //! Get node numbers of given element
    //! Returns the node numbers of a  given element.
    //! \param connect (out) contains global node numbers
    //! \param iElem (in) element number
    void GetElemNodes( StdVector<UInt> & connect, 
                       const UInt iElem ) override;

    //! Returns element neighbors of given node
    //! \param node number of interest
    inline const StdVector<Elem*>& GetElemsByNode(UInt node) override
    {
      if (!mappedNodeToElems_)
        SetNodesToElemsMap(); // is thread-safe
      return mapNodeToElems_[node];
    }

    void AddNamedNodes(const std::string& name, StdVector<UInt> & nodeNums) override;
    
    void AddNamedElems(const std::string& name, StdVector<UInt> & elemNums) override;


    //! Get coordinates of element nodes

    //! Returns the coordinates of all nodes of one element.
    //! \param coordMat (out) coordinates of the element nodes 
    //!                         (spaceDim \f$\times\f$ nrNodes);
    //! \param connect (in) global node numbers of element
    //! \param updated (in) flag indicating if updated geometry should be used
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect,
                            bool updated ) override;

    //! Get elements associated with given node

    //! Returns a list of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes of 
    //!                          nodeList
    //! \param node (in) node for which neighbouring elements
    //!                      are needed
   void GetElemsNextToNode( StdVector<const Elem*> & elemList, const UInt & node) override;
    
    //! Get elements associated with given node

    //! Returns a list of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes of 
    //!                          nodeList
    //! \param node (in) node for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the 
    //!                          neihgbouring elements are searched in
    void GetElemsNextToNode( StdVector<const Elem*> & elemList,
                               const UInt & node,
                               const StdVector<RegionIdType>& regionIds) override;
    
    //! Get elements associated with given nodes

    //! Returns a list of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! \param elemList (out) elements which have one or more nodes of 
    //!                          nodeList
    //! \param nodeList (in) list of nodes for which neighbouring elements 
    //!                         are needed
    //! \param regionIds (in) identifiers for the regions, where the 
    //!                          neihgbouring elements are searched in
    void GetElemsNextToNodes( StdVector<const Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> & regionIds) override;

    //! Get number of elements associated with given nodes

    //! Returns the number of elements, which have one or more of the given
    //! common. The elements are taken out of a given list of regions.
    //! IMPORTANT: Before using this method, SetNodesToElemsMap() has to be
    //! called first.
    //! \param num (out) number of elements which have one or more nodes
    //!                          of nodeList
    //! \param node (in) node for which neighbouring elements
    //!                      are needed
    //! \param regionIds (in) identifiers for the regions, where the
    //!                       neihgbouring elements are searched in
    void GetNumOfElemsNextToNodes( UInt & num,
        const UInt & node,
        const StdVector<RegionIdType>& regionIds) override;

    void ClearNodeToElemConnectivity() override;

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
                                &neighRegions ) override;

    //! Get list of volume regions attached to another region. e.g. get all the volume regions of a surface
    //! \param reg_id (in) query parameter
    //! \param volRegIds (out) list
    void GetListOfVolumeRegions( const RegionIdType reg_id, StdVector<RegionIdType> &volRegIds ) override;

    //! Get the surface element of which \param volElemNum is in ptVolElems[0,1]. Surf elem must be in region reg_id.
    //! If not found, Elem vector \parm surfEl (in) is empty.
    void GetAdjacentSurfElem( const UInt volElemNum, StdVector<Elem *> & surfEl, const RegionIdType reg_id = ALL_REGIONS) override;
    
    //! Extract center of element
    void GetElemCentroid(Vector<Double>& center, UInt eNUm, bool isupdated) override;

    //@}

    unsigned int GetNumFaces() const override { return faces_.GetSize(); };

    const Face& GetFace( UInt faceNr) const override {
      assert(faces_.GetSize() > 0);
      return faces_[faceNr-1];
    }
        
    unsigned int GetNumEdges() const override { return edges_.GetSize(); };
    
    const Edge& GetEdge(unsigned int edgeNr ) const override {
      assert(edges_.GetSize() > 0);
      return edges_[edgeNr-1];
   }

    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    //! Returns the volume of a given region

    //! This method returns the volume of a given region by iterating over
    //! all elements (volume / surface) and summing up their volume. 
    //! 'Volume' here means, that for 2D elements the third dimension is 
    //! assumed to be 1m.
    //! \param regionId (in) region identifier 
    //! \param updated (in) flag indicating if updated geometry should be used
    double CalcVolumeOfRegion( const RegionIdType regionId, bool updated = false) override;


    //! \copydoc Grid::CalcVolumeOfEntityList
    Double CalcVolumeOfEntityList( shared_ptr<EntityList> ent,
                                   bool updated = false ) override;

    /** Total volume of a sparse grid. Works only for parallelograms. */
    Double CalcHullVolume(bool updated = false) override;

    //! @copydoc Grid::CalcBoundingBoxOfRegion
    void CalcBoundingBoxOfRegion (const RegionIdType regId,
                                  Matrix<Double> & minMax,
                                  CoordSystem* cSys = nullptr) override;

    /** @see Grid::CalcRegulardGridDiscretization() */
    StdVector<unsigned int> CalcRegulardGridDiscretization() override;


    //@}


    // =======================================================================
    // MISCELLANEOUS
    // =======================================================================
    //@{ \name Miscellaneous  
 
    //! Returns node numbers of a list of Elements

    //! This method returns the unique node numbers of
    //! a list of given elements. There are no duplicate entries.
    //! \param nodeList (out) list of unique node numbers in elemList
    //! \param elemList (in) list of elements
    //! \param onlyLinNodes (in) if true, only the corner nodes are retrieved
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             StdVector<const Elem*> & elemList,
			     bool onlyLinNodes = false) override;
    
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             StdVector<Elem*> & elemList,
           bool onlyLinNodes = false) override;


    //! Set offset for coordinates due to updated Lagrangian formulation
    void SetNodeOffset( const StdVector<UInt>& nodes, 
                        const Vector<Double>& offsets ) override;

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    bool HasNodalOffset() override;
    //@}

    //! add node to the grid
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED 
    void AddNode( const Vector<Double> & coord, UInt & inode ) override;
    
    //! add multiple nodes to the grid
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED 
    //! \param coords (in) coordinates of points
    //! \param inode (out) node numbers
    void AddNodes( const StdVector< Vector<Double> > & coords,
                           StdVector< UInt > & inodes) override;

    //! Add surface elements
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED 
    //! \param regionId (in) elements will be added to region with this id
    //! \param surfelems (in) surface elements to be added
    //! \param elemids (out) element id numbers returned
    void AddSurfaceElems( const RegionIdType regionid,
                                  const StdVector< SurfElem* > & surfelems,
                                  StdVector< UInt > & elemids) override;

    //! Add volume elements 
    //! USAGE OF THIS FUNCTION CAN BE DANGEROUS NOT
    //! ALL NECCESARY FEATURES MAY BE IMPLEMENTED 
    //! \param regionId (in) elements will be added to region with this id
    //! \param volelems (in) volume elements to be added
    //! \param elemids (out) element id numbers returned
    void AddVolumeElems(  const RegionIdType regionid,
                                  const StdVector< Elem* > & volelems,
                                  StdVector< UInt > & elemids) override;

    //! Remove all elements from the given region
    //! \param regionid (in) id of the region
    void ClearRegion( const RegionIdType regionid ) override;

    //! DIRTY HACK WARNING: Deletes the nodes in the list regardless of
    //! whether they are still connected to elements or not!
    //! Used with rotating nonconforming interfaces for deleting the
    //! new nodes of intersection elements after each time step.
    //! \param name (in) name of the named nodes list
    void DeleteNamedNodes( const std::string &name ) override;
    
    /** gives the element with the lowest elemNum for a region.
     * Slow implementation with linear search */
    Elem* SearchFistRegionElement(RegionIdType reg) const;


  private:

    /** checks the domain in the xml file for a pattern region.
     * @param replace resized to elements+1 by element number. True if the
     * pattern region has an entry there. Untouched if no pattern region
     * @return in case of an pattern region it is added to regionData and the
     * id is returned, otherwise NO_REGION_ID */
    RegionIdType CheckPatternRegion(StdVector<bool>& replace);

    /** Helper for FinishInit(). Determines if there are only regular elements.
     * Can be expensive!
     * @return true means that the region is regular */
    bool CheckForRegularRegion(RegionIdType reg);

    //! Find for every node the neighbouring elements

    //! Methods fills the mapNodeToElems_ or mapNodeToElemsNew_ vector with a NodeNeighbourElems-
    //! entry for every volume-region
    void SetNodesToElemsMap() override;


    inline double CalcVolumeOfAllRegions(bool updated=false) {
      // Volume of all regions
      Double s = 0.0;
      for( UInt i = 0; i < volRegionIds_.GetSize(); i++ )
        s += CalcVolumeOfRegion(volRegionIds_[i], updated);
      return s;
    }

    //! helper struct for storing the number of neighbour-elements for every node
    struct NodeNeighbourElems{
      std::unordered_map<UInt, StdVector<Elem*> > nodeNeighElems;
      RegionIdType regID;
    };


    // =======================================================================
    // Automatic Layer Generation and Geometry Computation
    // =======================================================================
  public:
    //! Check if autoLayerGeneration parameters are specified for a region and call
    //! CreateExternalLayer if so. Return otherwise.
    void TriggerAutoLayerGeneration() override;

    //! use parent implementation
    using Grid::GetGeometryOnRegionNodes;

  private:
    //! Computes an external grid layer that is used as a PML region. 
    //! Assigns the new volume region to the grid. 
    //! Additionally, assigns one surface region for each iso-surface layer within the new volume region.
    //! \param surfaceRegion (in) pointer to the surfaceRegion where the layer should be built upon
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    void GenerateExternalLayer(const RegionIdType surfRegionId, const PtrParamNode layerGenNode) override;

    //! This function triggers the computation of the geometry (normal vectors, principal vectors, 
    //! and principal curvatures) on every node in a surface region.
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    void ComputeGeometryOnSurfaceRegionNodes(const PtrParamNode layerGenNode) override;

    //! This function triggers the computation of the geometry (normal vectors, principal vectors, 
    //! and principal curvatures) on every node in a surface region.
    //! It exploits the knowledge of a simple geometry type to speed up computation and improve robustness.
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    void ComputeGeometryOfSimpleType(const RegionIdType& surfRegionId, const PtrParamNode layerGenNode) override;

    //! This function triggers reading the geometry from a file (normal vectors, principal vectors, 
    //! and principal curvatures on every node in a surface region).
    //! Stores the computed data into the geometryRegionMap_
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    void ReadGeometryFromFile(const PtrParamNode layerGenNode) override;

    //! This function triggers writing the geometry that is computed by ComputeGeometryOnSurfaceRegionNodes()
    //! to a file.
    //! \param layerGenNode (in) the param node to the respective 'newRegion' of the layerGenerationList
    void WriteGeometryToFile(const PtrParamNode layerGenNode) override;


    // =======================================================================
    // Helper Methods
    // =======================================================================
  private:
    //@{ \name Helper Methods

    //! Creates the surface elements

    //! This method creates the surface elements, by assigning each surface 
    //! element one or two volume neighbours. Also the flag for indicating
    //! the direction of the surface normal is calculated.
    //! \param elems (input) set containing surface elements which are not
    //!                      yet converted to \a SurfElem*
    //! \param mappedElems (output) set containing mapped surface elements 
    void CreateSurfaceElements( StdVector<Elem*> & elems,
                                StdVector<SurfElem*>& mappedElems );

    //! This method creates surface elements for each "volume" element contained
    //! in the given regions. We assume, that a surface element is always one
    //! of lower space dimension than its vloume element. the new elements get a
    //! element number and are not! contained in the grid
    //! \param (in) regionList the volume regions we want to create surfaces for
    //! \param (out) interiorSurfElems list of surface elements which have at least one neighbor
    //! \param (out) exteriorSurfElems list of elements on the boundary of the volume given
    //!              by the regionList. i.e. they have no surfelem neighbors
    void GenerateDGSurfaceElemes(std::set<RegionIdType> regionList,
                                         StdVector<shared_ptr<NcSurfElem> > & interiorSurfElems,
                                         StdVector<shared_ptr<NcSurfElem> > & exteriorSurfElems) override;

    //! Helper method for GenerateDGSurfaceElemes. We assume here that the NcSurfElem
    //! list hold only valid surface elements. i.e. the normal and everything is already computed
    //! there are two modes. The first one for the conforming case checks neighbor information
    //! according to egde or face numbers, the second one tries to gues the correct neighbors by
    //! performing tests related to the surface normal and intersection areas
    //! \param (in) surfElemList list of elements to be tested
    //! \param (out) interiorSurfElems list of surface elements which have at least one neighbor
    //! \param (out) exteriorSurfElems list of elements on the boundary of the volume given
    //!              by the regionList. i.e. they have no surfelem neighbors
    void ComputeSurfElemNeighbors(StdVector<shared_ptr<NcSurfElem> > surfElemList,
                                          StdVector<shared_ptr<NcSurfElem> > & interiorSurfElems,
                                          StdVector<shared_ptr<NcSurfElem> > & exteriorSurfElems,
                                          bool conforming);

    //! Trigger projection of mid-side nodes to element interior
    void MapMidSideNodes() override;
        
    //! Prints information about the grid into the .info.xml file
    void ToInfo(PtrParamNode in);
   
    /** find entity with minimum distance. Either node or element barycenter
     * Either with nanoflann (USE_NANOFLANN=ON by default) or by brute force (comparing all elements).
     * @param eps if we find something withing this distance we return it. For nanoflann ignored.
     * @return first 0-based entity matching coord */
    UInt FindNearestEntity( bool isNode, Vector<Double>& coord, double eps = 1e-9 );

    /** initialize and setup node/elemGridKDTree_ for nodes or elements for nanoflann. Checks initialization. */
    void CreateKDTree(bool isNode);


    //! Correct the connectivity of elements in the grid.
    //! The function iterates over elements, checks for a negative Jacobi determinant, 
    //! and makes some attempts to correct the connectivity to achieve a positive determinant.
    //! \param regionId (input) allows to specify a region to check. Leaving at default will check all regions
    void CorrectElementConnectivities(RegionIdType regionId = -1);

    /** 
     * makes named nodes from line
     * gmesh can not creat named nodes from line as ansys so openCFS needs to do
     * it. It is also capable of exluding nodes which are found in a different
     * line. This may be needed if you want only the interior nodes of a line.
     */
    void makeNameNodesFromLines();

    //@}

    // =======================================================================
    // General attributes
    // =======================================================================
    //@{ \name General Attributes

    //! ID of grid
    std::string gridId_;

    //! Total number of nodes
    UInt numNodes_;

    //! Total number of elements
    UInt numElems_;

    //! Number of dynamically generated surfelems
    UInt numNcSurfElems_;
    
    //! Flag indicating use of quadratic elements
    bool isQuadratic_;


    //@}

    // =======================================================================
    // Mesh attributes
    // =======================================================================
    //@{ \name Mesh Attributes
  
    //! Vector with nodal coordinates
    StdVector<Vector<Double> > coords_;

    /** for USE_NANOFLANN the O(log n) FindNearestEntity() version. */
    NodeGridKDTree* nodeKdTree_ = nullptr; // either from GridKDTree.hh or dummy implementation in GridCFS.cc
    ElemGridKDTree* elemKdTree_ = nullptr;
    bool use_nanoflann_ = false;

    //! Vector with nodal coordinate offsets
    StdVector<Vector<Double> > deltCoords_;
  
    /** Main vector of elements (surface and volume), ordered by element number.
     * All elements will have the extension set */
    StdVector<Elem*> orderedElems_;
  
    //! Map containing number elements of each type
    std::map<Elem::FEType, UInt> numElemTypes_;

    //! Maps from a node number to all neighbor elements
    StdVector<StdVector<Elem*> > mapNodeToElems_;

    //! Flag to ensure that mapNodeToElems_ and mapNodeToElemsNew_ is only set up once
    //! Maximum number of nodes per element
    UInt maxNumElemNodes_;

    //! Indices to search for the elements containing on specific node number in nodeElemMap_
    //! The element connected to node number n are contained in the range:
    //! { nodeElemMap_[nodeElemMapIndices_[n]], nodeElemMap_[nodeElemMapIndices_[n+1]] {
    //! (last element excluded)
    StdVector<UInt> nodeElemMapIndices_;
    
    //! Contains the element numbers, contains the specific
    StdVector<UInt> nodeElemMap_;

    //! Flag to ensure that mapNodeToElems_ is only set up once
    bool mappedNodeToElems_ = false;
    //@}
  
    //! Vector containing all faces 
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
    
    //! Vector with nodes of named elements
    StdVector<StdVector<UInt> > namedElemNodes_;
  
    //! Vector with names of named elements
    StdVector<std::string> namedElemNames_;

    //@}
    shared_ptr<Timer> userDefinedTimer_;
    // the following timers are only added to .info.xml with -d (detailed info)
    shared_ptr<Timer> initKDTreeTimer_;
    shared_ptr<Timer> searchKDTreeTimer_;
    shared_ptr<Timer> checkRegularTimer_; // CheckForRegularRegion()
    shared_ptr<Timer> correctElemTimer_; // CorrectElementConnectivities()
    shared_ptr<Timer> finishInitTimer_;   // FinishInit()
    shared_ptr<Timer> volumeTimer_;  
    shared_ptr<Timer> mapFacesEdgesTimer_;  
  };

} // end of namespace
#endif // FILE_GRID_CFS
