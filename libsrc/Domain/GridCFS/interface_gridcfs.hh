#ifndef FILE_GRID_INTERFACE_CFS_2002
#define FILE_GRID_INTERFACE_CFS_2002

#include "Domain/grid.hh"
#include "grid_cfs.hh"
#include "grid_cfs.cc"

namespace CoupledField
{

  // Forward class declaration
  class FileType;


  //! Class for making an interface from base class grid to CFS grid

  //! This class serves as a adapter between the base class Grid and
  //! the implementation GridCFS
  template<UInt DIM> 
  class GridInterfaceCFS: public Grid
  {
  public:
    //! Constructor with parameter 
    GridInterfaceCFS(FileType * const aptFileType);

    //! Destructor
    virtual ~GridInterfaceCFS() { if (ptGridCFS) delete ptGridCFS;}
  
    //! Read in Mesh data
    void Read() {
      ptGridCFS->Read();

      // Transfer region Names to base class
      ptGridCFS->GetAllRegionNames(regionNames_);
    }

    //! Map element sub-entities
    void MapSubEntities() {
      ptGridCFS->MapSubEntities();
    }


    // ======================================================
    // GENERAL GRID INFORMATION
    // ======================================================
    //@{ \name General Grid Information

    //! Search and add directivity nodes, if defined in input file,
    //! to list of saved nodes     
    void GetNodesForDirectivity();

    //! Return if grid uses quadratic elements
    bool IsQuadratic() {
      return ptGridCFS->IsQuadratic();
    }

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    UInt GetNumElemOfType( FEType type ) {
      return ptGridCFS->GetNumElemOfType( type );
    }

    //! Return dimension of mesh
    UInt GetDim() {
      return ptGridCFS->GetDim(); 
    }
  
    //! Return maximum number of nodes
    UInt GetNumNodes() {
      return ptGridCFS->GetNumNodes(); 
    }
  
    //! Returns the number of nodes contained in given region
    UInt GetNumNodes( const StdVector<RegionIdType> & regions ) {
      return ptGridCFS->GetNumNodes(regions); 
    }
  
    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName ) {
      return ptGridCFS->GetNumNodes(nodesName);
    }

    //! Returns the total number of elements in the grid
    UInt GetNumElems() {
      return ptGridCFS->GetNumElems();
    }
      
    //! Return maximum number of volume elements 
    UInt GetNumVolElems() {
      return ptGridCFS->GetNumVolElems(); 
    }
    //! Return maximum number of surface elements 
    UInt GetNumSurfElems() {
      return ptGridCFS->GetNumSurfElems(); 
    }
  
    //! Returns number of element contained in given regions
    UInt GetNumElems( const StdVector<RegionIdType> & regions ) {
      return ptGridCFS->GetNumElems(regions);
    }
    
    //! Get vector with all region identifiers
    void GetRegionIds( StdVector<RegionIdType> & regions ) {
      ptGridCFS->GetRegionIds(regions); 
    }

    //! Get vector with all volume region identifiers
    void GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
      ptGridCFS->GetVolRegionIds(volRegions); 
    }
  
    //! Get vector with all surface region identifiers
    void GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {
      ptGridCFS->GetSurfRegionIds(surfRegions);
    }
    
    //! Get list with names of all named nodes
    void GetListNodeNames( StdVector<std::string> & nodeNames) {
      ptGridCFS->GetListNodeNames(nodeNames);
    }

    //! Get list with names of all named elements
    void GetListElemNames( StdVector<std::string> & elemNames) {
      ptGridCFS->GetListElemNames(elemNames);
    }
    //@}


    // ======================================================
    // NODE ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Node Access Functions

    //! Get list of nodes by their name
    void GetNodesByName( StdVector<UInt> & nodeList,
                         const std::string & name ) {
      ptGridCFS->GetNodesByName(nodeList, name); 
    }

    //! Get list of nodes contained in a region
    void GetNodesByRegion( StdVector<UInt> & nodeList,
                           const RegionIdType regionId ) {
      ptGridCFS->GetNodesByRegion(nodeList,regionId);
    }
    
    //! Get coordinates of node with global number inode
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Point<DIM> & rfPoint,
                            const UInt inode, bool updated = false ) {
      ptGridCFS->GetNodeCoordinate(rfPoint, inode, updated);
    }

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode,  bool updated = false ) {
      ptGridCFS->GetNodeCoordinate(rfPoint, inode, updated);
    }
    //@}
  
    // ======================================================
    // ELEMENT ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Element Access Functions
    
    
    //! Get element with given element number
    const Elem * GetElem( UInt elemNr ) {
      return ptGridCFS->GetElem(elemNr);
    }

    //! Get list of elements (surface / volumes)
    void GetElems( StdVector<Elem*> & elems, 
                   const RegionIdType regionId ) {
      ptGridCFS->GetElems(elems, regionId);
    }

    //! Get list of volume elements
    void GetVolElems( StdVector<Elem*> & elems, 
                      const RegionIdType regionId ) {
      ptGridCFS->GetVolElems(elems, regionId);
    }
  
    //! Get list of surface elements
    void GetSurfElems( StdVector<SurfElem*> & surfElems, 
                       const RegionIdType regionId ) {
      ptGridCFS->GetSurfElems(surfElems, regionId); 
    }

    //! Get list of elements by their names
    void GetElemsByName( StdVector<Elem*> & elems,
                         const std::string & elemsName ) {
      ptGridCFS->GetElemsByName(elems, elemsName);
    }
  
    //! Get node numbers of given element
    void GetElemNodes( StdVector<UInt> & connect, 
                       const UInt iElem ) {
      ptGridCFS->GetElemNodes(connect, iElem);
    }
  
    //! Get coordinates of element nodes
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect,
                            bool updated = false ) {
      ptGridCFS->GetElemNodesCoord(coordMat, connect, updated);
    }
  
    //! Get elements associated with given nodes
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds ) {
      ptGridCFS->GetElemsNextToNodes(elemList, nodeList, regionIds);
    }

    //! Get volume elements lying next to given surface elements
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                & neighRegions ) {
      ptGridCFS->GetElemsNextToSurface(neighbours, surfElems, 
                                       neighRegions );
    }
    
    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    
    //! Returns surface element normal without defined orientation

    void CalcSurfNormal( Vector<Double> & n, 
                         const Elem & surfElem,
                         bool updated = false ) {
      ptGridCFS->CalcSurfNormal(n, surfElem, updated );
    }
    
    //! Returns surface element normal with defined orientation
    
    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem,
                                 bool updated = false ) {
      ptGridCFS->CalcSurfNormalOutOfVol(n, surfElem, volElem, updated );
    }

    //! Returns the volume of a given region
    Double CalcVolumeOfRegion( const RegionIdType regionId, 
                               bool isaxi,
                               bool updated = false ) {
      return ptGridCFS->CalcVolumeOfRegion(regionId, isaxi, updated );
    }
      
    //@}
    

    // ======================================================
    // MISCELLANEOUS
    // ======================================================
    //@{ \name Miscellaneous  
 

    //! Returns node numbers of a list of Elements
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     bool onlyLinNodes = false) {
      ptGridCFS->GetNodesOfElemList(nodeList, elemList, onlyLinNodes);
    }

    //! Set offset for coordinates due to updated Lagrangian formulation
    void SetNodeOffset( const StdVector<UInt>& nodes, 
                        const Vector<Double>& offsets ) {
      ptGridCFS->SetNodeOffset( nodes, offsets );
    }

    //! Return status of presence of nodal coordinate offsets (up. Lagrange)
    bool HasNodalOffset() {
      return ptGridCFS->HasNodalOffset(); 
    }
  
  protected:

    void GetAllRegionNames( StdVector<std::string> & regionNames ) {
      ptGridCFS->GetAllRegionNames(regionNames);
    }

  private:
    GridCFS<DIM> * ptGridCFS;
    ///
  };

  template<UInt DIM>
  inline GridInterfaceCFS<DIM>::GridInterfaceCFS(FileType * aptFileType)
    : Grid(aptFileType)
  {
    ENTER_FCN( "GridInterfaceCFS<Dim>::GridInterfaceCFS<Dim>" );
    ptGridCFS=new GridCFS<DIM>(ptFileType); 
  }

#if defined(__GNUC__) || defined(__sgi)
  template class GridInterfaceCFS<3>;
  template class GridInterfaceCFS<2>;
#endif


} // end of namespace
#endif // FILE_GRID
