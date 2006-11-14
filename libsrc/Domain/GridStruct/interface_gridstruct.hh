#ifndef FILE_GRID_INTERFACE_STRUCT_2004
#define FILE_GRID_INTERFACE_STRUCT_2004

#include "Domain/grid.hh"
#include "grid_struct.hh"
#include "grid_struct.cc"

namespace CoupledField
{
  
  class FileType;
  
  /// Class for working with grid
  template<Integer DIM>
  class GridInterfaceStruct: public Grid
  {
  public:
    /// Constructor with parameter - pointer to FileType for reading initial grid
    GridInterfaceStruct(UInt aDIM);

    /// Deconstructor
    virtual ~GridInterfaceStruct() { if (ptGridStruct) delete ptGridStruct;}

    //! Read of mesh
    virtual void Read()
    { 
      ptGridStruct->Read();

      // Transfer region Names to base class
      ptGridStruct->GetAllRegionNames(regionNames_);
    }

    //! Map element faces
    void MapFaces() {
      ptGridStruct->MapFaces();
    }

    //! Map element edges
    void MapEdges() {
      ptGridStruct->MapEdges();
    }

    //==================================================================================
    //special functions for structured Grid
    //! Read of structured mesh
    virtual void GenGridStruct(const UInt elemx, const UInt elemy,const UInt elemz)
    { ptGridStruct->GenGridStruct(elemx, elemy, elemz);}

    //! Transform Grid
    //  virtual void TransformGridStruct(UInt& nodeShift, UInt& shiftFactor, const UInt flag)
    //     { ptGridStruct->TransformGridStruct(nodeShift, shiftFactor, flag);}
    virtual void TransformGridStruct(UInt & shiftFactor, UInt & nodeShift,
					UInt & elemgrid, Double &  meshsize, const UInt flag)
    { ptGridStruct->TransformGridStruct(shiftFactor, nodeShift,
					elemgrid, meshsize, flag);}

    //Get Maximum nuber of Elements in x,y,z-direction
    virtual Integer GetMaxElem(const std::string dir)
    {return ptGridStruct->GetMaxElem(dir);}

    //===================================================================================

    // ======================================================
    // GENERAL GRID INFORMATION
    // ======================================================
    //@{ \name General Grid Information

    //! Return if grid uses quadratic elements
    bool IsQuadratic() {
      return ptGridStruct->IsQuadratic();
    }

    //! Return number of elements of a given type
    //! \param type Type of finite element (LINE, TRIA, ...)
    UInt GetNumElemOfType( FEType type ) {
      return ptGridStruct->GetNumElemOfType( type );
    }
    
    //! Return dimension of mesh
    UInt GetDim() {
      return ptGridStruct->GetDim(); 
    }
  
    //! Return maximum number of nodes
    UInt GetNumNodes() {
      return ptGridStruct->GetNumNodes(); 
    }
  
    //! Returns the number of nodes contained in given region
    UInt GetNumNodes( const StdVector<RegionIdType> & regions ) {
      return ptGridStruct->GetNumNodes(regions); 
    }
  
    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName ) {
      return ptGridStruct->GetNumNodes(nodesName);
    }

    //! Returns the total number of elements in the grid
    UInt GetNumElems() {
      return ptGridStruct->GetNumElems();
    }

    //! Return maximum number of volume elements 
    UInt GetNumVolElems() {
      return ptGridStruct->GetNumVolElems(); 
    }
    //! Return maximum number of surface elements 
    UInt GetNumSurfElems() {
      return ptGridStruct->GetNumSurfElems(); 
    }
  
    //! Returns number of element contained in given regions
    UInt GetNumElems( const StdVector<RegionIdType> & regions ) {
      return ptGridStruct->GetNumElems(regions);
    }

    //! Get vector with all region identifiers
    void GetRegionIds( StdVector<RegionIdType> & regions ) {
      ptGridStruct->GetRegionIds(regions); 
    }
  
    //! Get vector with all volume region identifiers
    void GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
      ptGridStruct->GetVolRegionIds(volRegions); 
    }
  
    //! Get vector with all surface region identifiers
    void GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {
      ptGridStruct->GetSurfRegionIds(surfRegions);
    }

    //! Get list with names of all named nodes
    void GetListNodeNames( StdVector<std::string> & nodeNames) {
      ptGridStruct->GetListNodeNames(nodeNames);
    }
    
    //! Get list with names of all named elements
    void GetListElemNames( StdVector<std::string> & elemNames) {
      ptGridStruct->GetListElemNames(elemNames);
    }
    //@}


    // ======================================================
    // NODE ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Node Access Functions

    //! Get list of nodes by their name
    void GetNodesByName( StdVector<UInt> & nodeList,
                         const std::string & name ) {
      ptGridStruct->GetNodesByName(nodeList, name); 
    }

    //! Get list of nodes contained in a region
    void GetNodesByRegion( StdVector<UInt> & nodeList,
                           const RegionIdType regionId ) {
      ptGridStruct->GetNodesByRegion(nodeList,regionId);
    }
    
    //! Get coordinates of node with global number inode
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Point<DIM> & rfPoint,
                            const UInt inode,
                            bool updated = false ) {
      ptGridStruct->GetNodeCoordinate(rfPoint, inode, updated);
    }

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode,
                            bool updated = false) {
      ptGridStruct->GetNodeCoordinate(rfPoint, inode, updated);
    }
    //@}
  
    // ======================================================
    // ELEMENT ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Element Access Functions
  
    //! Get element with given element number
    const Elem * GetElem( UInt elemNr ) {
      return ptGridStruct->GetElem(elemNr);
    }

    
    //! Get list of elements (surface / volumes)
    void GetElems( StdVector<Elem*> & elems, 
                   const RegionIdType regionId ) {
      ptGridStruct->GetElems(elems, regionId);
    }


    //! Get list of volume elements
    void GetVolElems( StdVector<Elem*> & elems, 
                      const RegionIdType regionId ) {
      ptGridStruct->GetVolElems(elems, regionId);
    }
  
    //! Get list of surface elements
    void GetSurfElems( StdVector<SurfElem*> & surfElems, 
                       const RegionIdType regionId ) {
      ptGridStruct->GetSurfElems(surfElems, regionId); 
    }

    //! Get list of elements by their names
    void GetElemsByName( StdVector<Elem*> & elems,
                         const std::string & elemsName ) {
      ptGridStruct->GetElemsByName(elems, elemsName);
    }
  
    //! Get node numbers of given element
    void GetElemNodes( StdVector<UInt> & connect, 
                       const UInt iElem ) {
      ptGridStruct->GetElemNodes(connect, iElem);
    }
  
    //! Get coordinates of element nodes
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect,
                            bool updated = false) {
      ptGridStruct->GetElemNodesCoord(coordMat, connect, updated);
    }
  
    //! Get elements associated with given nodes
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds ) {
      ptGridStruct->GetElemsNextToNodes(elemList, nodeList, regionIds);
    }

    //! Get volume elements lying next to given surface elements
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                & neighRegions ) {
      ptGridStruct->GetElemsNextToSurface(neighbours, surfElems, 
                                          neighRegions );
    }
    
    //@}

    // =======================================================================
    // ELEMENT FACE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Surface Access Functions

    //! Get total number of faces in the grid
    UInt GetNumFaces() {
      return ptGridStruct->GetNumFaces();
    }
    
    //! Return face with given face
    const Face& GetFace( UInt faceNr) {
      return ptGridStruct->GetFace( faceNr );
    }

    // =======================================================================
    // EDGE ACCESS FUNCTIONS
    // =======================================================================
    //@{ \name Edge Access Functions
          

    //! Get number ofe edges
    UInt GetNumEdges() {
      return ptGridStruct->GetNumEdges();
    }
    
    const Edge& GetEdge( UInt edgeNr ) {
      return ptGridStruct->GetEdge( edgeNr );
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
      ptGridStruct->CalcSurfNormal(n, surfElem, updated );
    }
    
    //! Returns surface element normal with defined orientation
  
    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem,
                                 bool updated = false  ) {
      ptGridStruct->CalcSurfNormalOutOfVol(n, surfElem, volElem,updated);
    }
  
    //! Returns the volume of a given region
    Double CalcVolumeOfRegion( const RegionIdType regionId, 
                               bool isaxi = false,
                               bool updated = false ) {
      return ptGridStruct->CalcVolumeOfRegion(regionId,isaxi,updated);
    }
  
    //@}
    

    // ======================================================
    // MISCELLANEOUS
    // ======================================================
    //@{ \name Miscellaneo2us  
 

    //! Returns node numbers of a list of Elements
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     bool onlyLinNodes = false) {
      ptGridStruct->GetNodesOfElemList(nodeList, elemList, onlyLinNodes);
    }

    
    //! Set offset for coordinates due to updated Lagrangian formulation
    void SetNodeOffset( const StdVector<UInt>& nodes, 
                        const Vector<Double>& offsets ) {
      ptGridStruct->SetNodeOffset( nodes, offsets );
    }

    bool HasNodalOffset() {
      return ptGridStruct->HasNodalOffset();
    }
  
  protected:

    void GetAllRegionNames( StdVector<std::string> & regionNames ) {
      ptGridStruct->GetAllRegionNames(regionNames);
    }


  protected:


  private:
    GridStruct<DIM> * ptGridStruct;
    ///
  };

  template<Integer DIM>
  inline GridInterfaceStruct<DIM>::GridInterfaceStruct(UInt adim)
    : Grid(NULL)
  {
    ENTER_FCN( "GridInterfaceStruct<DIM>::GridInterfaceStruct<DIM>" );
    //  lastlevel_=0;
    ptGridStruct=new GridStruct<DIM>(adim);
  }

#if defined(__GNUC__) || defined(__sgi)
  template class GridInterfaceStruct<3>;
  template class GridInterfaceStruct<2>;
#endif

} // end of namespace
#endif // FILE_GRID
