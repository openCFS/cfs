#ifndef FILE_GRID_INTERFACE_ADAPTIVE_2005
#define FILE_GRID_INTERFACE_ADAPTIVE_2005

#include "Domain/grid.hh"
#include "MultilevelGrid.hh"

namespace CoupledField
{

  // Forward class declaration
  class FileType;


  //! Class for making an interface from base class grid to Adaptive grid

  //! This class serves as a adapter between the base class Grid and
  //! the implementation GridAdaptive
  template<UInt DIM> 
  class GridInterfaceAdaptive: public Grid
  {
  public:
    //! Constructor with parameter 
    GridInterfaceAdaptive(FileType * const aptFileType);

    //! Destructor
    virtual ~GridInterfaceAdaptive() { if (ptGridAdaptive) delete ptGridAdaptive;}
  
    //! Read in Mesh data
    void Read() {
      return;
      // ptGridAdaptive->Read();

      // Transfer region Names to base class
      //ptGridAdaptive->GetAllRegionNames(regionNames_);
    }


    // ======================================================
    // GENERAL GRID INFORMATION
    // ======================================================
    //@{ \name General Grid Information

    //! Return dimension of mesh
    UInt GetDim() {
      //return ptGridAdaptive->GetDim(); 
    }
  
    //! Return maximum number of nodes
    UInt GetNumNodes() {
      //return ptGridAdaptive->GetNumNodes(); 
    }
  
    //! Returns the number of nodes contained in given region
    UInt GetNumNodes( const StdVector<RegionIdType> & regions ) {
      //return ptGridAdaptive->GetNumNodes(regions); 
    }
  
    //! Returns the number of nodes in the given nodelist
    UInt GetNumNodes( const std::string & nodesName ) {
      //return ptGridAdaptive->GetNumNodes(nodesName);
    }

    //! Returns the total number of elements in the grid
    UInt GetNumElems() {
      //return ptGridAdaptive->GetNumElems();
    }
      
    //! Return maximum number of volume elements 
    UInt GetNumVolElems() {
      //return ptGridAdaptive->GetNumVolElems(); 
    }
    //! Return maximum number of surface elements 
    UInt GetNumSurfElems() {
      //return ptGridAdaptive->GetNumSurfElems(); 
    }
  
    //! Returns number of element contained in given regions
    UInt GetNumElems( const StdVector<RegionIdType> & regions ) {
      //return ptGridAdaptive->GetNumElems(regions);
    }
    
    //! Get vector with all region identifiers
    void GetRegionIds( StdVector<RegionIdType> & regions ) {
      //ptGridAdaptive->GetRegionIds(regions); 
    }

    //! Get vector with all volume region identifiers
    void GetVolRegionIds( StdVector<RegionIdType> & volRegions ) {
      //ptGridAdaptive->GetVolRegionIds(volRegions); 
    }
  
    //! Get vector with all surface region identifiers
    void GetSurfRegionIds( StdVector<RegionIdType> & surfRegions ) {
      //ptGridAdaptive->GetSurfRegionIds(surfRegions);
    }
    
    //! Get list with names of all named nodes
    void GetListNodeNames( StdVector<std::string> & nodeNames) {
      //ptGridAdaptive->GetListNodeNames(nodeNames);
    }

    //! Get list with names of all named elements
    void GetListElemNames( StdVector<std::string> & elemNames) {
      //ptGridAdaptive->GetListElemNames(elemNames);
    }
    //@}


    // ======================================================
    // NODE ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Node Access Functions

    //! Get list of nodes by their name
    void GetNodesByName( StdVector<UInt> & nodeList,
                         const std::string & name ) {
      //ptGridAdaptive->GetNodesByName(nodeList, name); 
    }

    //! Get list of nodes contained in a region
    void GetNodesByRegion( StdVector<UInt> & nodeList,
                           const RegionIdType regionId ) {
      //ptGridAdaptive->GetNodesByRegion(nodeList,regionId);
    }
    
    //! Get coordinates of node with global number inode
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Point<DIM> & rfPoint,
                            const UInt inode ) {
      //ptGridAdaptive->GetNodeCoordinate(rfPoint, inode);
    }

    //! Get coordinates of node with global number inode as vector
    //! \param rfPoint (output) coordinates of point 2D
    //! \param inode (input) node number
    void GetNodeCoordinate( Vector<Double> & rfPoint,
                            const UInt inode ) {
      //ptGridAdaptive->GetNodeCoordinate(rfPoint, inode);
    }
    //@}
  
    // ======================================================
    // ELEMENT ACCESS FUNCTIONS
    // ======================================================
    //@{ \name Element Access Functions
    
    
    //! Get element with given element number
    const Elem * GetElem( UInt elemNr ) {
      //return ptGridAdaptive->GetElem(elemNr);
    }

    //! Get list of elements (surface / volumes)
    void GetElems( StdVector<Elem*> & elems, 
                   const RegionIdType regionId ) {
      //ptGridAdaptive->GetElems(elems, regionId);
    }

    //! Get list of volume elements
    void GetVolElems( StdVector<Elem*> & elems, 
                      const RegionIdType regionId ) {
      //ptGridAdaptive->GetVolElems(elems, regionId);
    }
  
    //! Get list of surface elements
    void GetSurfElems( StdVector<SurfElem*> & surfElems, 
                       const RegionIdType regionId ) {
      //ptGridAdaptive->GetSurfElems(surfElems, regionId); 
    }

    //! Get list of elements by their names
    void GetElemsByName( StdVector<Elem*> & elems,
                         const std::string & elemsName ) {
      //ptGridAdaptive->GetElemsByName(elems, elemsName);
    }
  
    //! Get node numbers of given element
    void GetElemNodes( StdVector<UInt> & connect, 
                       const UInt iElem ) {
      //ptGridAdaptive->GetElemNodes(connect, iElem);
    }
  
    //! Get coordinates of element nodes
    void GetElemNodesCoord( Matrix<Double> & coordMat,  
                            const StdVector<UInt> & connect ) {
      //ptGridAdaptive->GetElemNodesCoord(coordMat, connect);
    }
  
    //! Get elements associated with given nodes
    void GetElemsNextToNodes( StdVector<Elem*> & elemList, 
                              const StdVector<UInt> & nodeList,
                              const StdVector<RegionIdType> 
                              & regionIds ) {
      //ptGridAdaptive->GetElemsNextToNodes(elemList, nodeList, regionIds);
    }

    //! Get volume elements lying next to given surface elements
    void GetElemsNextToSurface( StdVector<Elem*> & neighbours, 
                                const StdVector<Elem*> & surfElems,
                                const StdVector<RegionIdType> 
                                & neighRegions ) {
      //ptGridAdaptive->GetElemsNextToSurface(neighbours, surfElems, 
      //                                      neighRegions );
    }
    
    //@}

    // =======================================================================
    // GEOMETRY CALCULATION
    // =======================================================================
    //@{ \name Geometry Calculation
    
    
    //! Returns surface element normal without defined orientation

    void CalcSurfNormal( Vector<Double> & n, 
                         const Elem & surfElem ) {
      //ptGridAdaptive->CalcSurfNormal(n, surfElem);
    }
    
    //! Returns surface element normal with defined orientation
    
    void CalcSurfNormalOutOfVol( Vector<Double> & n,
                                 const Elem & surfElem,
                                 const Elem & volElem ) {
      //ptGridAdaptive->CalcSurfNormalOutOfVol(n, surfElem, volElem);
    }

    //! Returns the volume of a given region
    Double CalcVolumeOfRegion( const RegionIdType regionId, 
                               Boolean isaxi) {
      //return ptGridAdaptive->CalcVolumeOfRegion(regionId, isaxi);
    }
      
    //@}
    

    // ======================================================
    // MISCELLANEOUS
    // ======================================================
    //@{ \name Miscellaneous  
 

    //! Returns node numbers of a list of Elements
    void GetNodesOfElemList( StdVector<UInt> & nodeList,
                             const StdVector<Elem*> & elemList,
			     Boolean onlyLinNodes = FALSE) {
      //ptGridAdaptive->GetNodesOfElemList(nodeList, elemList, onlyLinNodes);
    }
    
  
  protected:

    void GetAllRegionNames( StdVector<std::string> & regionNames ) {
      //ptGridAdaptive->GetAllRegionNames(regionNames);
    }

  private:
    grd::MultilevelGrid * ptGridAdaptive;
    ///
  };

  template<UInt DIM>
  inline GridInterfaceAdaptive<DIM>::GridInterfaceAdaptive(FileType * aptFileType)
    : Grid(aptFileType)
  {
    ENTER_FCN( "GridInterfaceAdaptive<Dim>::GridInterfaceAdaptive<Dim>" );
    //ptGridAdaptive=new GridAdaptive<DIM>(ptFileType);
    ptGridAdaptive=new grd::MultilevelGrid;
  }

#if defined(__GNUC__) || defined(__sgi)
  template class GridInterfaceAdaptive<3>;
  template class GridInterfaceAdaptive<2>;
#endif

} // end of namespace
#endif // FILE_GRID
