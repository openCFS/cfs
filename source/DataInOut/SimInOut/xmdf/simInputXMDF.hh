// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUTXMDF_2006
#define FILE_SIMINPUTXMDF_2006

#include <set>

#include <DataInOut/simInput.hh>

#include "H5Cpp.h"

namespace CoupledField {

  //! Class for reading in a mesh created by the ANSYS mkmesh-extension.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by Ansys mkmesh-interface. 
  class SimInputXMDF: virtual public SimInput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimInputXMDF(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimInputXMDF();

    //@}

    virtual void InitModule(Grid *mi);

    virtual void ReadMesh();
  
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    UInt GetDim();
    
    //! Get total number of nodes in mesh
    UInt GetNumNodes();
    
    //! Get total number of elements in mesh
    UInt GetNumElems( const Integer );
    
    //! Get total number of regions
    UInt GetNumRegions();

    //! Get total number of named nodes
    UInt GetNumNamedNodes();

    //! Get total number of named elements
    UInt GetNumNamedElems();

    //@}
  
    // =======================================================================
    // ENTITY NAME ACCESS
    // =======================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
    
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the regionIdType is guaranteed to be defined by
    //! a number type (UInt, uint32), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( std::vector<std::string> & regionNames );
    
    //! Get vector with region names of given dimension

    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements.
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                      const UInt dim );

    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    virtual void GetNodeNames( StdVector<std::string> & nodeNames );
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    virtual void GetElemNames( StdVector<std::string> & elemNames );

    //! This method returns the root Group of the main HDF5 file.
    H5::Group& GetMainRoot() 
    {
      return mainRoot_;
    }
    
    void GetNumberedRegionName(std::string regionName,
                               std::string& numberedRegName);
    
  protected:
    
    typedef std::vector< std::vector<UInt> > regionElemType;
    typedef std::vector< std::set<UInt, std::less<UInt>, std::allocator<UInt> > > regionNodeType;

    void ReadRegions(const H5::Group& meshGroup);
    void ReadNamedNodes(const H5::Group& meshGroup);
    void ReadNamedElems(const H5::Group& meshGroup);

    //@}

    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //@}

  private:
    std::vector< UInt > regionDims_;
    std::vector< std::string > regionNames_;
    std::vector< std::string > nodeNames_;
    std::vector< std::string > elemNames_;
    std::vector< std::string > readRegions_;
    bool statsRead_;
    bool genRegionNodes_;
    UInt numRegions_;
    H5::Group mainRoot_;
    H5::H5File mainFile_;
    UInt multiStep_, step_;
    bool msChange_;

    void ReadMeshStats(const H5::Group& meshGroup);
    
    typedef struct region_desc_type { 
      char name[32]; 
      UInt dim; 
    };

    typedef struct named_entity_desc_type { 
      char name[32]; 
    };
    
    Integer fg_nElems;
    Integer fg_nNodes;
    Integer fg_nNodesPerElem;
    std::vector<Integer> fg_ElemTypes;
    std::vector<double> fg_XNodeLocs, fg_YNodeLocs, fg_ZNodeLocs;
    std::vector<Integer> fg_NodesInElem;

    typedef std::set<UInt,
                     std::less<UInt>,
                     std::allocator<UInt> > EntitySet;

    //! Set to store all nodes in readRegions_
    EntitySet readNodeSet_;

    //! Set to store all elements in readRegions_
    EntitySet readElemSet_;

    //! readNodeMap_ maps the node numbers from the file
    //! to actual node numbers in the grid.
    std::map<UInt, UInt> readNodeMap_;

    //! readElemMap_ maps the elem numbers from the file
    //! to actual element numbers in the grid.
    //! Needed for named elements.
    std::map<UInt, UInt> readElemMap_;
    
  };

}

#endif
