#ifndef FILE_MOD_XMDF_2006
#define FILE_MOD_XMDF_2006

#include <set>

#include <DataInOut/meshio.hh>

#include "Xmdf.h"

namespace CoupledField {

  //! Class for reading in a mesh created by the ANSYS mkmesh-extension.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by Ansys mkmesh-interface. 
  class XMDF: virtual public MeshIOModule {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    XMDF(std::string fileName);
    
    //! Destructor
    virtual ~XMDF();

    //@}

    virtual void RegisterModule();
  
    virtual void InitModule(const std::string & baseDir,
                            const std::string & baseName,
                            Grid *mi);

    virtual void ReadMesh();
    virtual void WriteMesh();
  
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    uint32_t GetDim();
    
    //! Get total number of nodes in mesh
    uint32_t GetNumNodes();
    
    //! Get total number of elements in mesh
    uint32_t GetNumElems( const uint32_t dim = 0 );
    
    //! Get total number of regions
    uint32_t GetNumRegions();

    //! Get total number of named nodes
    uint32_t GetNumNamedNodes();

    //! Get total number of named elements
    uint32_t GetNumNamedElems();

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
    //! a number type (uint32_t, uint32), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( std::vector<std::string> & regionNames );
    

  protected:
    
    //! Transform type of elem in pointer to base class BaseFE

    //! This method maps the type number of an element - as given in the 
    //! mesh file - to a pointer to a reference finite element.
    //! \param itype (input) element type number as read in from the mesh

    FEType XMDFElemType2ElemType( const int32_t type );
    int32_t ElemType2XMDFElemType( const FEType type );

    void ReorderConnectivity( const int32_t eType,
                              const bool toXMDF,
                              const uint32_t* in,
                              uint32_t* out);

    hid_t CreateGroup(const hid_t parentGroup, const std::string name);

    typedef std::vector< std::vector<uint32_t> > regionElemType;
    typedef std::vector< std::set<uint32_t, std::less<uint32_t>, std::allocator<uint32_t> > > regionNodeType;

    void ReadRegions(const hid_t meshGroup);
    void ReadNamedNodes(const hid_t meshGroup);
    void ReadNamedElems(const hid_t meshGroup);

    void WriteRegions(const hid_t meshGroup,
                      const regionElemType regionElems,
                      const regionNodeType regionNodes);

    void WriteNamedNodes(const hid_t meshGroup);
    void WriteNamedElems(const hid_t meshGroup);

    //@}

    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //@}


  };

}

#endif // FILE_MOD_XMDF_2006

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
