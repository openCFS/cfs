#ifndef FILE_SIMINPUT_2006
#define FILE_SIMINPUT_2006

#include <string>
#include <vector>

#include "General/environment.hh"
#include "Utils/tools.hh"
#include "Domain/grid.hh"
#include "General/exception.hh"

namespace CoupledField
{

#define MESHIO_ERROR_EX(STR, FILE, LINE){                       \
    (*error) << STR;                                            \
    CoupledField::Error( FILE, LINE);                           \
  }

#define MESHIO_ERROR(STR){                                      \
    (*error) << STR;                                            \
    CoupledField::Error( __FILE__, __LINE__ );                  \
  }
                
#define MESHIO_INFO_EX(STR, FILE, LINE){                       \
    (*debug) << STR;                                           \
  }

#define MESHIO_INFO(STR){                                       \
    (*debug) << STR;                                            \
  }

#define MESHIO_WARN_EX(STR, FILE, LINE){                         \
    (*warning) << STR;                                           \
    CoupledField::Warning( FILE, LINE );                         \
  }

#define MESHIO_WARN(STR){                                        \
    (*warning) << STR;                                           \
    CoupledField::Warning( __FILE__, __LINE__ );                 \
  }
                  
#define MESHIO_DEBUG_EX(STR, FILE, LINE){                        \
    (*debug) << STR;                                             \
  }                                                              \

#define MESHIO_DEBUG(STR){                                       \
    (*debug) << STR;                                             \
  }                                                              \
    

 


  //! Abstract base class for hanling exceptions and errors
  class ErrorHandler {

  public:
    //! Constructor
    ErrorHandler() {};

    //! Destructor
    virtual ~ErrorHandler() {};

    //! Error (non-recoverable)
    virtual void Error( const Exception &exc ) = 0;

    //! Warning (recoverable)
    virtual void Warning( const Exception &exc ) = 0;
  };

  extern std::auto_ptr<ErrorHandler> errHandler;  

  //! Abstract class for reading in mesh data

  //! This class defines an abstract interface for accessing 
  //! files containing geometric mesh information. 
  //!
  //! \note All mesh and geometric entities are counted one-based,
  //! whereas the acces to the C++ built in datatypes is zero-based!
  
  class SimInput
  {

  public:

    // ========================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // ========================================================================
    //@{ \name Constructor / Initialization
  

    //! Constructor with name of mesh-file
    SimInput(std::string fileName) : fileName_(fileName) {};

    //! Destructor
    virtual ~SimInput() {};

    //@}

    virtual void InitModule(Grid *mi) = 0;
      
    virtual void ReadMesh() = 0;

    // ========================================================================
    // GENERAL MESH INFORMATION
    // ========================================================================
    //@{ \name General Mesh Information

    //! Get dimension of the mesh
    virtual UInt GetDim() = 0;

    //! Get total number of nodes in mesh
    virtual UInt GetNumNodes() = 0;
 
    //! Get total number of elements in mesh
    virtual UInt GetNumElems( const int32_t dim = -1 ) = 0;

    //! Get total number of regions
    virtual UInt GetNumRegions() = 0;

    //! Get total number of named nodes
    virtual UInt GetNumNamedNodes() = 0;

    //! Get total number of named elements
    virtual UInt GetNumNamedElems() = 0;

    //@}
  
    // =========================================================================
    // ENTITY NAME ACCESS
    // =========================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
 
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the RegionIdType is guaranteed to be defined by
    //! a number type (UInt, uint32), the regionId of an element can
    //! be directly used as index to the regionNames-vector
    virtual void GetAllRegionNames( std::vector<std::string> & regionNames ) = 0;

    //! Get vector with region names of given dimension

    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements.
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    virtual void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                                      const UInt dim ) = 0;

    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    virtual void GetNodeNames( StdVector<std::string> & nodeNames ) = 0;
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    virtual void GetElemNames( StdVector<std::string> & elemNames ) = 0;

    //@}

  protected:

    //! Name of input file
    std::string fileName_;

    std::string baseDir_;
    std::string baseName_;
    Grid *mi_;

    UInt dim_;
    std::vector<UInt> numElemsOfDim_;
    std::map<UInt, StdVector<std::string> > regionNamesOfDim_;
    StdVector<std::string> namedElems_;
    StdVector<std::string> namedNodes_;
    
  };

}

#endif // FILE_FILETYPE

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
