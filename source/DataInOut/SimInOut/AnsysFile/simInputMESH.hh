// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MOD_ANSYSFILE_2002
#define FILE_MOD_ANSYSFILE_2002

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/simInput.hh"
#include "Domain/elem.hh"
#include "General/defs.hh"
#include "General/environment.hh"

namespace CoupledField {

  //! Class for reading in a mesh created by the ANSYS mkmesh-extension.

  //! Class, that is derived from class FileType for reading mesh-input data,
  //! which is produced by Ansys mkmesh-interface. 
class BaseResult;
class Grid;
template <class TYPE> class StdVector;

  class SimInputMESH: virtual public SimInput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimInputMESH(std::string fileName, PtrParamNode inputNode);
    
    //! Destructor
    virtual ~SimInputMESH();

    //@}

    virtual void InitModule();
    
    virtual void ReadHeader();

    virtual void ReadMesh(Grid *mi);
  
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    UInt GetDim();
    
    //! Get total number of nodes in mesh
    UInt GetNumNodes();
    
    //! Get total number of elements in mesh
    UInt GetNumElems( const Integer dim = 0 );
    
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
    //! \note Since the RegionIdType is guaranteed to be defined by
    //! a number type (UInt, UInt), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( StdVector<std::string> & regionNames );
    
    //! Get vector with region names of given dimension
    
    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    void GetRegionNamesOfDim( StdVector<std::string> & regionNames,
                              const UInt dim );
   
    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    void GetNodeNames( StdVector<std::string> & nodeNames );
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    void GetElemNames( StdVector<std::string> & elemNames );

    //@}

    // =========================================================================
    // GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information
    //! Fill pre-initialized results object with values of specified step
/*    virtual void GetResult( UInt sequenceStep,
                            UInt stepValue,
                            shared_ptr<BaseResult> result,
                            bool isHistory = false ) {
      
    }*/ //todo 
    //@}


  protected:
    
    // =======================================================================
    // AUXILLIARY METHODS
    // =======================================================================

    //@{
    //! \name Auxilliary Methods

    //! Returns a region identifier for a given region

    //! This method takes the name of a region and its dimension and returns
    //! the according identifier. The grid grd (if != NULL) is used to obtain
    //! a real regionId, if not given a local set of indices is used.
    //! \param regionName (input) name of region 
    //! \param dim (input) dimension of the elements (1, 2 or 3)
    //! \param grd (input) the grid structure used for getting the regionIds
    RegionIdType ObtainRegionId( const std::string & regionName,
                                 const UInt dim,
                                 Grid* grd);

    
    //! Transform type of elem in pointer to base class BaseFE

    //! This method maps the type number of an element - as given in the 
    //! mesh file - to a pointer to a reference finite element.
    //! \param itype (input) element type number as read in from the mesh
    Elem::FEType AnsysType2ElemType(const UInt itype);
    //@}
    
    //! Forward the file pointer to the next section just at the beginning of the first non-comment line
    std::string ForwardToSection();
    
    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //! Dimension of the mesh
    UInt dim_;

    //! Number of Elements per Dimension
    std::vector<UInt> numElemsDim_;

    //! Total number of nodes
    UInt numNodes_;
    
    //! number of node BCs
    UInt numNodeBC_;
    
    //! number of save nodes
    UInt numSaveNodes_;

    //! all named node names
    StdVector<std::string> nodeNames_;
    
    //! number of named elements
    UInt numNamedElems_;
    
    //! all named element names
    StdVector<std::string> elemNames_;

    //! struct with info for a region
    struct RegionInfo {
      std::string name;
      UInt dim;
      RegionIdType id;
    };
    
    //! Vector containing all regions of the mesh
    StdVector<RegionInfo> regions_;
    
    //! Pointer to input file
    std::ifstream inFile_;

    //@}

 #ifdef ADAPTGRID
    //! read 2d elements from the input mesh-file for the Grid_RG 
    void ReadEl4AdaptGrid2d( std::vector<grd::Element*> &elems,
                             std::vector<grd::Vertex*> *vertices,
                             const std::vector<std::string> sd );

    //! read 3d elements from the input mesh-file for the Grid_RG 
    void ReadEl4AdaptGrid3d( std::vector<grd::Element*> &elems,
                             std::vector<grd::Vertex*> *vertices,
                             const std::vector<std::string> sd );

    //! for each element of the mesh in format of Grid_RG set value:
    //! color of subdomain
    //! \param namesd name of the subdomain
    //! \param sd vector with colors of all subdomains
    void SetNumSD(grd::Element *ptEl, const std::string namesd,
                  const std::vector<std::string> sd );
 #endif

  };

}

#endif
