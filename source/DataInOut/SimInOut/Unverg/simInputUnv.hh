// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUTUNV_2006
#define FILE_SIMINPUTUNV_2006

#include <DataInOut/Logging/cfslog.hh>
#include <DataInOut/simInput.hh>


namespace CoupledField
{

  class MeshInterface;

  // declare logging stream
  DECLARE_LOG(simInputUNV)

  /**
   **/
  class SimInputUnv : public SimInput
  {
  public:
    // ========================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // ========================================================================
    //@{ \name Constructor / Initialization
  

    //! Constructor with name of mesh-file
    SimInputUnv(std::string fileName, ParamNode * inputNode );

    //! Destructor
    virtual ~SimInputUnv() {};

    //@}

    virtual void InitModule();
      
    virtual void ReadMesh(Grid *mi);

    virtual bool ReadData();

    // ========================================================================
    // GENERAL MESH INFORMATION
    // ========================================================================
    //@{ \name General Mesh Information

    //! Get dimension of the mesh
    virtual UInt GetDim();

    //! Get total number of nodes in mesh
    virtual UInt GetNumNodes();
 
    //! Get total number of elements in mesh
    virtual UInt GetNumElems( const Integer dim = -1 );

    //! Get total number of regions
    virtual UInt GetNumRegions();

    //! Get total number of named nodes
    virtual UInt GetNumNamedNodes();

    //! Get total number of named elements
    virtual UInt GetNumNamedElems();

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
    virtual void GetAllRegionNames( StdVector<std::string> & regionNames );

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

    //@}

  private:
    FEType UnvType2ElemType( const uint32_t elemType );

  }; 

} 

#endif 
