// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SIMINPUT_REFELEMS_2012
#define FILE_SIMINPUT_REFELEMS_2012

#include <set>

#include <DataInOut/SimInput.hh>

namespace CoupledField {

  //! Class for generating all type of reference elements.

  class SimInputRefElems: virtual public SimInput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimInputRefElems(std::string fileName, PtrParamNode inputNode,
                     PtrParamNode infoNode);
    
    //! Destructor
    virtual ~SimInputRefElems();

    //@}

    virtual void InitModule();

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

    // =======================================================================
    // ENTITY ACCESS
    // =======================================================================
    //@{ \name Entity Access
    
    //! Get all nodal coordinates from 3D grid
    
    //! This method reads all nodal coordinates into a vector of 3D-Points.
    //! \param nodeCoords (output) vector containing nodal coordinates
    void GetCoordinates( std::vector< double > & nodeCoords );
    
    //! Read all elements of given dimension

    //! This method reads all elements of a given dimension (1D, 2D or 3D).
    //! The output is a vector of vectors, where the outer index corresponds
    //! to the different regions and the inner one to the different elements
    //! per region.
    //! \param elems (output) vector containing vectors of pointers to elements
    //!                       per region. The access is like \c elems 
    //!                       \c [regionNr] \c [elemeNr]
    //! \param regionId (output) vector containing the region Ids of the
    //!                          elements corresponding to the outer index in 
    //!                          the elems vector
    //! \param dim (input) dimension of the elements to be read (1,2 or 3)
    void GetElements( std::vector< std::vector<UInt> > & elems,
                      std::vector< std::vector<Elem::FEType> > & elemTypes,
                      std::vector< std::vector<UInt> > & elemNums,
                      std::vector<RegionIdType> & regionIds,
                      const UInt dim );

    //! Read all named nodes
    
    //! This method reads in all named nodes with their according names.
    //! \param nodes (output) vector containing node numbers for each region.
    //!                       The access is like \c nodes \c [nameNr] 
    //!                       \c [nodeNr]
    //! \param nodeNames (output) vector containing the corresponding
    //!                           node names 
    void GetNamedNodes(StdVector<StdVector<UInt> > & nodes,
                       StdVector<std::string> & nodeNames );

    //! Read all named elements

    //! This method reads in all named elements with their according names.
    //! \param elems (output) vector containing node numbers for each region.
    //!                       The access is like \c elems \c [nameNr] 
    //!                       \c [elemNr]
    //! \param elemNames (output) vector containing the corresponding
    //!                           element names 
    void GetNamedElems(StdVector<StdVector<UInt> > &elems,
                       StdVector<std::string> &elemNames );
    //@}
    

    // =========================================================================
    // GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information
    //! Fill pre-initialized results object with values of specified step
    virtual void GetResult( UInt sequenceStep,
                            UInt stepValue,
                            shared_ptr<BaseResult> result,
                            bool isHistory = false ) {
      
    } 
    //@}


  protected:
    
    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //! Dimension of the mesh
    UInt dim_;

    //! Total number of elements
    UInt maxNumElems_;

    //! Total number of nodes
    UInt maxNumNodes_;

    //! Type of reference element
    Elem::FEType feType_;
    //@}

  };

}

#endif
