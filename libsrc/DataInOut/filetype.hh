#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include <string>

#include "Utils/StdVector.hh"

#ifdef ADAPTGRID
#include "Element.h"  
#endif

namespace CoupledField
{

  // Forward class declaration
  struct Elem;

  //! Abstract class for reading in mesh data

  //! This class defines an abstract interface for accessing 
  //! files containing geometric mesh information. 
  //!
  //! \note All mesh and geometric entities are counted one-based,
  //! whereas the acces to the C++ built in datatypes is zero-based!
  
  class FileType
  {

  public:

    // ========================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // ========================================================================
    //@{ \name Constructor / Initialization
  

    //! Constructor with name of mesh-file
    FileType(const Char * const afilename);

    //! Destructor
    virtual ~FileType();

    //@}
  
    // ========================================================================
    // GENERAL MESH INFORMATION
    // ========================================================================
    //@{ \name General Mesh Information

    //! Get dimension of the mesh
    virtual UInt GetDim() = 0;

    //! Get total number of nodes in mesh
    virtual UInt GetNumNodes() = 0;
 
    //! Get total number of elements in mesh
    virtual UInt GetNumElems( const UInt dim = 0 ) = 0;

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
    //! a number type (UInt, UInt), the regionId of an element can
    //! be directly used as index to the regionNames-vector
    virtual void GetAllRegionNames( StdVector<std::string> & regionNames ) = 0;

  
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

    // =========================================================================
    // ENTITY ACCESS
    // =========================================================================
    //@{ \name Entity Access

    //! Get all nodal coordinates from 3D grid

    //! This method reads all nodal coordinates into a vector of 3D-Points.
    //! \param nodeCoords (output) vector containing nodal coordinates
    virtual void GetCoordinates( StdVector<Point<3> > & nodeCoords ) = 0;

    //! Get all nodal coordinates from 2D grid

    //! This method reads all nodal coordinates into a vector of 2D-Points.
    //! \param nodeCoords (output) vector containing nodal coordinates
    virtual void GetCoordinates( StdVector<Point<2> > & nodeCoords ) = 0;

    //! Get vector of nodes for each region

    //! This method reads the node numbers of each region into a 
    //! separate vector. 
    //! \param nodes (output) vector containing the node numbers for each
    //!                       region. The access is like \c elems 
    //!                       \c [regionNr] \c [nodeNr]
    //! \param regionId (output) vector containing the region Ids of the
    //!                          nodes corresponding to the outer index in the
    //!                          nodes vector
    virtual void GetNodesOfRegions( StdVector<StdVector<UInt> > &nodes,
                                    const StdVector<RegionIdType>  
                                    & regionId ) = 0;
  
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
    virtual void GetElements( StdVector< StdVector<Elem*> > & elems, 
                              StdVector<RegionIdType> & regionId,
                              const UInt dim ) = 0;
  
    //! Read all named nodes

    //! This method reads in all named nodes with their according names.
    //! \param nodes (output) vector containing node numbers for each region.
    //!                       The access is like \c nodes \c [nameNr] 
    //!                       \c [nodeNr]
    //! \param nodeNames (output) vector containing the corresponding
    //!                           node names 
    virtual void GetNamedNodes( StdVector<StdVector<UInt> > & nodes,
                                StdVector<std::string> & nodeNames ) = 0;
    //! Read all named elements

    //! This method reads in all named elements with their according names.
    //! \param elems (output) vector containing node numbers for each region.
    //!                       The access is like \c elems \c [nameNr] 
    //!                       \c [elemNr]
    //! \param elemNames (output) vector containing the corresponding
    //!                           element names 
    virtual void GetNamedElems( StdVector<StdVector<UInt> > & elems,
                                StdVector<std::string> & elemNames ) = 0;
    //@}
                               
  protected:

    //! Name of input file
    Char * filename;

  };

}


#endif // FILE_FILETYPE
