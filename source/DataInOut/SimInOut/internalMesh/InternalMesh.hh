// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_INTERNAL_MESH
#define FILE_INTERNAL_MESH

#include <set>

#include "DataInOut/SimInput.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{
  using std::string;  

  /** class that generates a mesh inside cfs for regular 
   * rectangular/ cubic domains */ 
  class InternalMesh : virtual public SimInput 
  {
  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    InternalMesh(string fileName, PtrParamNode inputNode, PtrParamNode infoNode);
    
    //! Destructor
    virtual ~InternalMesh() {}

    //@}

    virtual void InitModule();

    virtual void ReadMesh(Grid *mi);

    void ReadMeshFile(Grid *mi);

    void ReadMeshNetwork(Grid *mi);
  
    // =======================================================================
    // GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    UInt GetDim() { return dim_; }
    
    //! Get total number of nodes in mesh
    UInt GetNumNodes() { return maxNumNodes_; }
    
    //! Get total number of elements in mesh
    UInt GetNumElems(const Integer dim = 0);
    
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
    void GetAllRegionNames(StdVector<string> &regionNames);
    
    //! Get vector with region names of given dimension
    
    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    void GetRegionNamesOfDim(StdVector<string> &regionNames,
                             const UInt dim);
   
    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    void GetNodeNames(StdVector<string> &nodeNames);
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of named elements.
    //! \param elemNames (output) vector with names of named elements
    void GetElemNames(StdVector<string> &elemNames);

    //@}

    // =======================================================================
    // ENTITY ACCESS
    // =======================================================================
    //@{ \name Entity Access
    
    //! Get vector of nodes for each region

    //! This method reads the node numbers of each region into a 
    //! separate vector. 
    //! \param nodes (output) vector containing the node numbers for each
    //!                       region. The access is like \c elems 
    //!                       \c [regionNr] \c [nodeNr]
    //! \param regionId (output) vector containing the region Ids of the
    //!                          nodes corresponding to the outer index in the
    //!                          nodes vector
    void GetNodesOfRegions(StdVector<StdVector<UInt> > &nodes,
                           const StdVector<RegionIdType> & regionId);
    
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
    void GetElements(StdVector<StdVector<UInt> > & elems,
                     StdVector<StdVector<Elem::FEType> > & elemTypes,
                     StdVector<StdVector<UInt> > & elemNums,
                     StdVector<RegionIdType> & regionIds,
                     const UInt dim );

    //! Read all named nodes
    
    //! This method reads in all named nodes with their according names.
    //! \param nodes (output) vector containing node numbers for each region.
    //!                       The access is like \c nodes \c [nameNr] 
    //!                       \c [nodeNr]
    //! \param nodeNames (output) vector containing the corresponding
    //!                           node names 
    void GetNamedNodes(StdVector<StdVector<UInt> > & nodes,
                       StdVector<string> & nodeNames );

    //! Read all named elements

    //! This method reads in all named elements with their according names.
    //! \param elems (output) vector containing node numbers for each region.
    //!                       The access is like \c elems \c [nameNr] 
    //!                       \c [elemNr]
    //! \param elemNames (output) vector containing the corresponding
    //!                           element names 
    void GetNamedElems(StdVector<StdVector<UInt> > & elems,
                       StdVector<string> & elemNames );
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


  private:
    
    // =======================================================================
    // AUXILLIARY METHODS
    // =======================================================================

    //@{
    //! \name Auxilliary Methods
   
    //! Transform type of elem in pointer to base class BaseFE

    //! This method maps the type number of an element - as given in the 
    //! mesh file - to a pointer to a reference finite element.
    //! \param itype (input) element type number as read in from the mesh
    Elem::FEType AnsysType2ElemType(const UInt itype);
    //@}

    /** Parses the boundary conditions defined in xml file */
    void ParseBoundary(PtrParamNode bdr);

    // =======================================================================
    // CLASS ATTRIBUTES
    // =======================================================================
    //@{
    //! \name Attributes

    //! Dimension of the mesh
    UInt dim_;

    //! Number of elements in x, y and z-direction
    UInt nelems_[3];

    //! Minimal coordinates
    Vector<Double> minimal_;

    //! Maximal coordinates
    Vector<Double> maximal_;
    
    //! Total number of elements
    UInt maxNumElems_;

    //! Total number of nodes
    UInt maxNumNodes_;
    
    //! Vector containing all region names of mesh
    StdVector<string> regionNames_;
    
    //! Vector containing dimension of corresponding regionNames_
    StdVector<UInt> regionDim_;
    
    //! Array indicating if elems of given dimension were read in
    StdVector<bool> elemDimReadIn_;

    //! Vector with nodal numbers for each region
    StdVector<std::set<UInt> > regionNodes_;
    
    //! Pointer to xml tree created from input file
    PtrParamNode xml_;

    /** Info Node base */
    PtrParamNode info_;
    //@}

    // =======================================================================
    // NETWORK ATTRIBUTES
    // =======================================================================
    
    //! Bool to differentiate mesh reading mode
    bool readNetwork_ = false;

    PtrParamNode networkNode_;

  };

}

#endif
