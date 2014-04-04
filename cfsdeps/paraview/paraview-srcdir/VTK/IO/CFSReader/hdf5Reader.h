// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HDF5_READER_HH
#define FILE_HDF5_READER_HH

#include <vector>
#include <cpp/H5Cpp.h>

#include "hdf5Common.h"

#ifdef WIN32
#define EXPORT_CLASS __declspec(dllexport)
#else
#define EXPORT_CLASS
#endif


namespace H5CFS {

  //! Class for reading in mesh and simulation data from hdf5 file

  //! Class for handling the reading of mesh and simulation data from 
  //! HDF5 files.
  // class EXPORT_CLASS Hdf5Reader { // just for Windows at the moment
  class Hdf5Reader {

  public:

    // =======================================================================
    //  CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    Hdf5Reader( );
    
    //! Destructor
    virtual ~Hdf5Reader();

    //! Initialize module 
    void LoadFile( const std::string& fileName );

    /** //! Trigger reading of the mesh
    virtual void ReadMesh(Grid *mi);
    //@}
  **/
    // =======================================================================
    //  GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    unsigned int GetDim();
    
    //! Get total number of nodes in mesh
    unsigned int GetNumNodes();
    
    //! Get total number of elements in mesh
    unsigned int GetNumElems( const int );
    
    //! Get total number of regions
    unsigned int GetNumRegions();

    //! Get total number of named nodes
    unsigned int GetNumNamedNodes();

    //! Get total number of named elements
    unsigned int GetNumNamedElems();
    //@}
    
    // =======================================================================
    //  MESH ENTITY ACCESS
    // =======================================================================
    //@{ \name Accessing mesh information
    
    //! Get all nodal coordinates
    void GetNodeCoords( std::vector<std::vector<double> >& coords );
    
    //! Get nodal coordinate of one node
    void GetNodeCoord( unsigned int nodeNum, std::vector<double>& coord );
    
    //! Get all element definitions
    void GetElems( std::vector<ElemType>& type, 
                   std::vector<std::vector<unsigned int> >& connect );
    
    //! Get element definition of single element
    void GetElem( unsigned int elemNum, H5CFS::ElemType& type, 
                  std::vector<unsigned int>& connect );
                  
    //! Get elements of specific region
    void GetElemsOfRegion( const std::string& regionName,
                           std::vector<unsigned int>& elemNums );
                           
    //! Get nodes of specific region
    void GetNodesOfRegion( const std::string& regionName,
                           std::vector<unsigned int>& nodeNums );

    //! Get nodes of named node group
    void GetNamedNodes(const std::string& name,
                       std::vector<unsigned int>& nodeNums);
                           
    //! Get elems of named elem group
    void GetNamedElems(const std::string& name,
                       std::vector<unsigned int>& elemNums);
  
    //@}
    // =======================================================================
    //  ENTITY NAME ACCESS
    // =======================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
    
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the regionIdType is guaranteed to be defined by
    //! a number type (unsigned int, unsigned int32), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( std::vector<std::string> & regionNames );
    
    //! Get vector with region names of given dimension

    //! Returns a vector with the names of regions of a given dimension.
    //! This makes it possible to get for example all names of 
    //! 3D, 2D or 1D elements.
    //! \param regionNames (output) vector containing names of regions
    //! \param dim (input) dimension of the region (1,2, or 3)
    void GetRegionNamesOfDim( std::vector<std::string> & regionNames,
                                      const unsigned int dim ) const;

    //! Get vector with all names of named nodes

    //! Returns a vector which contains all names of named nodes.
    //! \param nodeNames (output) vector with names of named nodes
    void GetNodeNames( std::vector<std::string> & nodeNames );
  
    //! Get vector with all names of named elements

    //! Returns a vector which contains all names of element groups
    //! \param elemNames (output) vector with names of element groups
    void GetElemNames( std::vector<std::string> & elemNames );
    
    
    //! Get entities (nodes, elements), on which a result is defined on
    void GetEntities( EntityType type, const std::string& name,
                      std::vector<unsigned int>& entities );
    
    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<unsigned int, AnalysisType>& analysis,
                                   std::map<unsigned int, unsigned int>& numSteps,
                                   bool isHistory = false );

    //! Obtain list with result types in each sequence step
    void GetResultTypes( unsigned int sequenceStep, 
                         std::vector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );
    
    //! Return list with time / frequency values and step for a given result
    virtual void GetStepValues( unsigned int sequenceStep,
                                shared_ptr<ResultInfo> info,
                                std::map<unsigned int, double>& steps,
                                bool isHistory = false );

    //! Fill pre-initialized results object with values of specified step
    void GetResult( unsigned int sequenceStep,
                    unsigned int stepNum,
                    shared_ptr<Result> result,
                    bool isHistory = false );
    //@}
  protected:

    // =======================================================================
    //  HELPER METHODS
    // =======================================================================
    //@{ \name Helper methods

    //! Read nodal and element definitions
    void ReadNodeElemData(const H5::Group& meshGroup);

    //! Read node groups
    void ReadNodeGroups(const H5::Group& meshGroup);

    //! Read element groups
    void ReadElemGroups(const H5::Group& meshGroup);
  
    //! Read meta information about grid
    void ReadMeshStats(const H5::Group& meshGroup);
    
    //! Read mesh result
    void GetMeshResult( unsigned int sequenceStep, 
                        unsigned int stepNum,
                        shared_ptr<Result> result );
    //! Read history result
    void GetHistResult( unsigned int sequenceStep, 
                        unsigned int stepNum,
                        shared_ptr<Result> result );
//    
//    void LinearizeElems(const std::vector<unsigned int>& readElems,
//                        std::vector<int>& elemTypes, 
//                        std::vector<unsigned int>& globConnect, 
//                        std::vector<unsigned int>& readNodes);
    //@}

    // =======================================================================
    //  HDF5 DATA MEMBERS
    // =======================================================================
    //@{ \name HDF5 Data Members 

    //! Main hdf5 file
    H5::H5File mainFile_;

    //! Root group of main file
    H5::Group mainRoot_;
    
    //! Root group for mesh section
    H5::Group meshRoot_;
    
    //@}

    // =======================================================================
    //  CLASS ATTRIBUTES
    // =======================================================================
    //@{ \name Attributes
    
    //! Filename
    std::string fileName_;
    
    //! Native directory path to hdf5 file
    std::string baseDir_;

    //! Flag inicating if mesh meta data is already read in
    bool statsRead_;
    
    //! Flag indicating use of external files for mesh results
    bool hasExternalFiles_;

    //! List with names of regions
    std::vector< std::string > regionNames_;

    //! Map with number of dimensions for each region
    std::map<std::string, unsigned int> regionDims_;

    //! Map with element numbers for each region
    std::map<std::string, std::vector<unsigned int> > regionElems_;

    //! Map with node numbers for each region
    std::map<std::string, std::vector<unsigned int> > regionNodes_;

    //! List with names of node groups
    std::vector< std::string > nodeNames_;

    //! List with names of element groups
    std::vector< std::string > elemNames_;

    //! Map with element numbers for each named element group
    std::map<std::string, std::vector<unsigned int> > entityElems_;

    //! Map with element numbers for each named element and node group
    std::map<std::string, std::vector<unsigned int> > entityNodes_;

    // Number of nodes in mesh file
    unsigned int numNodes_;

    // Number of elements in mesh file
    unsigned int numElems_;
    
    // Node coordinates
    std::vector<double> nodeCoords_;
    
  };

} // end of namespace

#endif
