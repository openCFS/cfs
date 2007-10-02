// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMINPUT_HDF5_HH
#define FILE_CFS_SIMINPUT_HDF5_HH

#include <DataInOut/simInput.hh>
#include "H5Cpp.h"

namespace CoupledField {

  //! Class for reading in mesh and simulation data from hdf5 file

  //! Class for handling the reading of mesh and simulation data from 
  //! HDF5 files.
  class SimInputHDF5: virtual public SimInput {

  public:

    // =======================================================================
    //  CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimInputHDF5(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimInputHDF5();

    //! Initialize module with pointer to grid
    virtual void InitModule();

    //! Trigger reading of the mesh
    virtual void ReadMesh(Grid *mi);
    //@}
  
    // =======================================================================
    //  GENERAL MESH INFORMATION
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
    //  ENTITY NAME ACCESS
    // =======================================================================
    //@{ \name Entity Name Access
  
    //! Get vector with all region names in mesh
    
    //! Returns a vector with the names of regions in the mesh of all
    //! dimensions.
    //! \param regionNames (output) vector containing names of regions
    //! \note Since the regionIdType is guaranteed to be defined by
    //! a number type (UInt, uint32), the regionId of an element can
    //! be directly used as index to the regions-vector
    void GetAllRegionNames( StdVector<std::string> & regionNames );
    
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
    
    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<UInt, AnalysisType>& analysis,
                                   std::map<UInt, UInt>& numSteps,
                                   bool isHistory = false );

    //! Obtain list with result types in each sequence step
    void GetResultTypes( UInt sequenceStep, 
                         StdVector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );
    
    //! Return list with time / frequency values and step for a given result
    virtual void GetStepValues( UInt sequenceStep,
                                    shared_ptr<ResultInfo> info,
                                    std::map<UInt, Double>& steps,
                                    bool isHistory = false );
    
    //! Return entitylist the result is defined on
    void GetResultEntities( UInt sequenceStep,
                               shared_ptr<ResultInfo> info,
                               StdVector<shared_ptr<EntityList> >& list,
                               bool isHistory = false );

    //! Fill pre-initialized results object with values of specified step
    void GetResult( UInt sequenceStep,
                      UInt stepNum,
                      shared_ptr<BaseResult> result,
                      bool isHistory = false );
    //@}

  protected:

    // =======================================================================
    //  HELPER METHODS
    // =======================================================================
    //@{ \name Helper methods

    //! Read elements of regions
    void ReadRegions(const H5::Group& meshGroup);

    //! Read node groups
    void ReadNodeGroups(const H5::Group& meshGroup);

    //! Read element groups
    void ReadElemGroups(const H5::Group& meshGroup);
  
    //! Read meta information about grid
    void ReadMeshStats(const H5::Group& meshGroup);
    
    //! Read mesh result
    void GetMeshResult( UInt sequenceStep, UInt stepNum,
                           shared_ptr<BaseResult> result );
    
    //! Read history result
    void GetHistResult( UInt sequenceStep, UInt stepNum,
                           shared_ptr<BaseResult> result );
    //@}

    // =======================================================================
    //  HDF5 DATA MEMBERS
    // =======================================================================
    //@{ \name HDF5 Data Members 

    //! Main hdf5 file
    H5::H5File mainFile_;

    //! Root group of main file
    H5::Group mainRoot_;
    //@}

    // =======================================================================
    //  CLASS ATTRIBUTES
    // =======================================================================
    //@{ \name Attributes
    
    //! Flag inicating if mesh meta data is already read in
    bool statsRead_;

    //! Flag for creating named nodes for each region
    bool genRegionNodes_;

    //! List of regions to be read in from the file
    std::vector< std::string > readRegions_;

    //! List with names of regions
    std::vector< std::string > regionNames_;

    //! Map with number of dimensions for each region
    std::map<std::string, UInt> regionDims_;

    //! List with names of node groups
    std::vector< std::string > nodeNames_;

    //! List with names of element groups
    std::vector< std::string > elemNames_;
  };

} // end of namespace

#endif
