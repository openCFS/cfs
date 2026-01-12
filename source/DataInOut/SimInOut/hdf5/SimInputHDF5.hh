#ifndef FILE_CFS_SIMINPUT_HDF5_HH
#define FILE_CFS_SIMINPUT_HDF5_HH

#include <unordered_map>
#include <DataInOut/SimInput.hh>
#include "H5Cpp.h"


namespace CoupledField {

  class CoordSystem;

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
    SimInputHDF5(std::string fileName, PtrParamNode inputNode,
                 PtrParamNode infoNode );
    
    //! Destructor
    virtual ~SimInputHDF5();

    //! Initialize module with pointer to grid
    virtual void InitModule();

    //! Trigger reading of the mesh
    virtual void ReadMesh(Grid *mi);
    
    //! Return file name including path
    std::string GetFileName() {
      return fileName_;
    }
    //@}
  
    // =======================================================================
    //  GENERAL MESH INFORMATION
    // =======================================================================
    //@{ \name General Mesh Information

    //! Return dimension of the mesh
    UInt GetDim();

    // =========================================================================
    //  GENERAL SOLUTION INFORMATION
    // =========================================================================
    //@{ \name General Solution Information

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
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
                    bool isHist = false);

    //! Read one the strings in the user data group.
    void ReadStringFromUserData(const std::string& dSetName,
                                std::string& str);
    //@}


    // =======================================================================
    //  DATABASE SECTION
    // =======================================================================
    //@{ \name Database Handling

    //! Query, if hdf5 file has database
    void DB_Init();
    
    //! Obtain content of XML file
    void DB_GetParamFileContent( std::string& params );
    
    //! Obtain content of material file
    void DB_GetMatFileContent( std::string& params );
    
    //! Get coefficients for given vector
    void DB_GetFeFctCoefs( UInt sequenceStep, UInt stepNum,
                        const std::string& pdeName,
                        const std::string& funcName,
                        SingleVector * coefs ); 

    //! Return multisequence steps and their analysistypes
    void DB_GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                      std::map<UInt, Double>& accTime,
                                      std::map<UInt, bool>& isFinished );
                                    
    //! Return PDE and CoefFunctions in given multisequence step
    void DB_GetAvailPdeCoefFcts( UInt msStep, 
                                 std::map<std::string, 
                                 std::set<std::string> >& coefFcts );
    
    //! Return list with time / frequency values and step for a given result
    void DB_GetStepValues( UInt sequenceStep,
                           const std::string& pdeName,
                           const std::string& resultName,
                           std::map<UInt, Double>& stepValues );

    //! Close database group
    void DB_Close();

    //@}


    void GetNamedNodeResult(const std::string& nodeName,
                            const std::string resultName,
                            StdVector<Complex>& result);    
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
    void GetMeshResult( UInt sequenceStep, UInt stepNum,
                           shared_ptr<BaseResult> result );
    
    //! Read history result
    void GetHistResult( UInt sequenceStep, UInt stepNum,
                           shared_ptr<BaseResult> result);
    
    void LinearizeElems(const StdVector<UInt>& readElems,
                        StdVector<Integer>& elemTypes, 
                        StdVector<UInt>& globConnect, 
                        StdVector<UInt>& readNodes);

    void TransformNodes(CoordSystem& coordSys, double scaleFac);

    //@}

    // =======================================================================
    //  HDF5 DATA MEMBERS
    // =======================================================================
    //@{ \name HDF5 Data Members 

    //! Main hdf5 file
    H5::H5File mainFile_;

    //! Root group of main file
    H5::Group mainRoot_;
    
    //! Root for database
    H5::Group dbRoot_;
    //@}

    // =======================================================================
    //  CLASS ATTRIBUTES
    // =======================================================================
    //@{ \name Attributes
    
    //! Flag indicating if mesh meta data is already read in
    bool statsRead_;
    
    //! Flag indicating use of external files for mesh results
    bool hasExternalFiles_;

    //! Flag for creating named nodes for each region
    bool genRegionNodes_;
    
    //! Flag indicating if the complete grid is to be loaded
    bool readAllEntities_;
    
    //! Native directory path to hdf5 file
    std::string baseDir_;

    //! List of entities (regions & named nodes/elems) to be read in from the file
    std::set< std::string > readEntities_;

    //! List of entities (regions & named elems) to be linearized (quadr -> lin)
    std::set< std::string > linearizeEntities_;

    //! List with names of regions
    StdVector< std::string > regionNames_;

    //! Map with number of dimensions for each region
    std::map<std::string, UInt> regionDims_;

    //! List with names of node groups
    StdVector< std::string > nodeNames_;

    //! List with names of element groups
    StdVector< std::string > elemNames_;

    // Number of nodes in mesh file
    UInt numNodes_;
    
    // Node coordinates
    StdVector<Double> nodeCoords_;
    
    // Map from mesh file node numbers to grid node numbers
    std::unordered_map<UInt, UInt> f2GNodeNumMap_;

    // Map from grid node numbers to mesh file numbers
    std::unordered_map<UInt, UInt> g2FNodeNumMap_;
    
    // Map from mesh file elem numbers to grid elem numbers
    std::unordered_map<UInt, UInt> f2GElemNumMap_;

    // Map from grid elem numbers to mesh file elem numbers
    std::unordered_map<UInt, UInt> g2FElemNumMap_;
    
    // Map from grid entity nodes to indices of mesh entity nodes
    std::unordered_map<std::string, StdVector<UInt> > entityNodeMap_;    

    // Coordinate system into which the node coordinates should be mapped
    std::string coordSysId_;
    
    // Scale factor for node coordinates
    Double scaleFac_;
  };

} // end of namespace

#endif
