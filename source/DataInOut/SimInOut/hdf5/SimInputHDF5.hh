#ifndef FILE_CFS_SIMINPUT_HDF5_HH
#define FILE_CFS_SIMINPUT_HDF5_HH

#include <unordered_map>
#include <DataInOut/SimInput.hh>
#include <H5Ipublic.h>

namespace CoupledField 
{
  class CoordSystem;

  //! Class for reading in mesh and simulation data from hdf5 file

  //! Class for handling the reading of mesh and simulation data from 
  //! HDF5 files.
  class SimInputHDF5: public SimInput 
  {
  public:
    
    //! Constructor with name of mesh-file
    SimInputHDF5(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode );
    
    virtual ~SimInputHDF5();

    //! Return dimension of the mesh
    unsigned int GetDim() override { return dim_; }

    //! Return multisequence steps and their analysistypes
    void GetNumMultiSequenceSteps( std::map<unsigned, BasePDE::AnalysisType>& analysis,
                                   std::map<unsigned, unsigned>& numSteps,
                                   bool isHistory = false ) override;

    //! Obtain list with result types in each sequence step
    void GetResultTypes( unsigned sequenceStep, 
                         StdVector<shared_ptr<ResultInfo>>& infos,
                         bool isHistory = false ) override;

    //! Return list with time / frequency values and step for a given result
    void GetStepValues( unsigned sequenceStep,
                                shared_ptr<ResultInfo> info,
                                std::map<unsigned, double>& steps,
                                bool isHistory = false ) override;

    //! Return entitylist the result is defined on
    void GetResultEntities( unsigned sequenceStep,
                             shared_ptr<ResultInfo> info,
                            StdVector<shared_ptr<EntityList> >& list,
                            bool isHistory = false ) override;

    //! Fill pre-initialized results object with values of specified step
    void GetResult( unsigned sequenceStep,
                    unsigned stepNum,
                    shared_ptr<BaseResult> result,
                    bool isHist = false) override;

    //! Read one the strings in the user data group.
    void ReadStringFromUserData(const std::string& dSetName,
                                std::string& str);
    
    //! Query, if hdf5 file has database
    void DB_Init();
    
    //! Obtain content of XML file
    void DB_GetParamFileContent( std::string& params );
    
    //! Obtain content of material file
    void DB_GetMatFileContent( std::string& params );
    
    //! Get coefficients for given vector
    void DB_GetFeFctCoefs( unsigned sequenceStep, unsigned stepNum,
                        const std::string& pdeName,
                        const std::string& funcName,
                        SingleVector * coefs ); 

    //! Return multisequence steps and their analysistypes
    void DB_GetNumMultiSequenceSteps( std::map<unsigned, BasePDE::AnalysisType>& analysis,
                                      std::map<unsigned, double>& accTime,
                                      std::map<unsigned, bool>& isFinished );
                                    
    //! Return PDE and CoefFunctions in given multisequence step
    void DB_GetAvailPdeCoefFcts( unsigned msStep, 
                                 std::map<std::string, 
                                 std::set<std::string> >& coefFcts );
    
    //! Return list with time / frequency values and step for a given result
    void DB_GetStepValues( unsigned sequenceStep,
                           const std::string& pdeName,
                           const std::string& resultName,
                           std::map<unsigned, double>& stepValues );

    //! Close database group
    void DB_Close();

    void GetNamedNodeResult(const std::string& nodeName,
                            const std::string resultName,
                            StdVector<Complex>& result);    

    //! Initialize module with pointer to grid
    void InitModule() override;

    //! Return file name including path
    std::string GetFileName() { return fileName_; }

  private:

    //! Trigger reading of the mesh
    void ReadMesh(Grid *mi) override;

    //! Read nodal and element definitions
    void ReadNodeElemData(hid_t meshGroup);

    //! Read node groups
    void ReadNodeGroups(hid_t meshGroup);

    //! Read element groups
    void ReadElemGroups(hid_t meshGroup);
  
    //! Read meta information about grid
    void ReadMeshStats(hid_t meshGroup);
    
    //! Read mesh result
    void GetMeshResult( unsigned sequenceStep, unsigned stepNum,
                           shared_ptr<BaseResult> result );
    
    //! Read history result
    void GetHistResult( unsigned sequenceStep, unsigned stepNum,
                           shared_ptr<BaseResult> result);
    
    void LinearizeElems(const StdVector<unsigned>& readElems,
                        StdVector<Integer>& elemTypes, 
                        StdVector<unsigned>& globConnect, 
                        StdVector<unsigned>& readNodes);

    void TransformNodes(CoordSystem& coordSys, double scaleFac);

    //! Open entity subgroup for history results, handling mesh<->file number remapping
    hid_t OpenHistEntityGroup(hid_t entityGroup, const std::string& entityTypeString,
                              const EntityIterator& it);

    //! Main hdf5 file
    hid_t mainFile_ = -1;

    //! Root group of main file
    hid_t mainRoot_ = -1;
    
    //! Root for database
    hid_t dbRoot_ = -1;
    
    //! Flag indicating if mesh meta data is already read in
    bool statsRead_ = false;
    
    //! Flag indicating use of external files for mesh results
    bool hasExternalFiles_ = false;

    //! Flag for creating named nodes for each region
    bool genRegionNodes_ = false;
    
    //! Flag indicating if the complete grid is to be loaded
    bool readAllEntities_ = true;
    
    //! Native directory path to hdf5 file
    std::string baseDir_;

    //! List of entities (regions & named nodes/elems) to be read in from the file
    std::set< std::string > readEntities_;

    //! List of entities (regions & named elems) to be linearized (quadr -> lin)
    std::set< std::string > linearizeEntities_;

    //! List with names of regions
    StdVector< std::string > regionNames_;

    //! List with names of node groups
    StdVector< std::string > nodeNames_;

    //! List with names of element groups
    StdVector< std::string > elemNames_;

    // Number of nodes in mesh file
    unsigned numNodes_ = 0;
    
    // Node coordinates
    StdVector<double> nodeCoords_;
    
    // Map from mesh file node numbers to grid node numbers
    std::unordered_map<unsigned, unsigned> f2GNodeNumMap_;

    // Map from grid node numbers to mesh file numbers
    std::unordered_map<unsigned, unsigned> g2FNodeNumMap_;
    
    // Map from mesh file elem numbers to grid elem numbers
    std::unordered_map<unsigned, unsigned> f2GElemNumMap_;

    // Map from grid elem numbers to mesh file elem numbers
    std::unordered_map<unsigned, unsigned> g2FElemNumMap_;
    
    // Map from grid entity nodes to indices of mesh entity nodes
    std::unordered_map<std::string, StdVector<unsigned> > entityNodeMap_;    

    // Coordinate system into which the node coordinates should be mapped
    std::string coordSysId_ = "default";
    
    // Scale factor for node coordinates
    double scaleFac_ = 1.0;
  };

} // end of namespace

#endif
