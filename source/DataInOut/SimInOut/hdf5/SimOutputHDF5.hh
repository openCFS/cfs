// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUTHDF5_HH
#define FILE_CFS_SIMOUTPUTHDF5_HH

#include <set>
#include <map>

#include <any>

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;
#include "H5Cpp.h"

namespace CoupledField {

  //! Forward class declaration
  class SimState;

  //! define CFS-HDF5 file format version
#define CFS_HDF5_FORMAT_MAJOR 0
#define CFS_HDF5_FORMAT_MINOR 9

  //! HDF5 output writer class
  class SimOutputHDF5: virtual public SimOutput {

    // declare SimState class as friend, as it
    // needs direct write access to the underlying H5 file
    friend class SimState;
    
  public:

    // =======================================================================
    //  CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputHDF5(std::string fileName, PtrParamNode inputNode,
                  PtrParamNode infoNode, bool isRestart );
    
    //! Destructor
    virtual ~SimOutputHDF5();

    //! Initialize class 
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly );
    
    //! Return file name including path
    fs::path GetFileName() {
      return  currFileName_;
    }
    //@}

    // =======================================================================
    //  RESULTS RELATED SECTION
    // =======================================================================
    //@{ \name Result Handling

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         BasePDE::AnalysisType type,
                                         UInt numSteps  );
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd,
                                 bool isHistory );

    //! Begin single analysis step
    virtual void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    virtual void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    virtual void FinishStep( );

    //! End multisequence step
    virtual void FinishMultiSequenceStep( );

    //! Finalize the output
    virtual void Finalize();

    //! Initialize 
    virtual void InitModule();

    //! Write grid
    virtual void WriteGrid();
    //@}
    
    // =======================================================================
    //  DATABASE SECTION
    // =======================================================================
    //@{ \name Database Handling

    //! Initialize database
    void DB_Init();

    //! Write parameter and xml file
    void DB_WriteXmlFiles( fs::path simFile, fs::path matFile );

    //! Write python file
    void DB_WritePythonFile( fs::path pythonFile);

    //! Begin new multisequence step for database section
    void DB_BeginMultiSequenceStep( UInt step,
                                    BasePDE::AnalysisType type );

    //! Begin single analysis step
    void DB_BeginStep( UInt stepNum, Double stepVal );

    //! Write coefficients of coefficient function
    void DB_WriteFeFunction( const std::string& pdeCplName,
                             const std::string& fctName,
                             SingleVector* coefs );

    //! End multisequence step for database section
    void DB_FinishMultiSequenceStep(bool completed, Double accTime );

    //@}

    void WriteStringToUserData(const std::string& dSetName, 
                               const std::string& str);
    
  private:

    // =======================================================================
    //  GRID HELPER FUNCTIONS
    // =======================================================================
    //! Write Nodes/Edges/Faces/Elements of Regions to file
    void WriteRegions(const H5::Group& meshGroup);

    //! Write list of node groups to file
    void WriteNodeGroups(const H5::Group& meshGroup);

    //! Write list of element groups to file
    void WriteElemGroups(const H5::Group& meshGroup);

    //! Write Meta-Information about results to file
    void WriteResultDescriptions( const H5::Group& descGroup,
                                  UInt numSteps,
                                  bool isHistory );

    //! Create separate external file for each time / frequency step
    void CreateExternalFile();

    //! Write file meta information
    void WriteFileInfoHeader(); 
    
    //! Create/Open the file
    void OpenFile(bool truncate);
    
    //! Close the file
    void CloseFile();
    

  private:
    
    // =======================================================================
    //  HELPER METHODS
    // =======================================================================
    
    //! Add mesh result
    void AddMeshResult( shared_ptr<BaseResult> sol );
    
    //! Add history result
    void AddHistResult( shared_ptr<BaseResult> sol );
    
    //! Write single result to file
    void WriteResults( H5::Group& resultGroup,
                       Vector<Double>& resultVals,
                       const UInt numDOFs,
                       const bool isImag );
    
    //! Get lock during low-level file operation
    void LockFile();
    
    //! Release lock after performing file operation and flush file
    void UnlockFile();

    // =======================================================================
    //  HDF5 DATA MEMBERS
    // =======================================================================
    
    //@{ \name H5 Data MEMBERS

    //! Dataset property list (used for chunking and compression )
    H5::DSetCreatPropList dPropList_;
    
    //! Main file containing grid and meta information
    H5::H5File mainFile_;

    //! In case we use
    H5::H5File currStepFile_;

    //! Main / Root Group
    H5::Group mainGroup_;

    //! Mesh Group
    H5::Group meshGroup_;

    //! Group for results
    H5::Group resultsGroup_;

    //! Group for mesh results
    H5::Group meshResultsGroup_;
    
    //! Group for history results
    H5::Group histResultsGroup_;

    //! Group for current multisequence step for mesh results
    H5::Group currMSMeshGroup_;
    
    //! Group for current multisequence step for history results
    H5::Group currMSHistGroup_;

    //! Group for current analysis step for mesh results
    H5::Group currMeshStepGroup_;
    
    //! Group for current analysis step for history results
    H5::Group currHistStepGroup_;
    
    //! Group for internal database
    H5::Group dbGroup_;
    
    //! Group for database entries of current sequence step
    H5::Group currMsDbGroup_;
    //@}

    // =======================================================================
    //  GENERIC DATA MEMBERS
    // =======================================================================
    
    //@{ \name Generic Data Members

    //! Set with used capabilities, i.e. types of content written to file
    std::set<Capability> usedCapabilities_;
    
    //! Flag if module is initialized
    bool isInitialized_;
    
    //! Flag indicating if grid is already written
    bool gridWritten_;

    //! Flag indicating if external file per analysis step is used
    bool externalFiles_;

    //! Flag indicating if only grid is to be printed
    bool printGridOnly_;
    
    //! Flag, if database capability is used
    bool useDataBase_;
    
    //! Flag, if mesh / history results capability is used
    bool useResults_;
    
    //! Current multisequence number
    UInt currMS_;
    
    //! Number of steps in this multisequence step
    UInt currMSNumSteps_;
    
    //! Current analysis step
    UInt currStep_;
    
    //! Current analysis step value (time / frequency);
    Double currStepValue_;

    //! Type of current analysis type
    BasePDE::AnalysisType currAnalysisType_;

    //! Type definition for registered results
    typedef std::map< std::string, std::vector< 
      shared_ptr<BaseResult> > > ResDescType;
    
    //! Set of registered mesh results for current MSstep
    ResDescType registeredMeshResults_;
    
    //! Set of registered history results for current MSstep
    ResDescType registeredHistResults_;
    
    //! Map attribute saveBegin for each result type (mesh)
    std::map<std::string, UInt> meshResultSaveBegin_;
    
    //! Map attribute saveEnd for each result type (mesh)
    std::map<std::string, UInt> meshResultSaveEnd_;
        
    //! Map attribute saveBegin for each result type (mesh)
    std::map<std::string, UInt> meshResultSaveInc_;
    
    //! Map attribute saveBegin for each result type (history)
    std::map<std::string, UInt> histResultSaveBegin_;
        
    //! Map attribute saveEnd for each result type (history)
    std::map<std::string, UInt> histResultSaveEnd_;
            
    //! Map attribute saveBegin for each result type (history)
    std::map<std::string, UInt> histResultSaveInc_;
    
    //! Map with step numbers for each mesh result
    std::map<std::string, StdVector<UInt> > meshResultStepNums_;
    
    //! Map with step numbers for each hist result
    std::map<std::string, StdVector<UInt> > histResultStepNums_;
    
    //! Map with step values for each mesh result
    std::map<std::string, StdVector<Double> > meshResultStepVal_;
    
    //! Map with step values for each hist result
    std::map<std::string, StdVector<Double> > histResultStepVal_;
    
    //! Filename, so it can be used for reopening the file
    std::string currFileName_;
    
    //@}
    
    // =======================================================================
    //  DATABASE SECTION
    // =======================================================================
    
    //! Current analysis step
    UInt currStepDb_;

    //! Current analysis step value (time / frequency);
    Double currStepValueDb_;

    
  };
}
#endif //FILE_CFS_SIMOUTPUTXMDF_HH
