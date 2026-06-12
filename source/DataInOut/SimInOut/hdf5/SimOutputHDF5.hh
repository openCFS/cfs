// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUTHDF5_HH
#define FILE_CFS_SIMOUTPUTHDF5_HH

#include <set>
#include <map>
#include <chrono>
#include <H5Ipublic.h>

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>
#include <filesystem>

namespace fs = std::filesystem;

namespace CoupledField {

  //! Forward class declaration
  class SimState;

  //! HDF5 output writer class
  class SimOutputHDF5: public SimOutput 
  {
    // declare SimState class as friend, as it
    // needs direct write access to the underlying H5 file
    friend class SimState;
    
  public:

    //! Constructor with name of mesh-file
    SimOutputHDF5(std::string fileName, PtrParamNode inputNode,
                  PtrParamNode infoNode, bool isRestart );
    
    //! Destructor
    virtual ~SimOutputHDF5();

    //! Initialize class 
    void Init( Grid* ptGrid, bool printGridOnly ) override;
    
    //! Return file name including path
    fs::path GetFileName() { return  currFileName_; }
    
    //! Begin multisequence step
    void BeginMultiSequenceStep( unsigned int step,
                                         BasePDE::AnalysisType type,
                                         unsigned int numSteps  ) override;
    
    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                                 unsigned int saveBegin, unsigned int saveInc,
                                 unsigned int saveEnd,
                                 bool isHistory ) override;

    //! Begin single analysis step
    void BeginStep( unsigned int stepNum, double stepVal ) override;

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol ) override;

    //! End single analysis step
    void FinishStep( ) override;

    //! End multisequence step
    void FinishMultiSequenceStep( ) override;

    //! Finalize the output
    void Finalize() override;

    //! Initialize database
    void DB_Init();

    //! Write parameter and xml file
    void DB_WriteXmlFiles( fs::path simFile, fs::path matFile );

    //! Begin new multisequence step for database section
    void DB_BeginMultiSequenceStep( unsigned int step, BasePDE::AnalysisType type );

    //! Begin single analysis step
    void DB_BeginStep( unsigned int stepNum, double stepVal );

    //! Write coefficients of coefficient function
    void DB_WriteFeFunction( const std::string& pdeCplName,
                             const std::string& fctName,
                             SingleVector* coefs );

    //! End multisequence step for database section
    void DB_FinishMultiSequenceStep(bool completed, double accTime );

    void WriteStringToUserData(const std::string& dSetName, const std::string& str);
    
  private:

    //! Initialize 
    void InitModule();

    //! Write grid
    void WriteGrid();

    //! Write Nodes/Edges/Faces/Elements of Regions to file
    void WriteRegions(hid_t meshGroup);

    //! Write list of node groups to file
    void WriteNodeGroups(hid_t meshGroup);

    //! Write list of element groups to file
    void WriteElemGroups(hid_t meshGroup);

    //! Write Meta-Information about results to file
    void WriteResultDescriptions( hid_t descGroup,
                                  unsigned int numSteps,
                                  bool isHistory );

    //! Create separate external file for each time / frequency step
    void CreateExternalFile();

    //! Write file meta information
    void WriteFileInfoHeader(); 
    
    /** Create/Open the file
     * @param truncate for warmstart */
    void OpenFile(bool truncate);
    
    //! Close the file
    void CloseFile();
    
    //! Add mesh result
    void AddMeshResult( shared_ptr<BaseResult> sol );
    
    //! Add history result
    void AddHistResult( shared_ptr<BaseResult> sol );
    
    //! Write single result to file
    void WriteResults( hid_t resultGroup,
                       Vector<double>& resultVals,
                       const unsigned int numDOFs,
                       const bool isImag );
    
    /** test if he shall have an intermediate flush. 
     * CloseFile() implicitly flushed, no need to calls this function.
     * for very slow simulations we flush and then ParaView can check intermediate results */
    void AutoFlush();

    //! Main file containing grid and meta information
    hid_t mainFile_ = -1;

    //! In case we use
    hid_t currStepFile_ = -1;

    //! Main / Root Group
    hid_t mainGroup_ = -1;

    //! Mesh Group
    hid_t meshGroup_ = -1;

    //! Group for results
    hid_t resultsGroup_ = -1;

    //! Group for mesh results
    hid_t meshResultsGroup_ = -1;
    
    //! Group for history results
    hid_t histResultsGroup_ = -1;

    //! Group for current multisequence step for mesh results
    hid_t currMSMeshGroup_ = -1;
    
    //! Group for current multisequence step for history results
    hid_t currMSHistGroup_ = -1;

    //! Group for current analysis step for mesh results
    hid_t currMeshStepGroup_ = -1;
    
    //! Group for current analysis step for history results
    hid_t currHistStepGroup_ = -1;
    
    //! Group for internal database
    hid_t dbGroup_ = -1;
    
    //! Group for database entries of current sequence step
    hid_t currMsDbGroup_ = -1;

    //! Set with used capabilities, i.e. types of content written to file
    std::set<Capability> usedCapabilities_;
    
    //! Flag if module is initialized
    bool isInitialized_ = false;
    
    //! Flag indicating if grid is already written
    bool gridWritten_ = false;

    //! Flag indicating if external file per analysis step is used
    bool externalFiles_ = false;

    //! Flag indicating if only grid is to be printed
    bool printGridOnly_ = false;
    
    //! Flag, if database capability is used
    bool useDataBase_ = false;
    
    //! Flag, if mesh / history results capability is used
    bool useResults_ = false;
    
    /** seconds for auto flush. 0.0 means we flush on every AutoFlush() call, which is not too bad */
    double autoFlushSeconds_ = 2.0;

    /** compression level: 0 disables, 1 is fast and good, 2 ... 9 is usually slow but not much better */
    int compressionLevel_ = 1;

    //! Current multisequence number
    unsigned int currMS_ = 0;
    
    //! Number of steps in this multisequence step
    unsigned int currMSNumSteps_ = 0;
    
    //! Current analysis step
    unsigned int currStep_ = 0;
    
    //! Current analysis step value (time / frequency);
    double currStepValue_ = 0.0;

    //! Type of current analysis type
    BasePDE::AnalysisType currAnalysisType_ = BasePDE::AnalysisType::NO_ANALYSIS;
    
    //! Set of registered mesh results for current MSstep
    std::map<std::string, StdVector<shared_ptr<BaseResult>>> registeredMeshResults_;
    
    //! Set of registered history results for current MSstep
    std::map<std::string, StdVector<shared_ptr<BaseResult>>> registeredHistResults_;
    
    //! Map with step numbers for each mesh result
    std::map<std::string, StdVector<unsigned int> > meshResultStepNums_;
    
    //! Map with step numbers for each hist result
    std::map<std::string, StdVector<unsigned int> > histResultStepNums_;
    
    //! Map with step values for each mesh result
    std::map<std::string, StdVector<double> > meshResultStepVal_;
    
    //! Map with step values for each hist result
    std::map<std::string, StdVector<double> > histResultStepVal_;
    
    //! Filename, so it can be used for reopening the file
    std::string currFileName_;
    
    shared_ptr<Timer> initTimer_;

    /** last action point for automatic flush */
    std::chrono::steady_clock::time_point lastFlush_;
        
    //! Current analysis step
    unsigned int currStepDb_ = 0;

    //! Current analysis step value (time / frequency);
    double currStepValueDb_ = 0.0;
  };
}
#endif // end FILE_CFS_SIMOUTPUTHDF5_HH
