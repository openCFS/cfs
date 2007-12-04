// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUTXMDF_HH
#define FILE_CFS_SIMOUTPUTXMDF_HH

#include <set>
#include <map>

#include <boost/any.hpp>

#include <Domain/grid.hh>
#include <Domain/resultInfo.hh>
#include <DataInOut/simOutput.hh>

#include "H5Cpp.h"

namespace CoupledField {

  //! define CFS-HDF5 file format version
#define CFS_HDF5_FORMAT_MAJOR 0
#define CFS_HDF5_FORMAT_MINOR 9


  //! HDF5 output writer class
  class SimOutputHDF5: virtual public SimOutput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputHDF5(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimOutputHDF5();

    //! Initialize class 
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly );
    //@}


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

    void WriteStringToUserData(const std::string& dSetName, 
                               const std::string& str);

    //! Check for open objects in the hdf5 file
    void CheckOpenObjects();
    
    //! Write file meta information
    void WriteFileInfoHeader(); 
    

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
    //@}

    // =======================================================================
    //  GENERIC DATA MEMBERS
    // =======================================================================
    
    //@{ \name Generic Data Members

    //! Set with used capabilities, i.e. types of content written to file
    std::set<Capability> usedCapabilities_;
    
    //! Flag indicating if grid is already written
    bool gridWritten_;

    //! Flag indicating if external file per analysis step is used
    bool externalFiles_;

    //! Flag indicating if only grid is to be printed
    bool printGridOnly_;
    
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
    
    
    
    //@}
    
  };
}
#endif //FILE_CFS_SIMOUTPUTXMDF_HH
