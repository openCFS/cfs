// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUTXDMF_HH
#define FILE_CFS_SIMOUTPUTXDMF_HH

#include <set>
#include <map>

#include <Domain/Mesh/Grid.hh>
#include <Domain/Results/ResultInfo.hh>
#include <DataInOut/SimOutput.hh>

#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

namespace CoupledField {

  //! XDMF output writer class
  class SimOutputXDMF: virtual public SimOutput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputXDMF(std::string fileName, PtrParamNode inputNode, 
                  PtrParamNode infoNode, bool isRestart );
    
    //! Destructor
    virtual ~SimOutputXDMF();

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

    //! Set pointer to HDF5 writer, so that  we may control it, when we run as
    //! cfstool. This is necessary since cfstool only supports a single output
    //! writer.
    void SetHDF5Writer(SimOutputHDF5* simOutHDF5, bool inCFSTool) {
      simOutHDF5_ = simOutHDF5;
      inCFSTool_ = inCFSTool;
    }
    
    

  private:

    // =======================================================================
    //  GRID HELPER FUNCTIONS
    // =======================================================================
    //! Write Nodes/Edges/Faces/Elements of Regions to file
    void WriteRegions(std::ofstream& gridFile);

    //! Write list of node groups to file
    void WriteNodeGroups();

    //! Write list of element groups to file
    void WriteElemGroups();

  private:
    
    // =======================================================================
    //  HELPER METHODS
    // =======================================================================
    
    //! Add mesh result
    void AddMeshResult( shared_ptr<BaseResult> sol );
    
    //! Add history result
    void AddHistResult( shared_ptr<BaseResult> sol );
    
    //! Write out XDMF DTD
    void WriteDTD();

    void MapElemType(Elem::FEType eType, UInt& et,
                     UInt& numNodes, bool& wNodeNodes);
      
    // =======================================================================
    //  GENERIC DATA MEMBERS
    // =======================================================================
    
    //@{ \name Generic Data Members

    //! Set with used capabilities, i.e. types of content written to file
    std::set<Capability> usedCapabilities_;
    
    //! Flag indicating if grid is already written
    bool gridWritten_;

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
    
    //! Filename, so it can be used for reopening the file
    std::string currFileName_;

    //! Pointer to HDF5 writer
    SimOutputHDF5* simOutHDF5_;

    //! Are we running as cfstool?
    bool inCFSTool_; 
    
    //@}
    
  };
}
#endif //FILE_CFS_SIMOUTPUTXMDF_HH
