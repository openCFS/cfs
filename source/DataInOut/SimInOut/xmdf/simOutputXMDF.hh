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

  //! HDF5 output writer class
  class SimOutputXMDF: virtual public SimOutput {

  public:

    // =======================================================================
    // CONSTRUCTION AND INTIIALIZATION
    // =======================================================================
    //@{ \name Constructor / Initialization
    
    //! Constructor with name of mesh-file
    SimOutputXMDF(std::string fileName, ParamNode * inputNode);
    
    //! Destructor
    virtual ~SimOutputXMDF();

    //! Initialize class 
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly );
    //@}


    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         AnalysisType type,
                                         UInt numSteps  );
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd );

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

    //! Write list of named nodes to file
    void WriteNamedNodes(const H5::Group& meshGroup);

    //! Write list of named elements to file
    void WriteNamedElems(const H5::Group& meshGroup);

    //! Write Meta-Information about results to file
    void WriteResultDescriptions( const H5::Group& descGroup );
    
    //! Write single result to file
    void WriteResults( H5::Group& resultGroup,
                       Vector<Double>& resultVals,
                       const UInt numDOFs,
                       const bool isImag );

    //! Create separate external file for each time / frequency step
    void CreateExternalFile();

    void WriteStringToUserData(const std::string& dSetName, 
                               const std::string& str);

    //! Check for open objects in the hdf5 file
    void CheckOpenObjects();

  private:

    // =======================================================================
    //  HDF5 DATA MEMBERS
    // =======================================================================
    
    //@{ \name H5 Data MEMBERS

    //! Main file containing grid and meta information
    H5::H5File mainFile_;

    //! In case we use
    H5::H5File currStepFile_;

    //! Main / Root Group
    H5::Group mainGroup_;

    //! Mesh Group
    H5::Group meshGroup_;

    //! Group for data / results
    H5::Group dataGroup_;

    //! Group for grid / volume data
    H5::Group volDataGroup_;

    //! Group for current multisequence ste
    H5::Group currMSGroup_;

    //! Group for current analysis step
    H5::Group currStepGroup_;
    //@}

    // =======================================================================
    //  GENERIC DATA MEMBERS
    // =======================================================================
    
    //@{ \name Generic Data Members

    //! Flag indicating if grid is already written
    bool gridWritten_;

    //! Flag indicating if external file per analysis step is used
    bool externalFiles_;

    //! Flag indicating if only grid is to be printed
    bool printGridOnly_;
    
    //! Current multisequence number
    UInt currMS_;
    
    //! Current analysis step
    UInt currStep_;

    //! Type of current analysis type
    AnalysisType currAnalysisType_;

    //! Type definition for registered results
    typedef std::map< std::string, std::vector< 
      boost::shared_ptr<BaseResult> > > ResDescType;
    ResDescType registeredResults_;
    //@}
    
  };
}
#endif //FILE_CFS_SIMOUTPUTXMDF_HH
