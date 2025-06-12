// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_TEXT_SIM_OUTPUT_HH
#define FILE_CFS_TEXT_SIM_OUTPUT_HH

#include "Domain/Mesh/Grid.hh"
#include "DataInOut/SimOutput.hh"
#include "Domain/ElemMapping/EntityLists.hh"

namespace CoupledField {


  // forward class declaration
  class CoordSystem;

  class SimOutputText : public SimOutput {

  public:

    //! Enum defining type of collection
    typedef enum { ENTITY, TIMEFREQ, ALTOGETHER} CollectionType;

    //! Constructor
    SimOutputText( const std::string& fileName, PtrParamNode outputNode,
                   PtrParamNode infoNode, bool isRestart  );

    //! Destructor
    ~SimOutputText();

    //! Initialize class
    void Init( Grid * ptGrid, bool printGridOnly );
    
    //! Register result (within one multisequence step)
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd,
                         bool isHistory );

    //! Begin multisequence step
    void BeginMultiSequenceStep( UInt step,
                                 BasePDE::AnalysisType type,
                                 UInt numSteps);
    
    //! Begin single analysis step
    void BeginStep( UInt stepNum, Double stepVal );

    //! Add result to current step
    void AddResult( shared_ptr<BaseResult> sol );

    //! End single analysis step
    void FinishStep( );
    
  private:

    //! Create new output stream
    void CreateFiles( shared_ptr<BaseResult> res,
                      UInt step,
                      Double stepVal );
    
    //! Write step for collectionType ENTITY
    void WriteStepCollectEntity();
    
    //! Write step for collecitonType TIMEFREQ
    void WriteStepCollectTimeFreq();
    
    //! Write step for collectionType ALTOGETHER
    void WriteStepCollectAltogether();

    //! Create header information depending on result
    std::string ResultDofString( shared_ptr<BaseResult> res );

    //! Create header information for CSV format depending on result
    std::string ResultDofStringCSV( shared_ptr<BaseResult> res,
                                    const std::string& entityNum,
                                    const std::string& entityType,
                                    const std::string& entityName );
    
    // =======================================================================
    //  Helper functions to iterate over a list of entities
    // =======================================================================    

    //! Return for an entityiterator the node number and coordinates
    StdVector<std::string> GetNodeInfo( EntityIterator & it ) const;
    
    //! Return for an entityiterator the node number and coordinates
    StdVector<std::string> GetElemInfo( EntityIterator & it ) const;
    
    //! Return for an entityiterator the node number
    StdVector<std::string> GetRegionInfo( EntityIterator & it ) const;

    // =======================================================================
    //  Helper function to extract complex values in correct format
    // =======================================================================

    //! Extract complex value as Amplitude-Phase
    std::string ComplexAsAmplPhase( const Complex& val ) const;

    //! Extract complex value as Real-Imag
    std::string ComplexAsRealImag( const Complex& val ) const;

    //! type of collection
    CollectionType collecType_;

    //! Map with result objects for each result type
    ResultMapType resultMap_;
    
    //! Map of result and vector of ofstreams
    std::map<shared_ptr<BaseResult>,
             StdVector<std::ofstream*> > outFiles_;
    
    //! Store all fileNames, which are created during this simulation run
    std::set<std::string> usedFileNames_;

    //! Comment char
    char cmChar_;

    //! Delimiter string
    std::string delim_;

    //! Format output as CSV
    bool csv_;

    //! Coordinate system
    CoordSystem* coordSys_ = nullptr;

    //! Current multi-sequence step
    UInt currMS_ = 0;
    
    //! Type of analysis in current multisequence step
    BasePDE::AnalysisType actAnalysis_ = BasePDE::NO_ANALYSIS;

    // from revision 12363
    // Offset for step number in case of multisequence analysis
    Integer stepNumOffset_;
    
    //! Flag, if node / element numbers in file names should be numbered
    //! based on global node / element number or consecutively
    bool globalNumbering_;
  };
}



#endif
