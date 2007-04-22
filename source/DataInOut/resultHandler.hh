// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_RESULTHANDLER_HH
#define FILE_CFS_RESULTHANDLER_HH

#include "Utils/result.hh"

namespace CoupledField {

  //! Forward class declarations
  class SimOutput;
  class PostProc;

  //! Class for managing several result objects and output classes
  class ResultHandler {

  public:

    //! Typedef for operation mode of resultHandler
    //! EMBEDDED: Denotes that the results are generated during a
    //!           simulation run of CFS++
    //! STANDALONE: All results are available at once by reading in
    //!             from a file
    typedef enum {EMBEDDED, STANDALONE} OpMode;

    //! Constructor
    ResultHandler( OpMode mode );

    //! Desctructor
    virtual ~ResultHandler();

    // =======================================================================
    // GENERAL METHODS
    // =======================================================================

    //! Trigger only writing of the grid
    void WriteGrid();

    //! Add output writer
    void AddOutputDest( shared_ptr<SimOutput> out, const std::string& id);
    
    //! Initialize class

    //! This triggers the creation of files and
    //! the writing of the grid information
    void Init( Grid* ptGrid );

    //! Set number of multisequence steps
    void SetNumMultiSequenceSteps( UInt numSteps );
    
    // =======================================================================
    // METHODS FOR EMBEDDED FUNCTIONALITY
    // =======================================================================

    //! Register new result

    //! Register new result object. The result object doest not have to be
    //! filled already with values. Only the contained entitylist and 
    //! resultInfo data has to be set.
    //! Additional parameters are the start, end and increment values 
    //! (in case of transient / harmonic analysis), for which the date should
    //! be written. Also an associated output destination has to be given, where
    //! the result will be written to.
    void RegisterResult( shared_ptr<BaseResult> sol,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd, 
                         const StdVector<std::string> & outDestNames,
                         const std::string & postProcName = "",
                         bool writeResult = true );


    //! Begin new multisequence step
    void BeginMultiSequenceStep( UInt step, AnalysisType type, UInt numSteps );

    //! Begin new step in analyisis
    void BeginStep( UInt stepNum, Double stepVal );

    //! Update previously registered result
    void UpdateResult( shared_ptr<BaseResult> sol );
    
    //! Finish single analysis step
    void FinishStep( );

    //! Finish multisequence step
    void FinishMultiSequenceStep( );

    //! Last method call, before simulation is finished
    void Finalize( );

    //! Query, if the given result needs to be recomputed in the ccurent step
    bool IsResultNeeded( shared_ptr<BaseResult> res ) 
    { return (isNeeded_.find( res ) != isNeeded_.end()) ; }

    // =======================================================================
    // METHODS FOR STANDALONE FUNCTIONALITY
    // =======================================================================

    // ... to be implemented ...

  private:
    
    //! Structure containing additional information about result objects
    struct ResultContext {

      //! Next context
      StdVector<shared_ptr<ResultContext> > nextContexts;

      //! Result type
      shared_ptr<BaseResult> result;

      //! Step begin for saving
      UInt saveBegin;
      
      //! Last step for saving result
      UInt saveEnd;
      
      //! Increment for saving result
      UInt saveInc;

      //! Flag, if result is written to file
      bool writeResult;

      //! Flag indicating, if a result should be written only as final result
      bool isFinal;

      //! List of postprocessing procedures
      StdVector<shared_ptr<PostProc> > postProcs;

      //! List of outputs, the result gets written to
      StdVector<shared_ptr<SimOutput> > outputs;

    };

    //! Checks, if a result has to be written in the current step
    bool IsOutput( ResultContext& context );
    
    //! Return suitable output class(es) for given result object
    void GetDefaultOutputs( shared_ptr<BaseResult> res,
                            StdVector<std::string>& out );

    //! Register results recursively for postprocs
    void RegisterResultRec( ResultContext& actContext, 
                            const std::string& postProcName );

    //! Update results recursivly (postprocessing results)
    void UpdateResultRec( ResultContext& actContext );    

    //! Finish steps recursively, for all results defined by postProcs
    void FinishStepRec( ResultContext& context );
  
    //!  Finalize  all results defined by postProcs
    void FinalizeRec( ResultContext& context );

    //! Operation mode of resultHandler
    OpMode opMode_;
    
    //! Set of ResultContexts
    std::set<shared_ptr<ResultContext> > contexts_;

    //! Map relating result objects to their context
    std::map<shared_ptr<BaseResult>, 
             shared_ptr<ResultContext> > resultContexts_;

    //! Map relating outputIds with the output classes
    std::map<std::string, shared_ptr<SimOutput> > outFiles_;

    //! Set containing all result types needed for the current simulation step
    std::set<shared_ptr<BaseResult> > isNeeded_;

    //! Set containing all results, which are updated in the current step
    std::set<shared_ptr<BaseResult> > isUpdated_;

    //! Current multiSequence step
    UInt sequenceStep_;

    //! Current time step
    UInt actStep_;

    //! Current time / frequency
    Double actStepVal_;
    
    //! Indicates if one result is to be written in final step
    bool finalResultExists_;

    //! Current analysistype
    AnalysisType analysisType_;

  };
    
}

#endif
