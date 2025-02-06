// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_RESULTHANDLER_HH
#define FILE_CFS_RESULTHANDLER_HH

#include "Domain/Results/BaseResults.hh"
#include "PDE/BasePDE.hh"
#include "FeBasis/FeFunctions.hh"
#include <set>
#include <string>
#include <map>

namespace CoupledField {

  //! Forward class declarations
  class SimOutput;
  class SimInput;
  class PostProc;
  class ResultFunctor;

  //! Class for managing several result objects and output classes
  class ResultHandler {

  public:

    //! Constructor
    
    //! General constructor for ResultHandler instance
    //! \param postProcNode Pointer to paramnode for <postProcList> (is allowed
    //!                     to be empty)
    ResultHandler( PtrParamNode rootNode );

    //! Desctructor
    virtual ~ResultHandler();

    // =======================================================================
    // GENERAL METHODS
    // =======================================================================

    //! Add output writer
    void AddOutputDest( shared_ptr<SimOutput> out,
                        const std::string& writerId,
                        const std::string& gridId );
    
    //! Initialize class

    //! This triggers the creation of files and
    //! the writing of the grid information
    void Init(  std::map<std::string, Grid* >& gridMap,
                bool printGridOnly );

    //! Set number of multisequence steps
    void SetNumMultiSequenceSteps( UInt numSteps );
    
    void SetSequenceStep( UInt sequenceStep ) { sequenceStep_ = sequenceStep; };

    // =======================================================================
    // METHODS FOR OUTPUT FUNCTIONALITY
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
                         shared_ptr<ResultFunctor> fnc,
                         UInt sequenceStep,
                         UInt saveBegin, UInt saveInc,
                         UInt saveEnd, 
                         const StdVector<std::string> & outDestNames,
                         const std::string & postProcName,
                         bool writeResult,
                         bool isHistory );


    //! Begin new multisequence step
    void BeginMultiSequenceStep( UInt step, BasePDE::AnalysisType type, UInt numSteps );

    //! Begin new step in analyisis
    void BeginStep( UInt stepNum, Double stepVal );

    //! Update previously registered result
    void UpdateResult( shared_ptr<BaseResult> sol );

    //! New update method, which is just called once
    void UpdateResults();

    //! Finish single analysis step
    void FinishStep( );

    //! Finish multisequence step
    void FinishMultiSequenceStep( );

    //! Last method call, before simulation is finished
    void Finalize( );

    //! Query, if the given result needs to be recomputed in the current step
    bool IsResultNeeded( shared_ptr<BaseResult> res ) 
    { return (isNeeded_.find( res ) != isNeeded_.end()) ; }

    // =======================================================================
    // METHODS FOR INPUT FUNCTIONALITY
    // =======================================================================

    //! Add input writer
    void AddInputReader( shared_ptr<SimInput> inClass, const std::string& readerId );

    //! Retrieve input reader with given id
    shared_ptr<SimInput> GetInputReader( const std::string& readerId );
    
    //! Return number of multisequence steps for a givne inpute read
    void GetNumMultiSequenceSteps( const std::string& readerId,
                                   std::map<UInt, BasePDE::AnalysisType>& analysis,
                                   std::map<UInt, UInt>& numSteps,
                                   bool isHistory = false );

    //! Return result types present in a given input reader
    void GetResultTypes( const std::string& readerId,
                         UInt sequenceStep,
                         StdVector<shared_ptr<ResultInfo> >& infos,
                         bool isHistory = false );
    
    //! Return list with time / frequency values and step for a given result
     void GetStepValues( const std::string& readerId,
                         UInt sequenceStep,
                         shared_ptr<ResultInfo> info,
                         std::map<UInt, Double>& steps,
                         bool isHistory = false );
    
    //! Return the entities for a given result
    void GetResultEntities( const std::string& readerId,
                            UInt sequenceStep,
                            shared_ptr<ResultInfo> info,
                            StdVector<shared_ptr<EntityList> >& list,
                            bool isHistory = false );
    
    //! Fill previously initialized result object
    void GetResult( const std::string& readerId,
                    UInt sequenceStep,
                    UInt stepValue,
                    shared_ptr<BaseResult> result,
                    bool isHistory = false );
    
    //! Commoditiy method for aquiring a result from an input reader class
    shared_ptr<BaseResult> GetResult( const std::string& readerId,
                                      UInt sequenceStep,
                                      UInt stepValue, 
                                      SolutionType solType,
                                      const std::string& regionName );
    
    //! Commoditiy method for aquiring a result from an input reader class
    //! as NodeStoreSol object
    template<typename TYPE>
    shared_ptr<FeFunction<TYPE> > GetFeFunction( const std::string& readerId,
                                                 UInt sequenceStep,
                                                 UInt stepValue,  
                                                 SolutionType solType,
                                                 std::set<std::string> & regionNames,
                                                 PtrParamNode rootNode );

    //! Return the current time/frequency step
    UInt GetCurrStepNum() {return actStep_;}

    //! Return the total number of time/frequency steps
    UInt GetNumSteps() {return numSteps_;}

    /** This dumps the content of the result handler for debugging */
    void Dump(); 

    /** this switch is needed for Optimization to force streaming on every step.
     * If this is set to try, only Outputs with "isStreaming() == true" are being called */
    bool streamOnly = false;


    //! Structure containing additional information about result objects
    struct ResultContext {

      //! Next context
      StdVector<shared_ptr<ResultContext> > nextContexts;

      //! Result type
      shared_ptr<BaseResult> result;
      
      //! ResultFunctor
      shared_ptr<ResultFunctor> functor;

      //! Multisequence step
      UInt sequenceStep;
      
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
      
      //! Flag indicating, if result should be treated as history value
      bool isHistory;

      //! List of postprocessing procedures
      StdVector<shared_ptr<PostProc> > postProcs;

      //! List of outputs, the result gets written to
      StdVector<std::string > outputIds;

    };

    //! Allows access to the ResultContexts
    std::map<shared_ptr<BaseResult>, shared_ptr<ResultContext> >* GetResultContexts() {
      return &resultContexts_;
    };

    
  private:
    
    // =======================================================================
    // METHODS FOR OUTPUT FUNCTIONALITY
    // =======================================================================


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
  
    //!  Finalize all results defined by postProcs
    void FinishMultiSequenceStepRec( ResultContext& context );

    //! Interpolate result to a different grid
    void InterpolateRes( ResultContext& context,
                         Grid* destGrid,
                         shared_ptr<BaseResult>& mappedResult );
    
    //! Pointer to root parameter element
    PtrParamNode param_;
        
    //! Set of ResultContexts
    std::set<shared_ptr<ResultContext> > contexts_;

    //! Map relating result objects to their context
    std::map<shared_ptr<BaseResult>, 
             shared_ptr<ResultContext> > resultContexts_;

    //! Map relating outputIds with the output classes
    std::map<std::string, shared_ptr<SimOutput> > simOutputHandlers_;
    
    //! Map relating outputIds with the gridIds
     std::map<std::string, std::string > outGridIds_;

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
    
    //! Number of steps in current analysis
    UInt numSteps_;
    
    //! Indicates if one result is to be written in final step
    bool finalResultExists_;

    //! Current analysistype
    BasePDE::AnalysisType analysisType_;
    
    // =======================================================================
    // METHODS FOR INPUT FUNCTIONALITY
    // =======================================================================

    //! Map relating inputIds with the simInput classes
    std::map<std::string, shared_ptr<SimInput> > inFiles_;

  };
    
}

#endif
