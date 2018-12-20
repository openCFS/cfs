#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "PDE/BasePDE.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Driver/AnalysisID.hh"

namespace CoupledField
{

  // forward class declarations
  class Domain;
  class WriteResults;
  class ResultHandler;

  //! Base class for driving classes where we implemented time-stepping
  class BaseDriver
  {
  public:

    //! Constructor
    BaseDriver(shared_ptr<SimState> simState, Domain* domain, 
               PtrParamNode paramNode, PtrParamNode infoNode );

    //! DestructorInverseSourceDriver
    virtual ~BaseDriver();

    //! Initialization method
    
    //! Initialize driver
    //! \param restart If true, the driver performs restart
    virtual void Init(bool restart) = 0;

    /** This solved the analysis problem and involved generation the global
     * system and solving. 
     * <p>The different analysis types might include multiple harmonic or 
     * transient steps.</p>
     * <p>For optimization this defines a single forward problem, therefore
     * one might skip the writing of the results and call StoreResults()
     * explicitly</p>
     * @see StoreResults(double) */
    virtual void SolveProblem() = 0;
    
    /** Only of interest for optimization, where one might not want to generate
     * output (gid, hdf5, gmv, ...) for every forward solution. 
     * <p>E.g. because this are linesearch steps of IPOPT. For non-optimization SolveProblem()
     * will have been called with write_results = true.<p>
     * <p>Note that you have to Wrap within a Multisequencestep and finalize the result handler explicitly, 
     * as this can be done only once for HDF5 -> this is done in Optimization::SolveProblem()</p> */
    virtual void StoreResults(UInt stepNum = 1, double step_val = -1.0) { assert(false); }

    //! Return current analysistype

    //! Returns the current analysistype. 
    virtual BasePDE::AnalysisType GetAnalysisType( ) { 
      return analysis_; }
    
    //! Return current (multi)sequenceStep
    virtual UInt GetActSequenceStep();
    
    //! Return current time / frequency step of simulation
    virtual UInt GetActStep ( const std::string& pdename ) = 0;
    
    /** Give the current analysis_id which describes the current algebraic problem. */
    AnalysisID& GetAnalysisId() { return analysis_id_; }

    /** shortcut the Bloch Mode Eigenfrequency Analysis */
    virtual bool DoBlochModeEigenfrequency() const { return false; }

    /** This is an factory pattern implementation. The result is the
     * proper driver based on the analysis type and adaptiviy setting.
     * set this object in the domain and take care for deletion! */
    static BaseDriver* CreateInstance(shared_ptr<SimState> state, Domain* domain,
                                      PtrParamNode paramNode, PtrParamNode infoNode);  
  
    /** We need to differentiate the SingleDrivers from the MultiSequenceDriver */
    typedef enum { SINGLE_DRIVER = 0, MULTI_SEQUENCE_DRIVER } DriverClass;
  
    /** Identification of the driver */
    virtual DriverClass GetDriverClass() = 0;
    
    /** Is called by optimization to know number of needed result vectors */
    virtual UInt GetNumSteps() { return 1; }
 
    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() = 0;

    PtrParamNode GetInfo() { return info_; }

    PtrParamNode GetParam() { return param_; }

    ResultHandler* GetResultHandler() { return handler_; }
  protected:
    
    //! type of analysis
    BasePDE::AnalysisType analysis_;

    //! Pointer to parameter node
    PtrParamNode param_ ;
    
    //! Pointer to information node (general information)
    PtrParamNode info_;
    
    /** Identifies the current algebraic problem. Used for info.xml and exportLinSys */
    AnalysisID analysis_id_;

    //! Pointer to simulation domain
    Domain* domain_;

    //! current analysis step in a multiSequence analysis
    UInt sequenceStep_;
 
    /** a shortcut to domain->GetResultHandler() */
    ResultHandler* handler_;
    
    //! Simulation state file
    shared_ptr<SimState> simState_;
    
    //! Static flag to HALT the simulation
    static bool abortSimulation_;
  };

}

#endif // FILE_DRIVER
