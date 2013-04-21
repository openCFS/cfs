#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "PDE/BasePDE.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  // forward class declarations
  class Domain;
  class WriteResults;
  class ResultHandler;
  class InfoNode;

  //! Base class for driving classes where we implemented time-stepping
  class BaseDriver
  {
  public:

    //! Constructor
    BaseDriver(shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~BaseDriver();

    //! Initialization method
    virtual void Init() = 0;

    /** This solved the analysis problem and involved generation the global
     * system and solving. 
     * <p>The different analysis types might include multiple harmonic or 
     * transient steps.</p>
     * <p>For optimization this defines a single forward problem, therefore
     * one might skip the writing of the results and call StoreResults()
     * explicitly</p>
     * @param write_results if false nothing is written to the output files
     * @param analysis_id if given is *set* as current and used.
     * @see StoreResults(double) */
    virtual void SolveProblem(bool write_results = true, PtrParamNode analysis_id = PtrParamNode()) = 0;
    
    /** Only of interest for optimization, where one might not want to generate
     * output (gid, hdf5, gmv, ...) for every forward solution. 
     * <p>E.g. because this are linesearch steps of IPOPT. For non-optimization SolveProblem()
     * will have been called with write_results = true.<p>
     * <p>Note that you have to Wrap within a Multisequencestep and finalize the result handler explicitly, 
     * as this can be done only once for HDF5 -> this is done in Optimization::SolveProblem()</p> */
    virtual void StoreResults(UInt stepNum = 1,
                              double step_val = -1.0) { assert(false); }

    //! Return current analysistype

    //! Returns the current analysistype. 
    virtual BasePDE::AnalysisType GetAnalysisType( ) { 
      return analysis_; }
    
    //! Return current (multi)sequenceStep
    virtual UInt GetActSequenceStep();
    
    //! Return current time / frequency step of simulation
    virtual UInt GetActStep ( const std::string& pdename ) = 0;
    
    /** Every analysis step has a unique id stored within an InfoNode.
     * It allows the identification of
     * Analysis information (timestep, frequency, optimization iteration)
     * with OLAS information (time, memory, iterations, residual).
     * @see CreateAnalysisId() and CreateAnalysisIdChild) */
    PtrParamNode GetAnalysisId() { return analysis_id_; }
    
    /** Helper function to create a child analysis step
     * Adds a child-element to base with "analysis_id" = the analysis_id of the base
     * plus ":" plus the child_name (plus ":" plus child_id)
     * @param base where to add to or if NULL then a info/analysis/process/step is created
     * @param child_name e.g. "nonLin", "adjoint" ...
     * @param child_id will be added after child_name. is optional (-1)
     * @return the child element */
    PtrParamNode CreateAnalysisIdChild(PtrParamNode base, const std::string& child_name, int child_id = -1,
        const std::string& child_2_name = "", int child_2_id = -1);

    /** Adds a new analysis id for this driver.
     * @see CreateAnalysisIdChild() */
    PtrParamNode CreateAnalysisId(const std::string& child_name, int child_id = -1,
                               const std::string& child_2_name = "", int child_2_id = -1);

    
    /** This is an factory pattern implementation. The result is the
     * proper driver based on the analysis type and adaptiviy setting.
     * set this object in the domain and take care for deletion! */
    static BaseDriver* CreateInstance(shared_ptr<SimState> state, Domain* domain );  
  
    /** We need to differentiate the SingleDrivers from the MultiSequenceDriver */
    typedef enum { SINGLE_DRIVER, MULTI_SEQUENCE_DRIVER } DriverClass;
  
    /** Identification of the driver */
    virtual DriverClass GetDriverClass() = 0;
    
    /** Is called by optimization to know number of needed result vectors */
    virtual UInt GetNumSteps() { return 1; }
 
    /** Helper method which determines if an AnalyisType is complex. */
    virtual bool IsComplex() = 0;

  protected:
    
    //! type of analysis
    BasePDE::AnalysisType analysis_;

    /** @see GetAnalysisId() */
    PtrParamNode analysis_id_;
    
    /** our report node */ 
    PtrParamNode driverNode; 
    
    //! Pointer to simulation domain
    Domain * domain_;

    //! current analysis step in a multiSequence analysis
    UInt actSequenceStep_;
 
    /** a shortcut to domain->GetResultHandler() */
    ResultHandler* handler_;
    
    //! Simulation state file
    shared_ptr<SimState> simState_;
    
  private:
    /** helper function. Items are separated via ':' to be replaces when using as filename! */
    std::string ConcatAnalysisId(PtrParamNode analysis_id, const std::string& child_name, int child_id, 
                                 const std::string& child_2_name, int child_2_id);
  };

}

#endif // FILE_DRIVER
