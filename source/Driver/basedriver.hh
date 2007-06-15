// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "General/environment.hh"

namespace CoupledField
{

  // forward class declarations
  class Domain;
  class WriteResults;

  //! Base class for driving classes where we implemented time-stepping
  class BaseDriver
  {
  public:

    //! Constructor
    BaseDriver();


    //! Destructor
    virtual ~BaseDriver();

    /** Main method for solving the driver problem, in optimization terms, the
     * state problem. The childs MUST either implement SolveProblem() or 
     * SolveProblem(int), there is just no compiler constraint possible! */        
    virtual void SolveProblem();
    
    //! Initialization method
    virtual void Init() = 0;

    /** In optimization case the problem gets the iteration counter for proper
     * writing the results. Repeatedly called in the optimization case.
     * @see SingleDriver */
    virtual void SolveProblem(int optimizationIteration);

    /** <p>For doing optimization we might solve state problems but not to write
     * the results always. Therefore some drivers will move the result writing
     * to an implementation of this method. Drivers not subject to optimization
     * will write the results in SolveProblem() and leave this method blank.</p>
     * <p>SolveProblem() has clearly to be executed in advance.</p> 
     * @param setp_val this names the "transient" step. If -1 this is same ase
     * the internal integer counter. */
    virtual void StoreResults(double step_val = -1.0)
    {
      // intentionally the default behaviour is that SolveProblem() did the job
    };


    //! Return current analysistype

    //! Returns the current analysistype. 
    virtual AnalysisType GetAnalysisType( ) { 
      return analysis_; }
    
    //! Return current (multi)sequenceStep
    virtual UInt GetActSequenceStep();
    
    //! Return current time / frequency step of simulation
    virtual UInt GetActStep ( const std::string& pdename ) = 0;
  
    /** This is an factory pattern implementation. The result is the
     * proper driver based on the analysis type and adaptiviy setting.
     * set this object in the domain and take care for deletion! */
    static BaseDriver* CreateInstance();  
  
    /** We need to differentiate the SingleDrivers from the MultiSequenceDriver */
    typedef enum { SINGLE_DRIVER, MULTI_SEQUENCE_DRIVER } DriverClass;
  
    /** Identification of the driver */
    virtual DriverClass GetDriverClass() = 0;
 
    /** Tells the driver that we run in an optimization context */
    void SetOptimization(bool value) { optimization_ = value; } 
  
  protected:
    
    //! type of analysis
    AnalysisType analysis_;

    //! --------------------- stuff for computation with adaptivity
    //! for printing a sequence of files in dir meshes in gmv-format

    //! counter of meshes for printing meshes
    UInt nummeshes_;  

    //! current analysis step in a multiSequence analysis
    UInt actSequenceStep_;

    //! print mesh in special file. this method is used in adaptive procedure for space.
    void PrintSeqMeshes();

    /** auxiliary function for computation with adaptivity: open files for
     *  printing sequence of refined meshes with error map */ 
    bool printMeshesOrNot();
  
    /** Here we store the last optimationIteration attribute of SolveProblem() to be
     * used in StoreResults() */
    int lastOptimizationIteration_;
    
    /** Tells the driver that we are in an optimization context. 
     * This is similar but still essentially different from a multi sequence in 
     * the sense that the (often static) state problem will be solved multiple times */
    bool optimization_;  
  };

}

#endif // FILE_DRIVER
