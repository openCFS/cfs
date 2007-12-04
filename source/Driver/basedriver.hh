// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "PDE/basePDE.hh"
#include "General/environment.hh"

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
    BaseDriver();


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
     * @param comment if ony uses export linear system in the xml file,
     *        this becomes part of the base filename. E.g. for naming the
     *        adjoint PDEs in optimization
     * @see StoreResults(double) */
    virtual void SolveProblem(bool write_results = true, const std::string& comment = "") = 0;
    
    /** Only of interest for optimization, where one might not want to generate
     * output (gid, hdf5, gmv, ...) for every forward solution. 
     * <p>E.g. because this are linesearch steps of IPOPT. For non-optimization SolveProblem()
     * will have been called with write_results = true.<p>
     * <p>Note that you have to Wrap within a Multisequencestep and finalize the result handler explicitly, 
     * as this can be done only once for HDF5 -> this is done in Optimization::SolveProblem()</p> */
    virtual void StoreResults(double step_val = -1.0) { assert(false); }

    //! Return current analysistype

    //! Returns the current analysistype. 
    virtual BasePDE::AnalysisType GetAnalysisType( ) { 
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
 
  protected:
    
    //! type of analysis
    BasePDE::AnalysisType analysis_;

    //! --------------------- stuff for computation with adaptivity
    //! for printing a sequence of files in dir meshes in gmv-format

    //! counter of meshes for printing meshes
    UInt nummeshes_;  

    //! current analysis step in a multiSequence analysis
    UInt actSequenceStep_;
 
    /** a shortcut to domain->GetResultHandler() */
    ResultHandler* handler_;
    
    //! print mesh in special file. this method is used in adaptive procedure for space.
    void PrintSeqMeshes();

    /** auxiliary function for computation with adaptivity: open files for
     *  printing sequence of refined meshes with error map */ 
    bool printMeshesOrNot();
    
  };

}

#endif // FILE_DRIVER
