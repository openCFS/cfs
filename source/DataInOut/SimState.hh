// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIM_STATE_HH
#define FILE_CFS_SIM_STATE_HH

#include <map>
#include <boost/filesystem/path.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/signals2.hpp>

#include "General/Environment.hh"
#include "PDE/BasePDE.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

  //! Forward class declarations
  class BaseFeFunction;
  class Domain;  
  class SimOutputHDF5;
  class SimInputHDF5;
  class MathParser;

  //! Class for data handling of internal simulation state
  
  /*! 
   * This class implements the functionality for saving and restoring
   * the internal state of a simulation at any given time / frequency step.
   * It internally utilizes the HDF5 classes to allow for an efficient and
   * transparent way of data handling.
   * 
   * The data necessary to store / load a simulation is:
   *  - Content of XML Parameter File
   *  - Content of Material File
   *  - Coefficient vectors of primary unknowns and their time-derivatives
   *    (e.g. mechanical displacement / velocity / acceleration).
   * 
   * If this information is available, a Domain object including all PDEs
   * and coupling objects can be re-created. But instead of solving for the
   * unknowns in the algebraic system, we use the coefficients from the
   * HDF5 file. 
   */
  class SimState  : public boost::enable_shared_from_this<SimState> {
  
  public:
    
    //! Define available interpolation methods
    //@{
    typedef enum { 
      NO_INTERPOLATION = 0, 
      CONSTANT = 1,
      NEAREST_NEIGHBOR = 2, 
      LINEAR } InterpolType;
    static Enum<InterpolType> InterpolTypeEnum;
    //@}

    //! Constructor
    SimState(bool useAsInput, Domain * parentDomain = NULL );

    //! Destructor
    virtual ~SimState();

    //! Query, if SimState provides input data
    bool HasInput() {
      return hasInput_;
    }
    
    // =======================================
    //  OUTPUT - FUNCTIONALITY 
    // =======================================

    //! Set pointer to hdf5 writer for saving internal data.
    //! If not set, a new writer is created internally
    void SetOutputHdf5Writer( shared_ptr<SimOutputHDF5> outFile );

    shared_ptr<SimOutputHDF5> GetOutputWriter() {
      return outFile_;
    }
    
    //! Pass parameter and material file names
    void SetMatParamFile( const fs::path& paramFile, 
                          const fs::path& matFile ); 

    //! Register new FeFunction to be saved
    void RegisterFeFct( shared_ptr<BaseFeFunction> feFct);
    
    //! Begin new multisequence step
    void BeginMultiSequenceStep( UInt step, BasePDE::AnalysisType type);

    //! Begin single analysis step
    void WriteStep( UInt stepNum, Double stepVal );
    
    //! Finish multisequence step
    void FinishMultiSequenceStep( bool completed, 
                                  Double accumulatedTime = 0.0 );
    
    
    // =======================================
    //  INPUT - FUNCTIONALITY 
    // =======================================
    typedef std::map<std::string, Grid* > GridMap;
    
    //! Retrieve a domain object from a previously set HDF5 file
    Domain* GetDomain(UInt sequenceStep,
                      const GridMap& map = GridMap() );

    //! Set pointer to input hdf5 reader
    void SetInputHdf5Reader( shared_ptr<SimInputHDF5> inFile );
    
    //! Generate input reader from output
    
    //! This method generates an input reader from an already set
    //! input reader.
    bool SetInputReaderToSameOutput();
    
    //! Return list of sequence steps 
    void GetSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                           std::map<UInt, Double>& accTime,
                           std::map<UInt, bool>& isFinished );
    
    //! Return true, if the given sequence step is already completed
    //! in the input file
    bool IsCompleted( UInt sequenceStep );
    
    //! Return number of available sequence steps in file
    UInt GetLastMsStepNum();
    
    //! Return last step number of input file
    void GetLastStepNum(UInt sequenceStep, UInt& lastStepNum,
                        Double& lastStepVal );
 
    //! Set interpolation strategy
    void SetInterpolation( InterpolType type, MathParser * parser,
                           BasePDE::AnalysisType analysis,
                           Double offset );
    
    //! Update internally to given time / frequency step
    void UpdateToStep( UInt sequenceStep, UInt stepNum );
   
    UInt GetSimStateSequenceStep(){return this->sequenceStep_;}
    
    //! Finalize state (delete registered fefunctions etc.)
    void Finalize();
    
  protected:
    
    //! Initialize simulation state object
    void Init();
    
    //! Callback method in case of interpolation of time 
    void UpdateTimeFreqStep();
    
    //! Pointer to HDF5 output class
    shared_ptr<SimOutputHDF5> outFile_;
    
    //! Pointer to HDF5 input class
    shared_ptr<SimInputHDF5> inFile_;
    
    //! Pointer to own domain
    Domain* domain_;
    
    //! Flag, if own writer is created
    bool ownWriter_;
    
    //! Path to parameter file
    fs::path paramFile_;
    
    //! Path to material file
    fs::path matFile_;

    //! Current sequencestep
    UInt sequenceStep_;
    
    //! Flag, if output is initialized yet
    bool isInitialized_;
    
    //! Flag, if simState is used for input
    bool hasInput_;
    
    //! Set of registered FeFunctions
    std::set<shared_ptr<BaseFeFunction> >feFcts_;

    //! Pointer to parent domain
    Domain * parentDomain_;
    // ===========================================
    //  Interpolation Related Data
    // ===========================================
    //! Connection to math parser instance 
    boost::signals2::connection conn_;
    
    //! Math parser for parent domain
    
    //! The math parser of the parent / main domain is neeeded to register
    //! a callback function, which is called everytime the time / frequency
    //! value of the main domain is changed.
    MathParser* parentParser_ = nullptr;
    
    //! Math parser handle for parser of parent domain
    unsigned int parentHandle_;
    
    //! Store stepNumbers
    StdVector<UInt> stepNums_;
    
    //! Store stepValues
    StdVector<Double> stepVals_;
    
    //! store type of interpolation
    InterpolType interpol_;
  };

}

#endif
