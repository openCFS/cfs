#ifndef FILE_TRANSIENTHARMONICDRIVER_2006
#define FILE_TRANSIENTHARMONICDRIVER_2006

#include "transientdriver.hh"
#include "harmonicDriver.hh"


namespace CoupledField {

  //! driver for transient problems.it is derived from BaseDriver;
  class TransientHarmonicDriver : public virtual TransientDriver, 
                                  public virtual HarmonicDriver {

  public:
    //! Constructor
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientHarmonicDriver( UInt stepOffset = 0,
                             Double timeOffset = 0.0,
                             std::string driverTag = "anyTag",
                             bool isPartOfSequence = false );
  
    //! Default destructor
    virtual ~TransientHarmonicDriver();

    //! Initialization method
    void Init();

    //! Main method, where time-stepping is implemented.
    void SolveProblem();

    //! Return current analysistype
    
    //! Returns the current analysistype. 
    //! \param pdeName Name of the pdename in case there is a coexistence
    //!                of two different analysistypes in one analysisstep
    //!                (e.g. transient-harmonic)
    AnalysisType GetAnalysisType( const std::string& pdename );

  //! Return current time / frequency step of simulation
  UInt GetActStep( const std::string& pdename );
  
  protected:

  private:

  };

}

#endif // FILE_TRANSIENTHARMONICDRIVER
