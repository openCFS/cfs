#ifndef FILE_TRANSIENTHARMONICDRIVER_2006
#define FILE_TRANSIENTHARMONICDRIVER_2006

#include "basedriver.hh"
#include "singleDriver.hh"

namespace CoupledField {

  //! driver for transient problems.it is derived from BaseDriver;
  class TransientHarmonicDriver : virtual public SingleDriver {

  public:
    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientHarmonicDriver( Domain * adomain,
                             UInt stepOffset = 0,
                             Double timeOffset = 0.0,
                             std::string driverTag = "anyTag",
                             Boolean isPartOfSequence = FALSE );
  
    //! Default destructor
    virtual ~TransientHarmonicDriver();

    //! Main method, where time-stepping is implemented.
    virtual void SolveProblem();

  protected:

  private:
    //! variables for transient driver
    UInt numstep_,isavebegin_,isaveincr_,isaveend_;
    Double firstdt_;

    //! variables for harmonic driver
    Double startFreq_,stopFreq_;
    UInt numFreq_;

    //! pointer to basePDE 
    //BasePDE * ptPDETransient_, ptPDEHarmonic_;


    UInt restartIncr_;



  };

}

#endif // FILE_TRANSIENTHARMONICDRIVER
