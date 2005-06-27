#ifndef FILE_TRANSIENT4SLICEDRIVER_2004
#define FILE_TRANSIENT4SLICEDRIVER_2004

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for transient problems usinf slicing technique
  class Transient4SliceDriver : virtual public SingleDriver {

  public:
    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    Transient4SliceDriver( Domain * adomain,
		     Integer stepOffset = 0,
		     Double timeOffset = 0.0,
		     std::string driverTag = "anyTag",
		     Boolean isPartOfSequence = FALSE );
  
    //! Default destructor
    virtual ~Transient4SliceDriver();

    //! main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem();

  protected:

  private:
    //!
    Integer numstep_,isavebegin_,isaveincr_,isaveend_;

    //for structured grid
    Integer elementsPerWavelength_;
    Double dimZ_;
    Double waveLength_;
    Double pulseTime_;
    Double soundSpeed_;
    Integer sik_;
    Integer tol_;
    Domain *domain_;
    //!
    Double firstdt_;
  };

}

#endif // FILE_TRANSIENTDRIVER
