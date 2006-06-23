#ifndef FILE_SINGLEDRIVER_2004
#define FILE_SINGLEDRIVER_2004

#include "basedriver.hh"

namespace CoupledField {

  // forward class declaration
  class BasePDE;

  //! Abstract base class for sinlge driver (static, transient, harmonic)
  class SingleDriver : public BaseDriver {
    
  public:
    
    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    SingleDriver(Domain * adomain, 
                 UInt stepOffset = 0, 
                 Double timeOffset = 0.0, 
                 std::string driverTag ="anyTag",
                 bool isPartOfSequence = false);
    
    //! Default destructor
    virtual ~SingleDriver();
    
    //! set the pde, which gets to be solved
    void SetPDE(BasePDE * pde);
    
    //! main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem()=0;
  
    //! return the actual frequency within a harmonic analysis
    virtual Double GetActFrequency() {return actFreq_; };

  protected:
  
    //! intialize all PDEs

    //! intialize all PDEs: A more detailed description what and how the PDEs
    //! are initialized would be much appreciated!
    void GetMyPDEs();

    //! pointer to basePDE 
    BasePDE * ptPDE_;

    //! true, if driver is part of a multiSequence
    bool isPartOfSequence_;

    //! tag for driver section
    std::string driverTag_;

    //! offset for first timestep
    UInt stepOffset_;

    //! offset for first time
    Double timeOffset_;

    //!
    Double actFreq_;

  };

}

#endif 
