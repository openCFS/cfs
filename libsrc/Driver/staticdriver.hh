#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public SingleDriver {

  public:
    //! constructor
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver section
    //! \param true, if driver is part of  multiSequence
    StaticDriver( UInt stepOffset = 0,
                  Double timeOffset = 0.0,
                  std::string driverTag = "anyTag",
                  bool isPartOfSequence = false );

    //! Destructor 
    ~StaticDriver();
  
    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) { return 1;}

    //! Initialization method
    void Init();

    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();

  private:
    
    //! offset for first timestep
    UInt stepOffset_;
    
    //! offset for first time
    Double timeOffset_;
    
  };

}

#endif // FILE_STATICDRIVER
