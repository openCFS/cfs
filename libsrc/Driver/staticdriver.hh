#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public SingleDriver {

  public:
    //! constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time
    //! \param driverTag tag for current driver section
    //! \param true, if driver is part of  multiSequence
    StaticDriver(Domain * adomain,
                 UInt stepOffset = 0,
                 Double timeOffset = 0.0,
                 std::string driverTag = "anyTag",
                 Boolean isPartOfSequence = FALSE);

    //! Destructor 
    ~StaticDriver();
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();

  };

}

#endif // FILE_STATICDRIVER
