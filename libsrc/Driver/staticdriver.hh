#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField {

  //! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public BaseDriver {

  public:
    //! Constructor

    //! \param adomain pointer to class Domain
    StaticDriver(Domain * adomain);

    //! Destructor 
    ~StaticDriver();
  
    //! Main method solution method

    //! This method constitutes the actual driving method which controls the
    //! solution process for the problem.
    void SolveProblem();

  };

}

#endif // FILE_STATICDRIVER
