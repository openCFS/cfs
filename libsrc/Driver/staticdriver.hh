#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

//! driver for static problems. it is derived from BaseDriver
  class StaticDriver : public BaseDriver
{
public:
  //! constructor
  /*!
    \param adomain pointer to class Domain
  */
  StaticDriver(Domain * adomain);

   //! deconstructor 
  virtual ~StaticDriver();
  
  //!  main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem();

protected:

private:

};

}

#endif // FILE_STATICDRIVER
