#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

//! driver for transient problems.it is derived from BaseDriver;
  class TransientDriver : virtual public BaseDriver
{
public:
  //!  constructor
  /*!
    \param adomain pointer to class Domain
  */
  TransientDriver(Domain * adomain);

   //! deconstructor 
  virtual ~TransientDriver();

  //!  main method, where time-stepping is implemented. it is for transient and static problem
  virtual void SolveProblem();

protected:

private:
  //!
  Integer numstep_,isavebegin_,isaveincr_,isaveend_;

  //!
  Double firstdt_;

};

}

#endif // FILE_TRANSIENTDRIVER
