#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

/// driver for static problems derived from basedriber;

  class TransientDriver : virtual public BaseDriver
{
public:
  //!
  TransientDriver(Domain<Point2D> * adomain);

   //!
  virtual ~TransientDriver();

  //!
  virtual void SolveProblem();

  //! with adaptivity in time
  virtual void SolveProblemAdapt();

  //!
  void SetupMatricesPDE(Integer pdenumber, Integer matrixtype);

protected:

private:
  //!
  Integer numstep_,isavebegin_,isaveincr_,isaveend_;

  //!
  Double firstdt_;

};

}

#endif // FILE_TRANSIENTDRIVER
