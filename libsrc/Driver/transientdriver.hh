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

  //!
  void SetupMatricesPDE(Integer pdenumber, Integer matrixtype);

protected:


private:
  //!
  Integer numstep,isavebegin,isaveincr,isaveend;

  //!
  Double firstdt;

};

}

#endif // FILE_STATICDRIVER
