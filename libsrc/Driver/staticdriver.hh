#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

/// driver for static problems derived from BaseDriver;

  class StaticDriver : virtual public BaseDriver
{
public:
  //!
  StaticDriver(Domain<Point2D> * adomain);

   //!
  virtual ~StaticDriver();
  
  //!
  void SolveProblem();

  //!
  void SetupMatricesPDE(Integer pdenumber, Integer matrixtype);

protected:


private:

};

}

#endif // FILE_STATICDRIVER
