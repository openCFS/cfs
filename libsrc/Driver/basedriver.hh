#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "domain.hh"
//#include "pde.hh"

namespace CoupledField
{

/// class where we implement Newmark method
class BaseDriver
{
public:
  //!
  BaseDriver(Domain<Point2D> * adomain);

   //!
  virtual ~BaseDriver();
  
  //!
  void SolveProblem();

  //!
  void SetupMatricesPDE(Integer pdenumber);

protected:
  //!
  Domain<Point2D> * ptdomain;

private:

};

}

#endif // FILE_DRIVER
