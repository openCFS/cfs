#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "domain.hh"

namespace CoupledField
{

/// there is a base class for driving classes where we implemented time-stepping

class BaseDriver
{
public:
  //!
  BaseDriver(Domain<Point2D> * adomain);

   //!
  virtual ~BaseDriver();
  
  //! there is a main method, where time-steping is implemented
  virtual void SolveProblem()=0;

  //!
  void SetupMatricesPDE(Integer pdenumber);

protected:
  //!
  Domain<Point2D> * ptdomain_;

private:
  //! options from input-file; if true, then we output first-der,second in output-file
//  Boolean SaveDer1_, SaveDer2_;
  
};

}

#endif // FILE_DRIVER
