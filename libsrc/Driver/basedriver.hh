#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "domain.hh"

namespace CoupledField
{

/// there is a base class for driving classes where we implemented time-stepping

template<class Dim>
class BaseDriver
{
public:
  //!
  BaseDriver(Domain<Dim> * adomain);

   //!
  virtual ~BaseDriver();
  
  //! there is a main method, where time-steping is implemented
  virtual void SolveProblem()=0;

  //! here adapt. time-stepping is implemented
  virtual void SolveProblemAdapt() 
  { Error("Not implemented"); }

  //!
  void SetupMatricesPDE(Integer pdenumber);

protected:
  //!
  Domain<Dim> * ptdomain_;

private:
  //! options from input-file; if true, then we output first-der,second in output-file
//  Boolean SaveDer1_, SaveDer2_;
  
};

#ifdef __GNUC__
template class BaseDriver<Point2D>;
template class BaseDriver<Point3D>;
#endif

}

#endif // FILE_DRIVER
