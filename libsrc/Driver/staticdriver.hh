#ifndef FILE_STATICDRIVER_2001
#define FILE_STATICDRIVER_2001

#include "basedriver.hh"

namespace CoupledField
{

/// driver for static problems derived from BaseDriver;

template<class Dim>
class StaticDriver : virtual public BaseDriver<Dim>
{
public:
  //!
  StaticDriver(Domain<Dim> * adomain);

   //!
  virtual ~StaticDriver();
  
  //!
  void SolveProblem();

  //!
  void SetupMatricesPDE(const Integer pdenumber, const Integer matrixtype);

protected:

private:

};

#ifdef __GNUC__
template class StaticDriver<Point2D>;
template class StaticDriver<Point3D>;
#endif

}

#endif // FILE_STATICDRIVER
