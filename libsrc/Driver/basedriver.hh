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
  virtual void SolveProblem()=0;

  //!
  void SetupMatricesPDE(Integer pdenumber);

  //! write solution in file
//  void WriteResultsInFile(BasePDE * ptPDE, const Integer step, const Double t);

protected:
  //!
  Domain<Point2D> * ptdomain;

private:
  //! options from input-file; if true, then we output first-der,second in output-file
  Boolean SaveDer1,SaveDer2;
  
};

}

#endif // FILE_DRIVER
