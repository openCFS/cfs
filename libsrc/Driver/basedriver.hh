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
  BaseDriver(Domain * adomain);

   //!
  virtual ~BaseDriver();
  
  //! there is a main method, where time-steping is implemented
  virtual void SolveProblem()=0;

  //! here adapt. time-stepping is implemented
  virtual void SolveProblemAdapt() 
  { Error("Not implemented",__FILE__,__LINE__); }

  virtual void SolveProblemAdaptSpace()
  { Error("Not implemented",__FILE__,__LINE__); }

  //!
  void SetupMatricesPDE(Integer pdenumber);

protected:
  //!
  Domain * ptdomain_;

  //! for printing a sequence of files in dir meshes in gmv-format
  WriteResults * ptMeshes_;
  Integer nummeshes_;  // counter of meshes

  //!
  void PrintSeqMeshes();

private:
  //! options from input-file; if true, then we output first-der,second in output-file
//  Boolean SaveDer1_, SaveDer2_;
  
};

}

#endif // FILE_DRIVER
