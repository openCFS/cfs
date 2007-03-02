#ifndef FILE_HARMONICADAPTDRIVER_2001
#define FILE_HARMONICADAPTDRIVER_2001

#include "harmonicDriver.hh"

namespace CoupledField
{

  //! driver for static problems. it is derived from BaseDriver
  class HarmonicAdaptSpaceDriver : virtual public HarmonicDriver
  {
  public:
    //! constructor
    /*!
      \param adomain pointer to class Domain
    */
    HarmonicAdaptSpaceDriver(Domain * adomain);

    //! deconstructor 
    virtual ~HarmonicAdaptSpaceDriver();
  
    //!  main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem();

  };

}

#endif // FILE_HARMONICADAPTDRIVER
