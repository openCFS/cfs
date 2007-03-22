// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
