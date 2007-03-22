// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_STATICADAPTIVESPACEDRIVER_2003
#define FILE_STATICADAPTIVESPACEDRIVER_2003

#include "staticdriver.hh"

namespace CoupledField
{

  //! driver for static problems. it is derived from BaseDriver
  class StaticAdaptSpaceDriver : public StaticDriver
  {
  public:
    //! constructor
    /*!
      \param adomain pointer to class Domain
    */
    StaticAdaptSpaceDriver(Domain * adomain);

    //! deconstructor 
    virtual ~StaticAdaptSpaceDriver();
  
    //!  main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem();

  protected:

  private:

  };

}

#endif
