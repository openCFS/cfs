// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_OLASPARAMS_HH
#define OLAS_OLASPARAMS_HH

/**********************************************************/

#include "OLAS/algsys/olascomm.hh"

namespace CoupledField {

  /**********************************************************/

  //! Class for objects that feed steering parameters to OLAS

  //! This class implements a communication object that can be used to feed
  //! steering parameters to the OLAS solution process.
  class OLAS_Params : public OLAS_BaseComm {

  public:

    //! default constructor (calls SetDefaultParams())
    OLAS_Params();
    //! Set mandatory default values
    void SetDefaultParams();

  private:

    //! Let the compiler instantiate interfaces for enumeration Pool
    void EnumInterfaces();

  };

  /**********************************************************/
} // namespace CoupledField

#endif // OLAS_OLASPARAMS_HH
