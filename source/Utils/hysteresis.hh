// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/defs.hh"
#include "General/exception.hh"

namespace CoupledField {

  class Hysteresis
  {
  public:
    Hysteresis(Integer numElem);

    ///
    virtual ~Hysteresis();

    //!
    virtual Double computeValue(Double xVal, Integer idxElem, bool overwrite = true) {
      EXCEPTION( "computeValue not implemented in base-Class" );
      return 0.0;
    };


    //!
    virtual Double computeValueAndUpdate(Double xVal, Integer idxElem, 
                                         bool overwrite = true ) {
      EXCEPTION( "computeValueAndUpdate not implemented in base-Class" );
      return 0.0;
    };


    //! 
    virtual Double getValue(  Integer idxElem) {
      EXCEPTION( "getValue not implemented in base-Class" );
      return 0.0;
    }

    //! 
    virtual UInt getStringLength( Integer idxElem ) {
      EXCEPTION( "getStringLength not implemented in base-Class" );
      return 0;
    }

    //!
    virtual Double GetIncX() {
      EXCEPTION( "GetIncX not implemented in base-Class" );
      return 1.0;
    };

    //!
    virtual Double updateMinMaxList(Double newX, Integer idxElem, 
                                  bool overwrite) {
      EXCEPTION( "updateMinMaxList not implemented in base-Class" );
      return 0;
    };

    //! 
    virtual void SetTimeStepVal(Double dt) {
      EXCEPTION( "SetTimeStepVal not implemented in base-Class" );
    };

    //!
    virtual Double EvalEverett(Double x1, Double x2, Integer idx ) {
      EXCEPTION( " EvalEverett not implemented in base-Class" );
      return 0.0;
    };

  protected:

  private:

    Integer numElements_;
  };


} //end of namespace


#endif

