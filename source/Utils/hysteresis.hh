// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_HYSTERESIS
#define FILE_HYSTERESIS

#include <string>
#include "General/environment.hh"

namespace CoupledField {

  class Hysteresis
  {
  public:
    Hysteresis(Integer numElem);

    ///
    virtual ~Hysteresis();

    //!
    virtual Double computeValue(Double xVal, Integer idxElem) {
      Error( "computeValue not implemented in base-Class",
             __FILE__, __LINE__ );
      return 0.0;
    };

    //! 
    virtual Double getValue(  Integer idxElem) {
      Error( "getValue not implemented in base-Class",
             __FILE__, __LINE__ );
      return 0.0;
    }

    //!
    virtual void updateMinMaxList(Double newX, Integer idxElem) {
      Error( "updateMinMaxList not implemented in base-Class",
             __FILE__, __LINE__ );
    };

    //! 
    virtual void SetTimeStepVal(Double dt) {
      Error( "SetTimeStepVal not implemented in base-Class"
             __FILE__, __LINE__ );
    };


  protected:

  private:

    Integer numElements_;
  };


} //end of namespace


#endif

