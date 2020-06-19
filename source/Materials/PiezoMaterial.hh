// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef PIEZOMATERIAL_DATA
#define PIEZOMATERIAL_DATA

#include <string>

#include "General/defs.hh"
#include "General/Environment.hh"
#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling piezo material data
  */

template <class TYPE> class Matrix;

  class PiezoMaterial : public BaseMaterial {

  public:

    //! Default constructor
    PiezoMaterial(MathParser* mp, CoordSystem * defaultCoosy);

    //! Destructor
    ~PiezoMaterial();

  };

} // end of namespace

#endif
