// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef MAGNETOSTRIMATERIAL_DATA
#define MAGNETOSTRIMATERIAL_DATA

#include "General/defs.hh"
#include "General/Environment.hh"
#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Magnetostrictive Materials
  /*! 
     Class for handling magnetostrictive material data
  */

template <class TYPE> class Matrix;

  class MagStrictMaterial : public BaseMaterial {

  public:

    //! Default constructor
    MagStrictMaterial(MathParser* mp, CoordSystem * defaultCoosy);

    //! Destructor
    ~MagStrictMaterial();

  };

} // end of namespace

#endif





