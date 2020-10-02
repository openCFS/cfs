// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     TestMaterial.hh
 *       \brief    Defines the material parameters for PDE TEST
 *
 *       \date     July 2, 2013
 *       \author   Manfred Kaltenbacher, TU Wien
 */
//================================================================================================

#ifndef TESTMATERIAL_DATA
#define TESTMATERIAL_DATA

#include "General/defs.hh"
#include "General/Environment.hh"
#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling material data for TEST PDE
  */

  class TestMaterial : public BaseMaterial {

  public:

    //! Default constructor
    TestMaterial(MathParser* mp, CoordSystem * defaultCoosy);

    //! Destructor
    ~TestMaterial();

  };

} // end of namespace

#endif
