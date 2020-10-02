// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ELECTRICCONDUCTIONMATERIAL_DATA
#define ELECTRICCONDUCTIONMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling electric current conduction material data
  */

  class ElectricConductionMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectricConductionMaterial(MathParser* mp, CoordSystem * defaultCoosy);

    //! Destructor
    ~ElectricConductionMaterial();

    //! Trigger finalization of material
    void Finalize();

    virtual PtrCoefFct GetScalCoefFncMultivariateNonLin(
        MaterialType matType,
        NonLinType nlType,
        Global::ComplexPart matDataType,
        StdVector<PtrCoefFct> dependencies,
        StdVector<RegionIdType> & regs);
  };

} // end of namespace

#endif
