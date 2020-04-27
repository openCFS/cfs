// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef HEATMATERIAL_DATA
#define HEATMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling heat material data
  */

  class HeatMaterial : public BaseMaterial {

  public:

    //! Default constructor
    HeatMaterial(MathParser* mp,
                 CoordSystem * defaultCoosy);

    //! Destructor
    ~HeatMaterial();

    //! Trigger finalization of material
    void Finalize();

  };

} // end of namespace

#endif
