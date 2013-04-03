// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FLOWMATERIAL_DATA
#define FLOWMATERIAL_DATA

#include "General/defs.hh"
#include "General/Environment.hh"
#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling flow material data
  */

  class FlowMaterial : public BaseMaterial {

  public:

    //! Default constructor
    FlowMaterial(MathParser* mp,
                 CoordSystem * defaultCoosy);

    //! Destructor
    ~FlowMaterial();

    //! Trigger finalization of flow material (calculation of kinematic viscosity)
    void Finalize();

    //! Expand scalar variable to tensor
    void GenerateTensor( MaterialType matType, PtrCoefFct scalar );
    
  private:
    void ComputeAllViscosities();
  };

} // end of namespace

#endif
