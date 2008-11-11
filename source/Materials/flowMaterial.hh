// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FLOWMATERIAL_DATA
#define FLOWMATERIAL_DATA

#include "baseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling flow material data
  */

  class FlowMaterial : public BaseMaterial {

  public:

    //! Default constructor
    FlowMaterial();

    //! Destructor
    ~FlowMaterial();

    //! Trigger finalization of flow material (calculation of kinematic viscosity)
    void Finalize();

    
    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    DataType dataType );


    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    DataType dataType ) const;

  private:
    void ComputeAllViscosities();
  };

} // end of namespace

#endif
