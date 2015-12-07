// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ACOUSTICMATERIAL_DATA
#define ACOUSTICMATERIAL_DATA

#include "General/defs.hh"
#include "General/Environment.hh"
#include "BaseMaterial.hh"

namespace CoupledField {

  //! Class for Material Data
  /*! 
     Class for handling acoustic material data
  */

  class AcousticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    AcousticMaterial(MathParser* mp,
                     CoordSystem * defaultCoosy);

    //! Destructor
    ~AcousticMaterial();

    //! set a scalar real material parameter
    void SetScalar( Double param, MaterialType matType, 
		    Global::ComplexPart dataType );
    //! set a string material parameter
    void SetScalar( const std::string& param, MaterialType matType,
        Global::ComplexPart dataType );


    //! get a scalar real material parameter
    void GetScalar( Double& param, MaterialType matType, 
		    Global::ComplexPart dataType ) const;
    //! get a scalar real material parameter
    void GetScalar( std::string& param, MaterialType matType,
        Global::ComplexPart dataType ) const;


  private:


  };

} // end of namespace

#endif
