#ifndef ELECTROSTATICMATERIAL_DATA
#define ELECTROSTATICMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {


  //! Class for Material Data
  /*! 
     Class for handling electrostatic material data
  */

  class Hysteresis;

  class ElectroStaticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElectroStaticMaterial(MathParser* mp,
                          CoordSystem * defaultCoosy);

    //! Destructor
    ~ElectroStaticMaterial();

    //! Finalize material
    void Finalize();

//
//    //======================== for hysteresis =======================================
//    //! compute scalar differential parameter
//    Double ComputeScalarDiffVal( UInt nrElem, Double Xval );
//
//
//    //set values for differential material approach
//    void SetPreviousHystVal( UInt nrElem, Double Xval );
//
//    Double ComputeScalarHystVal( UInt nrElem, Double Xval );

  };

} // end of namespace

#endif
