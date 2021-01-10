#ifndef ELECQUASISTATICMATERIAL_DATA
#define ELECQUASISTATICMATERIAL_DATA

#include "BaseMaterial.hh"

namespace CoupledField {


  //! Class for Material Data
  /*! 
     Class for handling electro-quasistatic material data
  */

    class ElecQuasistaticMaterial : public BaseMaterial {

  public:

    //! Default constructor
    ElecQuasistaticMaterial(MathParser* mp,
                          CoordSystem * defaultCoosy);

    //! Destructor
    ~ElecQuasistaticMaterial();

    //! Finalize material
    void Finalize();

  };

} // end of namespace

#endif
