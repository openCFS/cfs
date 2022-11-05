#include "ElecQuasistaticMaterial.hh"

#include "Domain/CoefFunction/CoefFunction.hh"

namespace CoupledField
{

  // ***********************
  //   Default Constructor
  // ***********************
  ElecQuasistaticMaterial::ElecQuasistaticMaterial(MathParser* mp,
                                               CoordSystem * defaultCoosy) 
  : BaseMaterial(ELECQUASISTATIC, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert( ELEC_CONDUCTIVITY_SCALAR );
    isAllowed_.insert( ELEC_CONDUCTIVITY_TENSOR );
    isAllowed_.insert( ELEC_CONDUCTIVITY_1 );
    isAllowed_.insert( ELEC_CONDUCTIVITY_2 );
    isAllowed_.insert( ELEC_CONDUCTIVITY_3 );

    isAllowed_.insert( ELEC_PERMITTIVITY_SCALAR );
    isAllowed_.insert( ELEC_PERMITTIVITY_TENSOR );
    isAllowed_.insert( ELEC_PERMITTIVITY_1 );
    isAllowed_.insert( ELEC_PERMITTIVITY_2 );
    isAllowed_.insert( ELEC_PERMITTIVITY_3 );
  }

  ElecQuasistaticMaterial::~ElecQuasistaticMaterial() {
  }

  void ElecQuasistaticMaterial::Finalize() {
    MaterialType orthoCond[3] = {
       ELEC_CONDUCTIVITY_1,
       ELEC_CONDUCTIVITY_2,
       ELEC_CONDUCTIVITY_3
    };
    CalcFull3x3Tensor(ELEC_CONDUCTIVITY_SCALAR, orthoCond, ELEC_CONDUCTIVITY_TENSOR);

    MaterialType orthoPerm[3] = {
        ELEC_PERMITTIVITY_1,
        ELEC_PERMITTIVITY_2,
        ELEC_PERMITTIVITY_3
    };
    CalcFull3x3Tensor(ELEC_PERMITTIVITY_SCALAR, orthoPerm, ELEC_PERMITTIVITY_TENSOR);
  }

}
