// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "AcousticMaterial.hh"

#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

namespace CoupledField {

  // ***********************
  //   Default Constructor
  // ***********************
  AcousticMaterial::AcousticMaterial(MathParser *mp, CoordSystem *defaultCoosy) : BaseMaterial(ACOUSTIC, mp, defaultCoosy)
  {
    //set the allowed material parameters
    isAllowed_.insert(DENSITY);
    isAllowed_.insert(ACOU_BULK_MODULUS);
    isAllowed_.insert(ACOU_SOUND_SPEED);
    isAllowed_.insert(FLUID_ADIABATIC_EXPONENT);
    isAllowed_.insert(FLUID_KINEMATIC_VISCOSITY);
    isAllowed_.insert(RAYLEIGH_ALPHA);
    isAllowed_.insert(RAYLEIGH_BETA);
    isAllowed_.insert(LOSS_TANGENS_DELTA);
    isAllowed_.insert(ACOU_BOVERA);
    isAllowed_.insert(ACOU_TDEF_INVDENS_CONST);
    isAllowed_.insert(ACOU_TDEF_INVDENS_A);
    isAllowed_.insert(ACOU_TDEF_INVDENS_ALPHA);
    isAllowed_.insert(ACOU_TDEF_INVDENS_B);
    isAllowed_.insert(ACOU_TDEF_INVDENS_BETA);
    isAllowed_.insert(ACOU_TDEF_INVDENS_C);
    isAllowed_.insert(ACOU_TDEF_INVDENS_GAMMA);
    isAllowed_.insert(ACOU_TDEF_INVBLK_CONST);
    isAllowed_.insert(ACOU_TDEF_INVBLK_A);
    isAllowed_.insert(ACOU_TDEF_INVBLK_ALPHA);
    isAllowed_.insert(ACOU_TDEF_INVBLK_B);
    isAllowed_.insert(ACOU_TDEF_INVBLK_BETA);
    isAllowed_.insert(ACOU_TDEF_INVBLK_C);
    isAllowed_.insert(ACOU_TDEF_INVBLK_GAMMA);
  }

  AcousticMaterial::~AcousticMaterial() {
  }

  void AcousticMaterial::Finalize() {
    if (scalarCoef_.find(ACOU_SOUND_SPEED) == scalarCoef_.end()) {
      CoefMap::const_iterator densIt = scalarCoef_.find(DENSITY),
                              bulkIt = scalarCoef_.find(ACOU_BULK_MODULUS);
      if (densIt != scalarCoef_.end() && bulkIt != scalarCoef_.end()) {
        SetCoefFct(ACOU_SOUND_SPEED, CoefFunction::Generate( mp_,  Global::REAL,
            CoefXprUnaryOp(mp_, CoefXprBinOp(mp_, bulkIt->second, densIt->second,
                                             CoefXpr::OP_DIV),
                           CoefXpr::OP_SQRT)));
      }
      else {
        EXCEPTION("Cannot not compute'" << MaterialTypeEnum.ToString(ACOU_SOUND_SPEED)
                  << "' for material '" << name_ << "' due to missing data");
      }
    }
  }
}
