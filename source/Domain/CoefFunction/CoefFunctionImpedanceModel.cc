#include "CoefFunctionImpedanceModel.hh"
#include "muParser.h"
#include "Materials/AcousticMaterial.hh"
#include "Materials/BaseMaterial.hh"

#include <boost/math/special_functions/hankel.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace CoupledField{

  // =========================================================================
  // CoefFunctionBCross
  // =========================================================================

  CoefFunctionImpedanceModel<Complex>::CoefFunctionImpedanceModel(MathParser* mp, \
      bool isNormalised)
  : CoefFunctionTimeFreq<Complex>(mp),
    isNormalised_(isNormalised),
    normalisedFactor_(1.0)
  {
    // this type of coefficient is always constant
    dependType_ = GENERAL;
    isAnalytic_ = false;
    isComplex_ = true;
    dimType_ = SCALAR;
    currFrequ_ = 0;
    c0_ = 0;

    impedanceType_ = IMP_NONE;
  }

  CoefFunctionImpedanceModel<Complex>::~CoefFunctionImpedanceModel() {

  }
  void CoefFunctionImpedanceModel<Complex>::GenerateSlitMpp(BaseMaterial* const material) {
    impedanceType_ = IMP_SLIT_MPP;
    Init(material);
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate_slitMpp() {
    Complex Z_mpp, Z_cav, Z_all;
    Complex tmp1, tmp2, tmp3;
    Double xi, omega;
    const Complex one_i(0,1);
    const Double eta = nu_ * density_;
    Double R_s; // surface resistance
    //Double maxSchnelle; // maximal particle velocity

    if (c0_ == 0) {
      EXCEPTION("GenerateSlitMpp(...) has to be called first");
    }

    if (dimType_ == SCALAR) {
      //calc Z_mpp
      //    i*2*pi*f*1/(1-tanh(2*pi*f/c*sqrt(i))/(2*pi*f/c*sqrt(i)))
      MathParser::HandleType h = mp_->GetNewHandle();
      const Double pi = mp_->GetExprVars(h, "_pi");
      const Double f = mp_->GetExprVars(h, "f");
      mp_->ReleaseHandle(h);

      if (f == 0)
      {
        std::cerr << "frequency zero!!!!!!!!!" << std::endl;
        return;
      }
      if (f == currFrequ_) // Already calculated value for this frequency
      {
        return;
      } else {
        currFrequ_ = f;
      }
      omega = 2 * pi * f;
      xi = holeDiam_*sqrt(omega / (4.0*nu_));
      R_s = 0.5*sqrt(2*density_*omega*eta);
      //maxSchnelle = 1.0 / ( c0_ * density_ ); // p / Z = p / (c * \rho)  // p = 1 <= Dirichlet BC

      tmp1 = one_i * omega * density_;
      tmp1 *= plateThick_/sigma_;
      tmp2 = xi * sqrt(one_i);
      tmp3 = -tanh(tmp2);
      tmp3 /= tmp2;
      tmp3 += 1;
      Z_mpp = tmp1/tmp3;

      Z_mpp += 4*R_s/sigma_;
      //Z_mpp += maxSchnelle * density_;
      //Z_mpp += beta_*flowMachNr_ * density_ * c0_;

#if 0
      // calc Z_cav
      tmp1 = boost::math::cyl_hankel_1(1,waveNum*outerR_);
      tmp1 /= boost::math::cyl_hankel_2(1,waveNum*outerR_);
      Z_cav = boost::math::cyl_hankel_2(0,waveNum*innerR_);
      Z_cav *= -tmp1;
      Z_cav += boost::math::cyl_hankel_1(0,waveNum*innerR_);

      tmp3 = boost::math::cyl_hankel_2(1,waveNum*innerR_);
      tmp3 *= -tmp1;
      tmp3 += boost::math::cyl_hankel_1(1,waveNum*innerR_);

      Z_cav *= one_i;
      Z_cav /= tmp3;
      Z_cav *= density_ * c0_;
#endif
      Calculate_cavityImpedance(Z_cav, omega);
      Z_all = Z_mpp + Z_cav;
      //constCoefScalar_ = - 1.0;
      //constCoefScalar_ *= one_i *  waveNum;
      constCoefScalar_ = normalisedFactor_ * density_ / Z_all;
    } else {
      EXCEPTION("CoefFunctionImpedanceMode only implemented for scalar")
    }

  }

  void CoefFunctionImpedanceModel<Complex>::GenerateCircMpp(BaseMaterial* const material) {
    impedanceType_ = IMP_CIRC_MPP;
    Init(material);
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate_circMpp() {
    EXCEPTION("bessel does not work!");
#if 0
    Complex Z_mpp, Z_all;
    Complex tmp1, tmp2, tmp3;
    Double omega, xi;
    const Complex one_i(0,1);
    const Double eta = nu_ * density_;
    Double R_s; // surface resistance
    //Double maxSchnelle; // maximal particle velocity

    if (c0_ == 0) {
      EXCEPTION("GenerateCircularMpp(...) has to be called first");
    }

    if (dimType_ == SCALAR) {
      //calc Z_mpp
      //    i*2*pi*f*1/(1-tanh(2*pi*f/c*sqrt(i))/(2*pi*f/c*sqrt(i)))
      MathParser::HandleType h = mp_->GetNewHandle();
      const Double pi = mp_->GetExprVars(h, "_pi");
      const Double f = mp_->GetExprVars(h, "f");
      mp_->ReleaseHandle(h);

      if (f == 0)
      {
        std::cerr << "frequency zero!!!!!!!!!" << std::endl;
        return;
      }
      if (f == currFrequ_) // Already calculated value for this frequency
      {
        return;
      } else {
        currFrequ_ = f;
      }
      omega = 2*pi*f;
      xi = holeDiam_*sqrt(omega/(4.0*nu_));
      R_s = 0.5*sqrt(2*density_*omega*eta);

      tmp1 = one_i * omega * plateThick_ * density_;
      tmp1 /= sigma_;
      tmp2 = (xi * sqrt(-one_i));
      tmp3 = boost::math::cyl_bessel_j(1, tmp2);
      tmp3 /= boost::math::cyl_bessel_j(0, tmp2);
      tmp3 /= tmp2;
      tmp3 *= 2.0;
      tmp3 += 1;
      Z_mpp = tmp1/tmp3;

      Z_mpp += 4*R_s/(sigma_ );

      //constCoefScalar_ = - 1.0;
      //constCoefScalar_ *= one_i *  waveNum;
      constCoefScalar_ = density_ / Z_mpp;
    } else {
      EXCEPTION("CoefFunctionImpedanceMode only implemented for scalar")
    }
#endif
  }


  void CoefFunctionImpedanceModel<Complex>::GenerateInterpolImpedance(BaseMaterial* const material) {
    impedanceType_ = IMP_INTERPOL;
    Double tmpBlkMod;

    LocPointMapped lp_dummy;
    // density_, DENSITY
    PtrCoefFct tmpPtFc = material->GetScalCoefFnc( DENSITY, Global::REAL );
    tmpPtFc->GetScalar(density_, lp_dummy);
    // tmpBlkMod, ACOU_BULK_MODULUS
    tmpPtFc = material->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    tmpPtFc->GetScalar(tmpBlkMod, lp_dummy);

    PtrCoefFct frequCoef = CoefFunction::Generate( mp_, Global::REAL, "f");
    impedanceCoef_real_ =
        material->GetScalCoefFncNonLin( ACOU_IMPEDANCE_REAL_VAL, Global::REAL, frequCoef);
    impedanceCoef_imag_ =
        material->GetScalCoefFncNonLin( ACOU_IMPEDANCE_IMAG_VAL, Global::REAL, frequCoef);

    // calc speed of sound
    c0_ = sqrt(tmpBlkMod/density_);
    if (isNormalised_) {
      normalisedFactor_ = 1.0 / (density_ * c0_);
    }
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate_interpol() {
    MathParser::HandleType h = mp_->GetNewHandle();
    const Double f = mp_->GetExprVars(h, "f");
    mp_->ReleaseHandle(h);
    if (f == currFrequ_) // Already calculated value for this frequency
    {
      return;
    } else {
      currFrequ_ = f;
    }

    Double Z_real, Z_imag;
    LocPointMapped lp_dummy;
    impedanceCoef_real_->GetScalar(Z_real, lp_dummy);
    impedanceCoef_imag_->GetScalar(Z_imag, lp_dummy);
    const Complex Z(Z_real, Z_imag);
    constCoefScalar_ = normalisedFactor_ * density_ / Z;
  }

  void CoefFunctionImpedanceModel<Complex>::GenerateImpedanceFct(BaseMaterial* const material) {
    impedanceType_ = IMP_FCT;
    LocPointMapped lp_dummy;
    Double tmpBlkMod;

    //Impedance function for real and imag part
    impedanceCoef_real_ = material->GetScalCoefFnc( ACOU_IMPEDANCE_REAL_VAL, Global::REAL );
    impedanceCoef_imag_ = material->GetScalCoefFnc( ACOU_IMPEDANCE_IMAG_VAL, Global::REAL );

    // density_, DENSITY
    PtrCoefFct tmpPtFc = material->GetScalCoefFnc( DENSITY, Global::REAL );
    tmpPtFc->GetScalar(density_, lp_dummy);
    // tmpBlkMod, ACOU_BULK_MODULUS
    tmpPtFc = material->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    tmpPtFc->GetScalar(tmpBlkMod, lp_dummy);

    // calc speed of sound
    c0_ = sqrt(tmpBlkMod/density_);
    if (isNormalised_) {
      normalisedFactor_ = 1.0 / (density_ * c0_);
    }
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate_impFct() {
    MathParser::HandleType h = mp_->GetNewHandle();
    const Double f = mp_->GetExprVars(h, "f");
    mp_->ReleaseHandle(h);
    if (f == currFrequ_) // Already calculated value for this frequency
    {
      return;
    } else {
      currFrequ_ = f;
    }

    Double Z_real, Z_imag;
    LocPointMapped lp_dummy;
    impedanceCoef_real_->GetScalar(Z_real, lp_dummy);
    impedanceCoef_imag_->GetScalar(Z_imag, lp_dummy);

    const Complex Z(Z_real, Z_imag);
    constCoefScalar_ = normalisedFactor_ * density_ / Z;
  }

  void CoefFunctionImpedanceModel<Complex>::Init(BaseMaterial* const material)
  {
    LocPointMapped lp_dummy;
    Double tmpBlkMod;
    // density_, DENSITY
    PtrCoefFct tmpPtFc = material->GetScalCoefFnc( DENSITY, Global::REAL );
    tmpPtFc->GetScalar(density_, lp_dummy);
    // nu_, KINEMATIC_VISCOSITY
    tmpPtFc = material->GetScalCoefFnc( KINEMATIC_VISCOSITY, Global::REAL );
    tmpPtFc->GetScalar(nu_, lp_dummy);
    // holeDiam_, HOLEDIAM
    tmpPtFc = material->GetScalCoefFnc( IMP_HOLE_DIAM, Global::REAL );
    tmpPtFc->GetScalar(holeDiam_, lp_dummy);
    // plateThick_, IMP_PLATE_THICKNESS
    tmpPtFc = material->GetScalCoefFnc( IMP_PLATE_THICKNESS, Global::REAL );
    tmpPtFc->GetScalar(plateThick_, lp_dummy);
    // sigma_, POROSITY
    tmpPtFc = material->GetScalCoefFnc( POROSITY, Global::REAL );
    tmpPtFc->GetScalar(sigma_, lp_dummy);
    // mppVolDepth_, MPP_VOL_DEPTH
    tmpPtFc = material->GetScalCoefFnc( MPP_VOLUME_DEPTH, Global::REAL );
    tmpPtFc->GetScalar(mppVolDepth_, lp_dummy);
    // flowMachNr_, FLOW_MACH_NUMBER
    tmpPtFc = material->GetScalCoefFnc( FLOW_MACH_NUMBER, Global::REAL );
    tmpPtFc->GetScalar(flowMachNr_, lp_dummy);
    // beta_, BETA
    tmpPtFc = material->GetScalCoefFnc( BETA, Global::REAL );
    tmpPtFc->GetScalar(beta_, lp_dummy);
    // tmpBlkMod, ACOU_BULK_MODULUS
    tmpPtFc = material->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    tmpPtFc->GetScalar(tmpBlkMod, lp_dummy);
    // calc speed of sound
    c0_ = sqrt(tmpBlkMod/density_);
    if (isNormalised_) {
      normalisedFactor_ = 1.0 / (density_ * c0_);
    }
  }

  inline void CoefFunctionImpedanceModel<Complex>::Calculate_cavityImpedance(Complex& Z_cav, const Double omega) {
    const Complex one_i(0,1);
    Z_cav = -one_i/tan(omega*mppVolDepth_/c0_);
  }

  void CoefFunctionImpedanceModel<Complex>::
  GetStrScalar( std::string& real, std::string& imag ) {
    real = "";
    imag = "i*2*pi*f*1/(1-tanh(2*pi*f/c0*sqrt(i))/(2*pi*f/c0*sqrt(i))) + ";
  }
}
