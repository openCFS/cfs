#include "CoefFunctionImpedanceModel.hh"
#include "muParser.h"
#include "Materials/AcousticMaterial.hh"
#include "Materials/BaseMaterial.hh"

#include <boost/math/special_functions/hankel.hpp>


namespace CoupledField{

  // =========================================================================
  // CoefFunctionBCross
  // =========================================================================

  CoefFunctionImpedanceModel<Complex>::CoefFunctionImpedanceModel(MathParser* mp)
  : CoefFunctionTimeFreq<Complex>(mp) {
    // this type of coefficient is always constant
    dependType_ = GENERAL;
    isAnalytic_ = false;
    isComplex_ = true;
    dimType_ = SCALAR;
    currFrequ_ = 0;
    c0_ = 0;
  }

  CoefFunctionImpedanceModel<Complex>::~CoefFunctionImpedanceModel() {

  }
  void CoefFunctionImpedanceModel<Complex>::GenerateZylindricMpp(BaseMaterial* const material,
                                                                 Double innerR, Double outerR) {
    LocPointMapped lp_dummy;
    Double tmpBlkMod;
    // density_, DENSITY
    PtrCoefFct tmpPtFc = material->GetScalCoefFnc( DENSITY, Global::REAL );
    tmpPtFc->GetScalar(density_, lp_dummy);
    // nu_, KINEMATIC_VISCOSITY
    tmpPtFc = material->GetScalCoefFnc( KINEMATIC_VISCOSITY, Global::REAL );
    tmpPtFc->GetScalar(nu_, lp_dummy);
    // slitWidth_, SLITWIDTH
    tmpPtFc = material->GetScalCoefFnc( SLITWIDTH, Global::REAL );
    tmpPtFc->GetScalar(slitWidth_, lp_dummy);
    // mppThick_, MPP_THICKNESS
    tmpPtFc = material->GetScalCoefFnc( MPP_THICKNESS, Global::REAL );
    tmpPtFc->GetScalar(mppThick_, lp_dummy);
    // sigma_, POROSITY
    tmpPtFc = material->GetScalCoefFnc( POROSITY, Global::REAL );
    tmpPtFc->GetScalar(sigma_, lp_dummy);
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
    innerR_ = innerR;
    outerR_ = outerR;
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate() {
    Complex Z_mpp, Z_cav, Z_wall;
    Complex tmp1, tmp2, tmp3;
    Double waveNum, shearWaveNum;
    const Complex i(0,1);
    const Double mu = nu_*density_;
    Double R_s; // surface resistance
    Double maxSchnelle; // maximal particle velocity

    if (c0_ == 0) {
      EXCEPTION("GenerateZylindricMpp(...) has to be called for");
    }

    if (dimType_ == SCALAR) {
      //calc Z_mpp
      //    i*2*pi*f*1/(1-tanh(2*pi*f/c*sqrt(i))/(2*pi*f/c*sqrt(i)))
      MathParser::HandleType h = mp_->GetNewHandle();
      const Double pi = mp_->GetExprVars(h, "_pi");
      Double f = mp_->GetExprVars(h, "f");
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
      waveNum = 2*pi*f/c0_;
      shearWaveNum = slitWidth_*sqrt(2*pi*f/(4.0*nu_));
      R_s = 0.5*sqrt(2*density_*2*pi*f*mu);
      maxSchnelle = 1.0 / ( c0_ * density_ ); // p / Z = p / (c * \rho)  // p = 1 <= Dirichlet BC

      tmp1 = i;
      tmp1 *= waveNum;
      tmp1 *= mppThick_/sigma_;
      tmp2 = shearWaveNum*sqrt(i);
      tmp3 = -tanh(tmp2);
      tmp3 /= tmp2;
      tmp3 += 1;
      Z_mpp = tmp1;
      Z_mpp /= tmp3;

      Z_mpp += 4*R_s/(sigma_*density_*c0_);
      Z_mpp += maxSchnelle/(sigma_ * c0_);
      Z_mpp += beta_*flowMachNr_/sigma_;

      // calc Z_cav
      tmp1 = boost::math::cyl_hankel_1(1,waveNum*outerR_);
      tmp1 /= boost::math::cyl_hankel_2(1,waveNum*outerR_);
      Z_cav = boost::math::cyl_hankel_2(0,waveNum*innerR_);
      Z_cav *= -tmp1;
      Z_cav += boost::math::cyl_hankel_1(0,waveNum*innerR_);

      tmp3 = boost::math::cyl_hankel_2(1,waveNum*innerR_);
      tmp3 *= -tmp1;
      tmp3 += boost::math::cyl_hankel_1(1,waveNum*innerR_);

      Z_cav *= i;
      Z_cav /= tmp3;

      Z_wall = Z_mpp;
      Z_wall += Z_cav;
      constCoefScalar_ = 1;
      constCoefScalar_ *= -i;
      constCoefScalar_ /= (Z_wall*density_*c0_);
    } else {
      EXCEPTION("CoefFunctionImpedanceMode only implemented for scalar")
    }

  }

  void CoefFunctionImpedanceModel<Complex>::
  GetStrScalar( std::string& real, std::string& imag ) {
    real = "";
    imag = "i*2*pi*f*1/(1-tanh(2*pi*f/c0*sqrt(i))/(2*pi*f/c0*sqrt(i))) + ";
  }
}
