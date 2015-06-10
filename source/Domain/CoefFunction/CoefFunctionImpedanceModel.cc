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
    PtrCoefFct tmpPtFc = material->GetScalCoefFnc( DENSITY, Global::REAL );
    tmpPtFc->GetScalar(density_, lp_dummy);
    //material->GetScalar(density_, DENSITY, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( KINEMATIC_VISCOSITY, Global::REAL );
    tmpPtFc->GetScalar(nu_, lp_dummy);
    //material->GetScalar(nu_, KINEMATIC_VISCOSITY, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( SLITWIDTH, Global::REAL );
    tmpPtFc->GetScalar(slitWidth_, lp_dummy);
    //material->GetScalar(slitWidth_, SLITWIDTH, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( MPP_THICKNESS, Global::REAL );
    tmpPtFc->GetScalar(mppThick_, lp_dummy);
    //material->GetScalar(mppThick_, MPP_THICKNESS, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( POROSITY, Global::REAL );
    tmpPtFc->GetScalar(sigma_, lp_dummy);
    //material->GetScalar(sigma_, POROSITY, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( FLOW_MACH_NUMBER, Global::REAL );
    tmpPtFc->GetScalar(flowMachNr_, lp_dummy);
    //material->GetScalar(flowMachNr_, FLOW_MACH_NUMBER, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( BETA, Global::REAL );
    tmpPtFc->GetScalar(beta_, lp_dummy);
    //material->GetScalar(beta_, BETA, Global::REAL);
    tmpPtFc = material->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    tmpPtFc->GetScalar(tmpBlkMod, lp_dummy);
    //material->GetScalar(tmpBlkMod, ACOU_BULK_MODULUS, Global::REAL);
    c0_ = sqrt(tmpBlkMod/density_);
    innerR_ = innerR;
    outerR_ = outerR;
  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate() {
//std::cerr << "Recalculate Impedance!!!!!!!!!" << std::endl;
    Complex Z_mpp, Z_cav, Z_wall;
    Complex tmp1, tmp2, tmp3;
    Double waveNum, shearWaveNum;
    const Complex i(0,1);
    /* const Double nu = 0.00001511;    // TODO: hard coded kinematic viscosity (air)
    const Double slitWidth = 0.000095; // TODO: hard coded slit width
    const Double thick = 0.001;
    const Double sigma = 0.043;
    const Double rho = 1.21;
    const Double machNr = 0.05;
    const Double beta = 0.1; */
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
      //StdVector<std::string> varNames;
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
      tmp1 *= mppThick_/sigma_; // TODO: sigma is missing, the porosity of the perforated surface. And thickness is missing
      tmp2 = shearWaveNum*sqrt(i);
      tmp3 = -tanh(tmp2);
      tmp3 /= tmp2;
      tmp3 += 1;
      Z_mpp = tmp1;
      Z_mpp /= tmp3;

      Z_mpp += 4*R_s/(sigma_*density_*c0_);
      Z_mpp += maxSchnelle/(sigma_ * c0_);
      Z_mpp += beta_*flowMachNr_/sigma_;

      // TODO: Z_mpp + velocity terms + Mach number term

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
std::cerr << std::endl;
std::cerr << __LINE__ << " Z_mpp " <<Z_mpp <<std::endl;
std::cerr << __LINE__ << " Z_cav " <<Z_cav <<std::endl;
std::cerr << __LINE__ << " constCoefScalar " <<constCoefScalar_ <<std::endl;
std::cerr << __LINE__ << " Resistance " <<Z_wall.real() <<std::endl;
std::cerr << __LINE__ << " Reactance " <<Z_wall.imag() <<std::endl;
std::cerr << std::endl;
std::cerr << __LINE__ << " sigma_*c0_ " <<sigma_*c0_ <<std::endl;
std::cerr << __LINE__ << " maxSchnelle " <<maxSchnelle <<std::endl;


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
