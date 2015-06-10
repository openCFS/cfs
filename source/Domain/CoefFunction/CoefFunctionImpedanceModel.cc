#include "CoefFunctionImpedanceModel.hh"
#include "muParser.h"

#include <boost/math/special_functions/hankel.hpp>


namespace CoupledField{

  // =========================================================================
  // CoefFunctionBCross
  // =========================================================================

  CoefFunctionImpedanceModel<Complex>::CoefFunctionImpedanceModel(MathParser* mp, Double c0, Double density,
                                                                  Double innerR, Double outerR)
  : CoefFunctionTimeFreq<Complex>(mp) {
    // this type of coefficient is always constant
    dependType_ = GENERAL;
    isAnalytic_ = false;
    isComplex_ = true;
    dimType_ = SCALAR;
    c0_ = c0;
    density_ = density;
    innerR_ = innerR;
    outerR_ = outerR;
  }

  CoefFunctionImpedanceModel<Complex>::~CoefFunctionImpedanceModel() {

  }

  void CoefFunctionImpedanceModel<Complex>::Recalculate() {
//std::cerr << "Recalculate Impedance!!!!!!!!!" << std::endl;
    Complex Z_mpp, Z_cav, Z_wall;
    Complex tmp1, tmp2, tmp3;
    Double waveNum, shearWaveNum;
    const Double nu = 0.00000935;    // TODO: hard coded kinemtic viscosity (air)
    const Double slitWidth = 0.000095; // TODO: hard coded slit width
    const Complex i(0,1);
    const Double thick = 0.001;
    const Double sigma = 0.043;
    const Double rho = 1.205;

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
        constCoefScalar_ = Complex(1.0,0.0);
        std::cerr << "frequency zero!!!!!!!!!" << std::endl;
        return;
      }
      waveNum = 2*pi*f/c0_;
      shearWaveNum = slitWidth*sqrt(2*pi*f/(4.0*nu));

      tmp1 = i;
      tmp1 *= waveNum;
      tmp1 *= thick/sigma; // TODO: sigma is missing, the porosity of the perforated surface. And thickness t is missing
      tmp2 = shearWaveNum*sqrt(i);
      tmp3 = -tanh(tmp2);
      tmp3 /= tmp2;
      tmp3 += 1;
      Z_mpp = tmp1;
      Z_mpp /= tmp3;

      Z_mpp += 4*outerR_/(sigma*rho*c0_);

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
//std::cerr << __LINE__ << " Z_mpp " <<Z_mpp <<std::endl;
//std::cerr << __LINE__ << " Z_cav " <<Z_cav <<std::endl;
//std::cerr << __LINE__ << " constCoefScalar " <<constCoefScalar_ <<std::endl;
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
