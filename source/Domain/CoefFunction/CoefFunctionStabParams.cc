//================================================================================================
/*!
 *       \file     CoefFunctionStabParams.cc
 *       \brief    Coefficient function which computes the stabilization parameters for SUPG
 *
 *       \date     07/04/2013
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#include "CoefFunctionStabParams.hh"

namespace CoupledField{

CoefFunctionStabParams::CoefFunctionStabParams(PtrCoefFct density,
							PtrCoefFct viscosity) : CoefFunction(){

	dimType_ = SCALAR;
	isAnalytic_ = false;
	dependType_ = GENERAL;

	density_ = density;
	viscosity_ = viscosity;

}

void CoefFunctionStabParams::GetScalar(Double& scal, const LocPointMapped& lpm ) {

	Double dens, vis;
	density_->GetScalar(dens,lpm);
	viscosity_->GetScalar(vis,lpm);

	Double vol = lpm.shapeMap->CalcVolume();
	UInt dim = lpm.jac.GetNumRows();
	if ( dim ==2 )
		scal = ( dens / (12.0*vis) ) * vol;
	else
		scal = ( dens / (12.0*vis) ) * std::pow(vol, 2.0/3.0);

//	std::cout << "\n StabVal:" << scal << std::endl;
}

}// end of namespace
