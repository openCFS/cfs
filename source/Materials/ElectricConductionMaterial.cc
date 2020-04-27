// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "ElectricConductionMaterial.hh"
#include <Utils/LinInterpolate.hh>
#include "Domain/CoefFunction/CoefFunctionApprox.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"
#include "Utils/BiLinInterpolate.hh"
#include "Utils/TriLinInterpolate.hh"


namespace CoupledField
{

// ***********************
//   Default Constructor
// ***********************
ElectricConductionMaterial::ElectricConductionMaterial(MathParser* mp,
                                                       CoordSystem * defaultCoosy)
: BaseMaterial(ELECTRICCONDUCTION, mp, defaultCoosy)
{
  //set the allowed material parameters
  isAllowed_.insert( ELEC_CONDUCTIVITY_SCALAR );
  isAllowed_.insert( ELEC_CONDUCTIVITY_TENSOR );
  isAllowed_.insert( ELEC_CONDUCTIVITY_1 );
  isAllowed_.insert( ELEC_CONDUCTIVITY_2 );
  isAllowed_.insert( ELEC_CONDUCTIVITY_3 );
  isAllowed_.insert( NONLIN_DEPENDENCY );
}

ElectricConductionMaterial::~ElectricConductionMaterial() {
}

void ElectricConductionMaterial::Finalize() {
  MaterialType orthoProps[3] = {
      ELEC_CONDUCTIVITY_1,
      ELEC_CONDUCTIVITY_2,
      ELEC_CONDUCTIVITY_3
  };
  CalcFull3x3Tensor(ELEC_CONDUCTIVITY_SCALAR, orthoProps, ELEC_CONDUCTIVITY_TENSOR);
}

PtrCoefFct ElectricConductionMaterial::
GetScalCoefFncMultivariateNonLin(MaterialType matType,
                                 NonLinType nlType,
                                 Global::ComplexPart matDataType,
                                 StdVector<PtrCoefFct> dependencies,
                                 StdVector<RegionIdType> & regs)
{
  UInt ndeps = dependencies.GetSize(); // 1 no temp dep, 2: temp dep
  // for a tripole, we've got a dependency more because two differences can be calculated
  if (nlType == NLELEC_TRIPOLE || nlType == NLELEC_TRIPOLE_TEMP_DEP) {
    ndeps += 1;
  }
  UInt nregs = regs.GetSize(); // 2 regs dipole, 3 regs tripole

  // Ensure that only real-valued parameters are used
  if( matDataType != Global::REAL ) {
    EXCEPTION( "Only real-valued nonlinear parameters are supported");
  }
  
  // check if isotropic material type is defined
  if( nonlinIsoParams_.find(matType) == nonlinIsoParams_.end() ) {
    EXCEPTION( "No nonlinear definition found for material type '"
        << MaterialTypeEnum.ToString(matType) << "'");
  }
  
  // check, if nonlinear curve was already calculated ??
  // using matType...
  MatDescriptorNl & matNl = nonlinIsoParams_[matType];

  // Check, if smooth spline approximation was already created 
  // and initialized
  if( !matNl.approxData ) {
    if ( matNl.approxType == SMOOTH_SPLINES ) {
      SmoothSpline * sp = new SmoothSpline( matNl.fileName, matType );
      sp->SetAccuracy( matNl.measAccuracy );
      sp->SetMaxY( matNl.maxVal );
      sp->CalcBestParameter();
      sp->CalcApproximation();
      sp->Print();
      matNl.approxData = sp;
    }
    else if ( matNl.approxType == LIN_INTERPOLATE ) {
      LinInterpolate * sp = new LinInterpolate( matNl.fileName, matType );
      //sp->SetAccuracy( matNl.measAccuracy );
      //sp->SetMaxY( matNl.maxVal );
      sp->SetFactor(matNl.factor); 
      sp->Print();
      matNl.approxData = sp;
     }
     else if ( matNl.approxType == BILIN_INTERPOLATE ) {

      if ( ndeps != 2 )
        EXCEPTION("bilinear interpoltation requires two ptrCoefFcts");
      std::string val = stringParams_[NONLIN_DEPENDENCY];
      if  (  val == "voltage-temperature" || val == "voltage-voltage" ) {
        BiLinInterpolate * sp = new BiLinInterpolate( matNl.fileName, matType );
	sp->SetFactor(matNl.factor); 
        sp->Print();
        matNl.approxData = sp;
      } else {
        EXCEPTION("'dependency' string in material xml must be specified with value 'voltage-temperature' (or 'voltage-voltage')" << 
	" to force user check column ordering: | voltage | temperature | sigma |");
      }
    }
     else if ( matNl.approxType == TRILIN_INTERPOLATE ) {

      if ( ndeps != 3 )
        EXCEPTION("trilinear interpoltation requires two ptrCoefFcts");
      std::string val = stringParams_[NONLIN_DEPENDENCY];
      if  (  val == "voltage-voltage-temperature" || 
             val == "voltage-voltage-temperature-Vds" ||
             val == "voltage-voltage-temperature-Vds-ElemArea" ||
             val == "voltage-voltage-temperature-Vds-nearestVs" ||
             val == "voltage-voltage-temperature-Vds-nearestVd"
          ) {
        TriLinInterpolate * sp = new TriLinInterpolate( matNl.fileName, matType );
	sp->SetFactor(matNl.factor); 
        sp->Print();
        matNl.approxData = sp;
      } else {
        EXCEPTION("'dependency' string in material xml must be specified with value 'voltage-voltage-temperature' " << 
	" to force user check column ordering: | voltage | voltage | temperature | sigma |");
      }
    }
    else {
      EXCEPTION( "No nonlinear approx/interpolate type not available '"
          << ApproxCurveTypeEnum.ToString(matNl.approxType) << "'");
    }
  }

  ApproxData * sp = matNl.approxData;
  // get linear starting value
  Double startVal = 0.0;
  this->GetScalar( startVal, matType, Global::REAL );

  shared_ptr<CoefFunctionComposite> coef;
  if (nlType == NLELEC_BIPOLE) {
    if (ndeps != 1 || nregs != 2)
      EXCEPTION("For a nonlinear electric dipole, two regions and one dependent coef function must be provided");
    coef.reset(new CoefFunctionBipole());
    coef->SetRegion(ANODE,regs[0]);
    coef->SetRegion(CATHODE,regs[1]);
    coef->SetDependCoef(NLELEC_BIPOLE, dependencies[0] );
    coef->Init( startVal, sp);

  }
  else if (nlType == NLELEC_BIPOLE_TEMP_DEP) {
    if (ndeps != 2 || nregs != 2)
      EXCEPTION("For a nonlinear temperature dependent electric dipole, two regions and two dependent coef functions must be provided");
    coef.reset(new CoefFunctionHeatBipole());
    coef->SetRegion(ANODE,regs[0]);
    coef->SetRegion(CATHODE,regs[1]);
    coef->SetDependCoef(NLELEC_BIPOLE, dependencies[0] );
    coef->SetDependCoef(NLELEC_CONDUCTIVITY, dependencies[1] ); // this means temperature dep
    coef->Init( startVal, sp);

  }
  else if (nlType == NLELEC_TRIPOLE) {
    if (ndeps - 1 != 1 || nregs != 3) // - 1 because of extra dependency of multiple terminal
      EXCEPTION("For a nonlinear electric tripole, three regions and one dependent coef function must be provided");
    coef.reset(new CoefFunctionTripole());
    coef->SetRegion(GATE,regs[0]);
    coef->SetRegion(DRAIN,regs[1]);
    coef->SetRegion(SOURCE,regs[2]);
    coef->SetDependCoef(NLELEC_TRIPOLE, dependencies[0] );
    coef->Init( startVal, sp);
  }
  else if (nlType == NLELEC_TRIPOLE_TEMP_DEP) {
    if (ndeps -1 != 2 || nregs != 3)
      EXCEPTION("For a nonlinear temperature dependent electric tripole, three regions and two dependent coef functions (elec and heat) must be provided");
    coef.reset(new CoefFunctionHeatTripole());
    coef->SetRegion(GATE,regs[0]);
    coef->SetRegion(DRAIN,regs[1]);
    coef->SetRegion(SOURCE,regs[2]);
    coef->SetDependCoef(NLELEC_TRIPOLE, dependencies[0] );
    coef->SetDependCoef(NLELEC_CONDUCTIVITY, dependencies[1] ); // this means temperature dep
    coef->Init( startVal, sp);

  }
  else {
    EXCEPTION("non linear type not recognized");
  }
  std::string val = stringParams_[NONLIN_DEPENDENCY];
  if (val.find("Vds") != std::string::npos) {
    coef->SetDivideByVds(true);
  }
  if (val.find("nearestVd") != std::string::npos) {
    coef->SetLocValue(DRAIN);
  }
  if (val.find("nearestVs") != std::string::npos) {
    coef->SetLocValue(SOURCE);
  }
  if (val.find("ElemArea") != std::string::npos) {
    coef->SetElemAreaMult(DRAIN); // refer to drain area
  }

    
  return coef;
}

}
