#include "CoefFunctionApprox.hh"

// classes for function / spline approximation
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField {

  DECLARE_LOG(coeffctapprox)
  DEFINE_LOG(coeffctapprox, "coeffctapprox")
// ============================================================================
//  ISOTROPIC VERSIONS
// ============================================================================

CoefFunctionApprox::CoefFunctionApprox() : CoefFunction() {
  // this type of coefficient is nonlinear (i.e. solution dependend)
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_ = false;
}

CoefFunctionApprox::~CoefFunctionApprox(){
  ;
}

void CoefFunctionApprox::Init( Double coefScalar, ApproxData * nLinFnc,
                               PtrCoefFct dependCoef ) {

  // set type to scalar
  dimType_ = SCALAR;
  nLinFnc_ = nLinFnc;
  coefScalar_ = coefScalar;
  dependCoef_ = dependCoef;
}

//! \see CoefFunction::GetScalar
void CoefFunctionApprox::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {

  
  // evaluate vector of dependency
  Vector<Double> elemSol;
  dependCoef_->GetVector( elemSol, lpm);
  
  if ( nLinFnc_->GetMatType() == MAG_PERMEABILITY ) {
    Double fieldAbs = elemSol.NormL2();

    if( fieldAbs == 0 ) { 
      coefScalar = coefScalar_;
    } else {
      coefScalar = nLinFnc_->EvaluateFuncNu(fieldAbs);
    }
  }
  else {
    coefScalar = nLinFnc_->EvaluateFunc(elemSol[0]);
  }
  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal = '" << elemSol[0] << ". IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;
}

bool IsComplex(){
  return false;
}

std::string CoefFunctionApprox::ToString() const {
  return "";
}
// ============================================================================
//  coef function composite
// ============================================================================

CoefFunctionComposite::CoefFunctionComposite() : CoefFunction() {

  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_ = false;
}

void CoefFunctionComposite::SetRegion(TerminalConnector tc, RegionIdType reg){

  terminals_[tc] = reg;
  nRegions_ = terminals_.size();

}
void CoefFunctionComposite::SetDependCoef(NonLinType nl, PtrCoefFct dep ) {

  dependCoefs_[nl] = dep;
  nDepCoefs_ = dependCoefs_.size();

}

Double CoefFunctionComposite::GetAvgTerminalValue(TerminalConnector tc, NonLinType nl, const LocPointMapped & lpm) {

  Grid * ptGrid = lpm.GetShapeMap()->GetGrid();
  StdVector <Elem*> elems;
  ptGrid->GetElems(elems, terminals_[tc]);
  Vector<Double> ElemAvg;
  ElemAvg.Resize(elems.GetSize());
  Double avg = 0.;
  for (UInt i = 0; i < elems.GetSize(); i++) {
    dependCoefs_[nl]->GetAvgElemValue(ElemAvg[i], elems[i]);
    avg += ElemAvg[i];
  }
  avg /= elems.GetSize();
  return avg;
}

void CoefFunctionComposite::Init( Double coefScalar, ApproxData * nLinFnc) {

  dimType_ = SCALAR;
  nLinFnc_ = nLinFnc;
  coefScalar_ = coefScalar;

}

// ========================================================================
//
// ============================================================================
//  Coef Function Bipole
// ============================================================================
//
void CoefFunctionBipole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {

  Double aAvg = GetAvgTerminalValue(ANODE, NLELEC_BIPOLE, lpm);
  Double bAvg = GetAvgTerminalValue(CATHODE, NLELEC_BIPOLE, lpm);
  Double diff = bAvg - aAvg;
  
  coefScalar = nLinFnc_->EvaluateFunc(diff);
  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal (diff) = '" << diff << ". IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;

}
//
// ============================================================================
//  Coef Function Heat Bipole
// ============================================================================
//
void CoefFunctionHeatBipole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {
  Double aAvg = GetAvgTerminalValue(ANODE, NLELEC_BIPOLE, lpm);
  Double bAvg = GetAvgTerminalValue(CATHODE, NLELEC_BIPOLE, lpm);
  Double diff = bAvg - aAvg;

  // get temperature
  Vector<Double> elemSol;
  dependCoefs_[NLELEC_CONDUCTIVITY]->GetVector( elemSol, lpm);


  coefScalar = nLinFnc_->EvaluateFunc(diff, elemSol[0]);
  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal (diff) = '" << diff << " and temperature = '" << elemSol[0] <<"'. IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;

}

// ========================================================================
//
// ============================================================================
//  Coef Function Tripole
// ============================================================================
//
void CoefFunctionTripole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {

  Double Vdrain = GetAvgTerminalValue(DRAIN, NLELEC_TRIPOLE, lpm);
  Double Vsource = GetAvgTerminalValue(SOURCE, NLELEC_TRIPOLE, lpm);
  Double Vgate = GetAvgTerminalValue(GATE, NLELEC_TRIPOLE, lpm);
  Double Vds = Vdrain-Vsource;
  Double Vgs = Vgate-Vsource;
   
  coefScalar = nLinFnc_->EvaluateFunc(Vgs, Vds);
  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal vgs = '" << Vgs << " and vds = '" << Vds <<"'. IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;

}
//
// ============================================================================
//  Coef Function Heat Tripole
// ============================================================================
//
void CoefFunctionHeatTripole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {
  Double Vdrain = GetAvgTerminalValue(DRAIN, NLELEC_TRIPOLE, lpm);
  Double Vsource = GetAvgTerminalValue(SOURCE, NLELEC_TRIPOLE, lpm);
  Double Vgate = GetAvgTerminalValue(GATE, NLELEC_TRIPOLE, lpm);
  Double Vds = Vdrain-Vsource;
  Double Vgs = Vgate-Vsource;

  // get temperature
  Vector<Double> elemSol;
  dependCoefs_[NLELEC_CONDUCTIVITY]->GetVector( elemSol, lpm);

  // special case: my material actually returns current and thus has to be divided by Vds 
  if (fabs(Vds) < 1e-9) { // Volts
    coefScalar = 1e-9;
  } else {
    coefScalar = 1./Vds*nLinFnc_->EvaluateFunc(Vds, Vgs, elemSol[0]);
  }
  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal vgs = '" << Vgs << " and vds = '" << Vds <<"' and temperature = '" << elemSol[0] <<"'. IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;

}

// ========================================================================

CoefFunctionApproxDeriv::CoefFunctionApproxDeriv() : CoefFunction() {
  // this type of coefficient is nonlinear, i.e. spatial and time dependent
  dependType_ = GENERAL;
  isComplex_ = false;
}

CoefFunctionApproxDeriv::~CoefFunctionApproxDeriv(){
}

void CoefFunctionApproxDeriv::Init( ApproxData * nLinFnc,
                                    UInt dimDMat,
                                    PtrCoefFct dependCoef ) {

  // set type to TENSOR
  dimType_ = TENSOR;
  nLinFnc_ = nLinFnc;
  dimDMat_ = dimDMat;
  this->dependCoef_ = dependCoef;
}

void CoefFunctionApproxDeriv::GetTensor(Matrix<Double>& coefMat, 
                                        const LocPointMapped& lpm ) {


  // evaluate vector of dependency
  Vector<Double> elemB;
  dependCoef_->GetVector( elemB, lpm);
  Double fieldAbs = elemB.NormL2();
  coefMat.Resize( dimDMat_, dimDMat_ );
  if( fieldAbs == 0 ) {
    coefMat.Init();
  } else {

    // evaluate derivative of reluctivity
    Double nuPrime = nLinFnc_->EvaluatePrimeNu(fieldAbs);
    // coefMat = B^T [ e_B^T * nu' * |B| * e_B] B 
    Vector<Double> eB(dimDMat_);
    eB = elemB / fieldAbs;
    coefMat.DyadicMult( eB );
    coefMat *= nuPrime * fieldAbs;
  }
}


std::string CoefFunctionApproxDeriv::ToString() const {
  EXCEPTION( "Implement me");
  return "";
}

// ============================================================================
//  ANISOTROPIC VERSIONS
// ============================================================================

CoefFunctionApproxAniso::CoefFunctionApproxAniso() : CoefFunction() {
  // this type of coefficient is nonlinear (i.e. solution dependend)
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_ = false;
}

//! Destructor
CoefFunctionApproxAniso::~CoefFunctionApproxAniso(){
  ;
}

//! Initialize with data
void CoefFunctionApproxAniso::Init( Double coefScalar, 
                                    StdVector<ApproxData*>  nLinFnc,
                                    StdVector<Double> angles,
                                    PtrCoefFct dependCoef ) {

  // set type to scalar
  dimType_ = SCALAR;
  coefScalar_ = coefScalar;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  dependCoef_ = dependCoef;
}

void CoefFunctionApproxAniso::GetScalar(Double& coefScalar, 
                                        const LocPointMapped& lpm ) {

  // evaluate vector of dependency
  Vector<Double> elemSol;
  dependCoef_->GetVector( elemSol, lpm);

  Double fieldAbs = elemSol.NormL2();

  if( fieldAbs == 0 ) { 
    coefScalar = coefScalar_;
  } else {


    //compute angle 
    Double angleB;
    if ( abs(elemSol[0]) > 1e-5 ) {
      angleB = abs( std::atan( elemSol[1] / elemSol[0] ) );
      angleB *= 180.0/3.1416;
    }
    else {
      angleB = 90.0;
    }
    

    // -----------------------------
    //  Angular based interpolation
    // -----------------------------
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleB > angles_[kend] || kend == 0) {
      coefScalar = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs);
    }
    else if ( angleB < angles_[0] ) {
      coefScalar = nLinFnc_[0]->EvaluateFuncNu(fieldAbs);
    }
    else {
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of xEntry
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleB)
          khi=k;
        else
          klo=k;
      }

      // size of x interval
//      std::cerr << "ANGLE: " << angleB;
//      std::cerr << "startAngle: " << angles_[klo]; 
//      std::cerr << ", stopAngle: " << angles_[khi];
      Double dxVal=angles_[khi] - angles_[klo];

      // The x-values must be distinct!
      if (dxVal == 0.0) {
        EXCEPTION("You cannot have two equal angular values!" );
      }

      // relative distance of xEntry to x-Value bounds
      Double a = ( angles_[khi] - angleB )/dxVal;
      Double b = ( angleB - angles_[klo] )/dxVal;
//      std::cerr << ", a = " << a << ", b=" << b << std::endl;
      coefScalar =   a * nLinFnc_[klo]->EvaluateFuncNu(fieldAbs) 
                   + b * nLinFnc_[khi]->EvaluateFuncNu(fieldAbs);
    }

    // ---------------------------------
    //  Nearest neighbour interpolation
    // ---------------------------------
    
//    UInt pos = 0;
//    Double dist, minDist;
//    minDist = abs( angles_[0] - angleB );
//    for (UInt i=1; i<angles_.GetSize(); i++ ) {
//      dist = abs( angles_[i] - angleB );
//      if ( dist < minDist ) {
//        pos = i;
//        minDist = dist;
//      }
//    }
//    std::cerr << "elem: " << lpm.ptEl->elemNum << "absB: " << fieldAbs << ", angle: " << angleB << ", pos: " << pos << std::endl;
//    coefScalar = nLinFnc_[pos]->EvaluateFuncNu(fieldAbs);

  }
}

std::string CoefFunctionApproxAniso::ToString() const {
  EXCEPTION( "Implement me");
  return "";
}

// ==========================================================================
CoefFunctionApproxDerivAniso::CoefFunctionApproxDerivAniso() : CoefFunction() {
  // this type of coefficient is nonlinear (i.e. solution dependend)
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_ = false;
}

//! Destructor
CoefFunctionApproxDerivAniso::~CoefFunctionApproxDerivAniso(){
  ;
}

//! Initialize with data
void CoefFunctionApproxDerivAniso::Init( StdVector<ApproxData*>  nLinFnc,
                                         StdVector<Double> angles,
                                         UInt dimDMat,
                                         PtrCoefFct dependCoef ) {
  // set type to TENSOR
  dimType_ = TENSOR;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  dimDMat_ = dimDMat;
  dependCoef_ = dependCoef;
}

void CoefFunctionApproxDerivAniso::GetTensor(Matrix<Double>& coefMat, 
                                             const LocPointMapped& lpm ) {

  // evaluate vector of dependency
  Vector<Double> elemB;
  dependCoef_->GetVector( elemB, lpm);

  Double fieldAbs = elemB.NormL2();
  coefMat.Resize( dimDMat_, dimDMat_ );
  if( fieldAbs == 0 ) {
    coefMat.Init();
  } else {

    //get pointers to the nlFunctions

    //compute angle 
    Double angleB;
    if ( abs(elemB[0]) > 1e-5 ) {
      angleB = abs( std::atan( elemB[1] / elemB[0] ) );
      angleB *= 180.0/3.1416;
    }
    else {
      angleB = 90.0;
    }

    Double nuPrime = 0.0;
    // -----------------------------
    //  Angular based interpolation
    // -----------------------------
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleB > angles_[kend] || kend == 0) {
      nuPrime = nLinFnc_[kend]->EvaluatePrimeNu(fieldAbs);
    }
    else if ( angleB < angles_[0] ) {
      nuPrime = nLinFnc_[0]->EvaluatePrimeNu(fieldAbs);
    }
    else {
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of xEntry
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleB)
          khi=k;
        else
          klo=k;
      }

//      std::cerr << "ANGLE: " << angleB;
//      std::cerr << "startAngle: " << angles_[klo]; 
//      std::cerr << ", stopAngle: " << angles_[khi];
      Double dxVal=angles_[khi] - angles_[klo];

      // The x-values must be distinct!
      if (dxVal == 0.0) {
        EXCEPTION("You cannot have two equal x values!" );
      }

      // relative distance of xEntry to x-Value bounds
      Double a = ( angles_[khi] - angleB )/dxVal;
      Double b = ( angleB - angles_[klo] )/dxVal;
//      std::cerr << ", a = " << a << ", b=" << b << std::endl;
      nuPrime =   a * nLinFnc_[klo]->EvaluatePrimeNu(fieldAbs) 
                + b * nLinFnc_[khi]->EvaluatePrimeNu(fieldAbs);
    }
    
    // ---------------------------------
    //  Nearest neighbour interpolation
    // ---------------------------------
//    UInt pos = 0;
//       Double dist, minDist;
//       minDist = abs( angles_[0] - angleB );
//       for (UInt i=1; i<angles_.GetSize(); i++ ) {
//         dist = abs( angles_[i] - angleB );
//         if ( dist < minDist ) {
//           pos = i;
//           minDist = dist;
//         }
//       }
//    // evaluate derivative of reluctivity
//    Double nuPrime = nLinFnc_[pos]->EvaluatePrimeNu(fieldAbs);
    
    
    
    
    // coefMat = B^T [ e_B^T * nu' * |B| * e_B] B 
    Vector<Double> eB(dimDMat_);
    eB = elemB / fieldAbs;
    coefMat.DyadicMult( eB );
    coefMat *= nuPrime * fieldAbs;
  }
}

std::string CoefFunctionApproxDerivAniso::ToString() const {
  EXCEPTION( "Implement me");
  return "";
}

} // namespace
