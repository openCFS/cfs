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

Double CoefFunctionComposite::GetTerminalValue(TerminalConnector tc, NonLinType nl, const LocPointMapped & lpm) {
  if (surfElems_.find(tc) != surfElems_.end())
    return GetLocalTerminalValue(tc, nl, lpm);
  else
    return GetAvgTerminalValue(tc, nl, lpm);
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

Double CoefFunctionComposite::GetMaxTerminalValue(TerminalConnector tc, NonLinType nl, const LocPointMapped & lpm) {

  EXCEPTION("Easy to implement");
  Double max = 0.;
  return max;
}

void CoefFunctionComposite::SetLocValue(TerminalConnector tc) {
  if ( surfElems_.find(tc) != surfElems_.end() ) {
    EXCEPTION("Terminal " << tc << " was already defined.");
  }
  surfElems_[tc][0] = NULL; // force seg fault at access of element 0, since element number start from 1, this should be ok
  }

Double CoefFunctionComposite::GetLocalTerminalValue(TerminalConnector tc, NonLinType nl, const LocPointMapped & lpm) {
  Double val = 0.;
  UInt volElemNum = lpm.ptEl->elemNum;
  
  if ( surfElems_[tc].find(volElemNum) == surfElems_[tc].end() ) {
    // element not yet found
    Grid * ptGrid = lpm.GetShapeMap()->GetGrid();
    StdVector <Elem*> surfElem;
    ptGrid->GetAdjacentSurfElem(volElemNum, surfElem, terminals_[tc]);
    if (surfElem.GetSize() != 1) {
      EXCEPTION("Could not find unique adjacent surface element for volume element '" << volElemNum
      << "' in region '" << terminals_[tc] << "'. What I got was: " << surfElem.ToString());
    }
    // add to map
    surfElems_[tc].insert(std::make_pair(volElemNum, surfElem[0]));
    LOG_DBG(coeffctapprox) << "Volume element '" << volElemNum << "' has adjacent surface element: " << surfElem[0]->elemNum;
  }

  dependCoefs_[nl]->GetAvgElemValue(val, surfElems_[tc][volElemNum]);
  return val;
}

void CoefFunctionComposite::MultiplyByElemArea( Double & value, const LocPointMapped & lpm) {
  if (!multElemArea_)
    return;

  UInt volElemNum = lpm.ptEl->elemNum;

  if ( elemAreas_.find(volElemNum) == elemAreas_.end() ) {
    // element not yet calculated
    Grid * ptGrid = lpm.GetShapeMap()->GetGrid();
    StdVector <Elem*> surfElem;
    ptGrid->GetAdjacentSurfElem(volElemNum, surfElem, terminals_[tcElemArea_]);
    if (surfElem.GetSize() != 1) {
      EXCEPTION("Could not find unique adjacent surface element for volume element '" << volElemNum
      << "' in region '" << terminals_[tcElemArea_] << "'. What I got was: " << surfElem.ToString());
    }
    // calculate area of this element
    shared_ptr<ElemShapeMap> esm = ptGrid->GetElemShapeMap( surfElem[0]);
    Double area = esm->CalcVolume();
    // add to map
    elemAreas_.insert(std::make_pair(volElemNum, area));
    LOG_DBG(coeffctapprox) << "Volume element '" << volElemNum << "' has adjacent surface element '" << surfElem[0]->elemNum
      << "', which has an area of: " << area;
  } 

  value *= elemAreas_[volElemNum];

}

void CoefFunctionComposite::SetElemAreaMult(TerminalConnector tc) {
  multElemArea_ = true;
  tcElemArea_ = tc;
}
void CoefFunctionComposite::Init( Double coefScalar, ApproxData * nLinFnc) {

  dimType_ = SCALAR;
  nLinFnc_ = nLinFnc;
  coefScalar_ = coefScalar;
  divideByVds_ = false;
  multElemArea_ = false;

}

// ========================================================================
//
// ============================================================================
//  Coef Function Bipole
// ============================================================================
//
void CoefFunctionBipole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {

  Double aAvg = GetTerminalValue(ANODE, NLELEC_BIPOLE, lpm);
  Double bAvg = GetTerminalValue(CATHODE, NLELEC_BIPOLE, lpm);
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
  Double aAvg = GetTerminalValue(ANODE, NLELEC_BIPOLE, lpm);
  Double bAvg = GetTerminalValue(CATHODE, NLELEC_BIPOLE, lpm);
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

  Double Vdrain = GetTerminalValue(DRAIN, NLELEC_TRIPOLE, lpm);
  Double Vsource = GetTerminalValue(SOURCE, NLELEC_TRIPOLE, lpm);
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
CoefFunctionHeatTripole::CoefFunctionHeatTripole() : CoefFunctionComposite() {
}

void CoefFunctionHeatTripole::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {
  Double Vdrain = GetTerminalValue(DRAIN, NLELEC_TRIPOLE, lpm);
  Double Vsource = GetTerminalValue(SOURCE, NLELEC_TRIPOLE, lpm);
  Double Vgate = GetAvgTerminalValue(GATE, NLELEC_TRIPOLE, lpm);
  Double Vds = Vdrain-Vsource;
  Double Vgs = Vgate-Vsource;

  // get temperature
  Vector<Double> elemSol;
  dependCoefs_[NLELEC_CONDUCTIVITY]->GetVector( elemSol, lpm);

  // special case: my material actually returns current and thus has to be divided by Vds 
  if (divideByVds_) {
    if (fabs(Vds) < 1e-9) { // Volts
      coefScalar = 1e-4;
    } else {
      coefScalar = 1./Vds*nLinFnc_->EvaluateFunc(Vds, Vgs, elemSol[0]);
    }
  }
  else {
    coefScalar = nLinFnc_->EvaluateFunc(Vds, Vgs, elemSol[0]);
  }

  MultiplyByElemArea(coefScalar, lpm);

  LOG_DBG(coeffctapprox) << "Returning approximated scalar '" << coefScalar << "' for dependVal vgs = '" << Vgs << " and vds = '" << Vds <<"' and temperature = '" << elemSol[0] <<"'. IP '" << lpm.lp.number << "', '" << lpm.lp.coord.ToString() << "' in element :" << lpm.ptEl->elemNum;

}

//
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

//! \see CoefFunction::GetScalar
void CoefFunctionApproxDeriv::GetScalar(Double& coefScalar,
                                       const LocPointMapped& lpm ) {


  // evaluate vector of dependency
  Vector<Double> elemB;
  dependCoef_->GetVector( elemB, lpm);
  Double fieldAbs = elemB.NormL2();

  if( fieldAbs == 0 ) {
    coefScalar = 0.0;
  }
  else {
    // evaluate derivative of reluctivity
    coefScalar = nLinFnc_->EvaluatePrimeNu(fieldAbs);
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
                                    StdVector<shared_ptr<CoefFunction> > nLinFnc,
                                    StdVector<Double> angles,
                                    StdVector<Double> zScalings,
                                    PtrCoefFct dependCoef ) {

  // set type to scalar
  dimType_ = SCALAR;
  coefScalar_ = coefScalar;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  zScalings_ = zScalings;
  dependCoef_ = dependCoef;
}

void CoefFunctionApproxAniso::GetScalar(Double& coefScalar, 
                                        const LocPointMapped& lpm ) {

  // evaluate vector of dependency
  Vector<Double> elemSol;
  dependCoef_->GetVector( elemSol, lpm);

  Double fieldAbs = elemSol.NormL2();

  if( fieldAbs < 1.0e-04 ) {
    coefScalar = coefScalar_;
  } else {

    // x-values have to be positive! 
//    if ( elemSol[0] <= 0 ) { 
//      EXCEPTION("CoefFunctionApproxAniso::GetScalar(): x value has to be positive!" );
//    }

    // -------------------------------
    //  compute angle phi of B vector 
    // -------------------------------
    Double angleBPhi;
    if ( abs(elemSol[0]) > 1e-5 ) { 
      angleBPhi = abs( std::atan( elemSol[1] / elemSol[0] ) );
      angleBPhi *= 180.0/3.141592654; // conversion rad to deg
    }
    else {
      angleBPhi = 90.0;
    }
    
    // ---------------------------------
    //  compute angle theta of B vector
    // ---------------------------------
    Double angleBTheta;
    angleBTheta = std::acos( elemSol[2] / sqrt(elemSol[0]*elemSol[0] + 
        elemSol[1]*elemSol[1] + elemSol[2]*elemSol[2] ));
    angleBTheta *= 180.0/3.141592654; // conversion rad to deg
    
    // theta in spherical coordinates is defined as the angle between the 
    // zenith direction and the line segment (origin->point) with a range [0°;180°]
    // therefore change range to [90°;-90°] and only use absolute value of it (-> symmetry!)
    angleBTheta = abs(90 - angleBTheta);   

    // take care that theta does not exceed its limits 
    // TODO-avolk: Improve handling like it is done for phi!
    if (angleBTheta > 90) {
      angleBTheta = 90;
      WARN("GetScalar(): theta of element " << lpm.ptEl->elemNum << " was greater than 90 degrees! :-(")
    }
    
    
    // ---------------------------- 
    //  Nearest neighbour approach
    // ----------------------------
//    UInt pos = 0;
//    Double dist, minDist;
//    minDist = abs( angles_[0] - angleBPhi );
//    for (UInt i=1; i<angles_.GetSize(); i++ ) {
//      dist = abs( angles_[i] - angleBPhi );
//      if ( dist < minDist ) {
//        pos = i;
//        minDist = dist;
//      }
//    }
//    std::cerr << "elem: " << lpm.ptEl->elemNum << "absB: " << fieldAbs << ", angle: " << angleBPhi << ", pos: " << pos << std::endl;
//    coefScalar = nLinFnc_[pos]->EvaluateFuncNu(fieldAbs);

    
    // ------------------------------------   
    // Interpolation method of coefficients
    // ------------------------------------
    //
    //                90°
    //                 o
    //                / \
    //               /   \
    //            4 o     o
    //     theta   /       \
    //            /         \
    //         3 o           o
    //          /     x       \
    //         /               \
    //     0° o-----o-----o-----o 90°
    //              1     2
    //                phi
    //
    // - locally interpolate over phi (points 1 and 2) to get coef_xy
    // - locally interpolate over theta (points 3 and 4) to get coef_z
    // - interpolate coef_xy and coef_z over theta to get final coefficient
    
    // ---------------------------------------------------
    //  Angular based interpolation in the xy-plane (phi)
    // ---------------------------------------------------
    
    Double coefScalarXY = 0.0;
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleBPhi > angles_[kend] || kend == 0) {
      //coefScalarXY = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs);
      nLinFnc_[kend]->GetScalar(coefScalarXY, lpm );
    }
    else if ( angleBPhi < angles_[0] ) {
      //coefScalarXY = nLinFnc_[0]->EvaluateFuncNu(fieldAbs);
      nLinFnc_[0]->GetScalar(coefScalarXY, lpm);
    }
    else {  // now that the angle is in range, search for the "best" given angle 
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      // klo and khi bracket the input value of phi-value
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleBPhi)
          khi=k;
        else
          klo=k;
      }
      
      // size of phi interval
      Double dPhiVal = angles_[khi] - angles_[klo]; 

      // The phi-values must be distinct!
      if (dPhiVal == 0.0) {
        EXCEPTION("You cannot have two equal angular values!" );
      }

      // relative distance of phi-value to phi-bounds
      Double ahi = ( angles_[khi] - angleBPhi ) / dPhiVal;
      Double alo = ( angleBPhi - angles_[klo] ) / dPhiVal;
      
      // value of coefficient interpolated within the xy-plane and theta=0°
//      coefScalarXY =   ahi * nLinFnc_[klo]->EvaluateFuncNu(fieldAbs)
//                      + alo * nLinFnc_[khi]->EvaluateFuncNu(fieldAbs);
      Double VALklo, VALkhi;
      nLinFnc_[klo]->GetScalar(VALklo, lpm);
      nLinFnc_[khi]->GetScalar(VALkhi, lpm);
      coefScalarXY =   ahi * VALklo + alo * VALkhi;
    }
       
    // --------------------------------------------------------------------
    //  Angular based interpolation according to z direction (angle theta)
    // --------------------------------------------------------------------
    
    // ASSUMPTION: material behavior of in thickness direction (theta) is 
    //             equivalent to behavior in transverse direction (phi).
    //             therefore we use the same BH-curves but scale them by 
    //             an appropriate scaling factor. 
    
    Double coefScalarZ = 0.0;

    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e. first or last) 
    if ( angleBTheta > angles_[kend] || kend == 0) {
      //coefScalarZ = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs);
      nLinFnc_[kend]->GetScalar(coefScalarZ, lpm);
    }
    else if ( angleBTheta < angles_[0] ) {
//      coefScalarZ = nLinFnc_[0]->EvaluateFuncNu(fieldAbs);
      nLinFnc_[0]->GetScalar(coefScalarZ, lpm);
    }
    else {  // now that the angle is in range, search for the "best" given angle 
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      // klo and khi bracket the input value of theta-value
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleBTheta)
          khi=k;
        else
          klo=k;
      }

      // size of theta interval
      Double dThetaVal = angles_[khi] - angles_[klo]; 

      // The theta-values must be distinct!
      if (dThetaVal == 0.0) {
        EXCEPTION("You cannot have two equal angular values!" );
      }

      // relative distance of theta-value to theta-bounds
      Double ahi = ( angles_[khi] - angleBTheta ) / dThetaVal;
      Double alo = ( angleBTheta - angles_[klo] ) / dThetaVal;

      // value of coefficient interpolated along z-direction and phi=0°
      // Note: scaling of mu by factor c leads to scaling of nu by 1/c since: 
      //       nu = 1/mu; c*mu => 1/c*mu
      Double VALkhi, VALklo;
      nLinFnc_[klo]->GetScalar(VALklo, lpm);
      nLinFnc_[khi]->GetScalar(VALkhi, lpm);
      coefScalarZ =   ahi * VALklo * (1/zScalings_[klo])
                     + alo * VALkhi * (1/zScalings_[khi]);
//      coefScalarZ =   ahi * nLinFnc_[klo]->EvaluateFuncNu(fieldAbs) * (1/zScalings_[klo])
//                    + alo * nLinFnc_[khi]->EvaluateFuncNu(fieldAbs) * (1/zScalings_[khi]);
    }
    
    // ----------------------------------------------------------------------------------------
    //  Interpolation between calculated coefficients in xy-plane and z-direction (angle theta)
    // ----------------------------------------------------------------------------------------
    
    Double bhi = ( 90 - angleBTheta ) / 90;
    Double blo = ( angleBTheta - 0.0 ) / 90;
    coefScalar = bhi * coefScalarXY + blo * coefScalarZ;


#ifndef NDEBUG
    // some debug output... TODO-avolk: weg damit!
    std::cerr << "GetScalar(): statistics of element " << lpm.ptEl->elemNum << ":";
//    std::cerr << " B_x=" << elemSol[0];
//    std::cerr << " B_y=" << elemSol[1];
//    std::cerr << " B_z=" << elemSol[2];
    std::cerr << " angleBPhi=" << angleBPhi;
    std::cerr << " angleBTheta=" << angleBTheta;
    std::cerr << " zScaling(end)=" << zScalings_[kend];
//    std::cerr << " startAngle=" << angles_[klo]; 
//    std::cerr << " stopAngle=" << angles_[khi];
//    std::cerr << " ahi=" << ahi << ", alo=" << alo << ", ahi+alo=" << ahi+alo;
    std::cerr << " coefScalarXY=" << coefScalarXY;
    std::cerr << " coefScalarZ=" << coefScalarZ;
    std::cerr << " coefScalar=" << coefScalar;
//    std::cerr << std::setprecision(4); 
    std::cerr << std::endl;
#endif
    
  }
}

std::string CoefFunctionApproxAniso::ToString() const {
  EXCEPTION( "Implement me");
  return "";
}

// ==========================================================================
CoefFunctionApproxDerivAniso::CoefFunctionApproxDerivAniso() : CoefFunction() {
  // this type of coefficient is nonlinear (i.e. solution dependent)
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_ = false;
}

//! Destructor
CoefFunctionApproxDerivAniso::~CoefFunctionApproxDerivAniso(){
  ;
}

//! Initialize with data
void CoefFunctionApproxDerivAniso::Init( StdVector<shared_ptr<CoefFunction> > nLinFnc,
                                         StdVector<Double> angles,
                                         StdVector<Double> zScalings,
                                         UInt dimDMat,
                                         PtrCoefFct dependCoef ) {
  // set type to TENSOR
  dimType_ = TENSOR;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  zScalings_ = zScalings;
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
  if( fieldAbs == 0.0 ) {
    coefMat.Init();
  } else {

    // -------------------------------
    //  compute angle phi of B vector
    // -------------------------------
    Double angleBPhi;
    if ( abs(elemB[0]) > 1e-5 ) {
      angleBPhi = abs( std::atan( elemB[1] / elemB[0] ) );
      angleBPhi *= 180.0/3.141592654; // conversion rad to deg
    }
    else {
      angleBPhi = 90.0;
    }

    // ---------------------------------
    //  compute angle theta of B vector
    // ---------------------------------
    Double angleBTheta;
    angleBTheta = std::acos( elemB[2] / sqrt(elemB[0]*elemB[0] + 
                              elemB[1]*elemB[1] + elemB[2]*elemB[2] ));
    angleBTheta *= 180.0/3.141592654; // conversion rad to deg

    // theta in spherical coordinates is defined as the angle between the 
    // zenith direction and the line segment (origin->point) with a range [0°;180°]
    // therefore change range to [90°;-90°] and only use absolute value of it (-> symmetry!)
    angleBTheta = abs(90 - angleBTheta);

    // take care that theta does not exceed its limits 
    // TODO-avolk: Improve handling like it is done for phi!
    if (angleBTheta > 90) {
      angleBTheta = 90;
      WARN("GetTensor(): theta of element " << lpm.ptEl->elemNum << " was greater than 90 degrees! :-(")
    }
    
    
    // ---------------------------- 
    //  Nearest neighbour approach
    // ----------------------------
//    UInt pos = 0;
//    Double dist, minDist;
//    minDist = abCoefFunctionApproxAniso::GetScalars( angles_[0] - angleB );
//    for (UInt i=1; i<angles_.GetSize(); i++ ) {
//      dist = abs( angles_[i] - angleB );
//      if ( dist < minDist ) {
//        pos = i;
//        minDist = dist;
//      }
//    }
//    // evaluate derivative of reluctivity
//    Double nuPrime = nLinFnc_[pos]->EvaluatePrimeNu(fieldAbs);
    
    
    // ------------------------------------   
    // Interpolation method of coefficients
    // ------------------------------------
    //
    //                90°
    //                 o
    //                / \
    //               /   \
    //            4 o     o
    //     theta   /       \
    //            /         \
    //         3 o           o
    //          /     x       \
    //         /               \
    //     0° o-----o-----o-----o 90°
    //              1     2
    //                phi
    //
    // - locally interpolate over phi (points 1 and 2) to get coef_xy
    // - locally interpolate over theta (points 3 and 4) to get coef_z
    // - interpolate coef_xy and coef_z over theta to get final coefficient    
    
    // ---------------------------------------------------
    //  Angular based interpolation in the xy-plane (phi)
    // ---------------------------------------------------
    Double nuPrime, nuPrimeXY, nuPrimeZ;
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleBPhi > angles_[kend] || kend == 0) {
      //nuPrimeXY = nLinFnc_[kend]->EvaluatePrimeNu(fieldAbs);
      nLinFnc_[kend]->GetScalar(nuPrimeXY, lpm);
    }
    else if ( angleBPhi < angles_[0] ) {
      //nuPrimeXY = nLinFnc_[0]->EvaluatePrimeNu(fieldAbs);
      nLinFnc_[0]->GetScalar(nuPrimeXY, lpm);
    }
    else {
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of phi-value
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleBPhi)
          khi=k;
        else
          klo=k;
      }
      
      // size of phi interval
      Double dPhiVal = angles_[khi] - angles_[klo];

      // The phi-values must be distinct!
      if (dPhiVal == 0.0) {
        EXCEPTION("You cannot have two equal angular values!" );
      }

      // relative distance of phi-value to phi-bounds
      Double ahi = ( angles_[khi] - angleBPhi ) / dPhiVal;
      Double alo = ( angleBPhi - angles_[klo] ) / dPhiVal;

      // value of nuPrime interpolated within the xy-plane and theta=0°
      Double VALklo, VALkhi;
      nLinFnc_[klo]->GetScalar(VALklo, lpm);
      nLinFnc_[khi]->GetScalar(VALkhi, lpm);
      nuPrimeXY =   ahi * VALklo + alo * VALkhi;
//      nuPrimeXY =   ahi * nLinFnc_[klo]->EvaluatePrimeNu(fieldAbs)
//                  + alo * nLinFnc_[khi]->EvaluatePrimeNu(fieldAbs);
    }
    
    // --------------------------------------------------------------------
    //  Angular based interpolation according to z direction (angle theta)
    // --------------------------------------------------------------------
        
    // ASSUMPTION: material behavior of in thickness direction (theta) is 
    //             equivalent to behavior in transverse direction (phi).
    //             therefore we use the same BH-curves but scale them by 
    //             an appropriate scaling factor. 
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    if ( angleBTheta > angles_[kend] || kend == 0) {
      nLinFnc_[kend]->GetScalar(nuPrimeZ, lpm);
      //nuPrimeZ = nLinFnc_[kend]->EvaluatePrimeNu(fieldAbs);
    }
    else if ( angleBTheta < angles_[0] ) {
      nLinFnc_[0]->GetScalar(nuPrimeZ, lpm);
//      nuPrimeZ = nLinFnc_[0]->EvaluatePrimeNu(fieldAbs);
    }
    else {
      UInt klo,khi,k;
      klo=0;
      khi=kend;
      // We will find the right place in the table by means of bisection.
      //  klo and khi bracket the input value of theta-value
      while (khi-klo > 1) {
        k=(khi+klo) >> 1; // binary right shift
        if (angles_[k] > angleBTheta)
          khi=k;
        else
          klo=k;
      }

      // size of theta interval
      Double dThetaVal = angles_[khi] - angles_[klo];

      // The theta-values must be distinct!
      if (dThetaVal == 0.0) {
        EXCEPTION("You cannot have two equal angular values!" );
      }

      // relative distance of phi-value to phi-bounds
      Double ahi = ( angles_[khi] - angleBTheta ) / dThetaVal;
      Double alo = ( angleBTheta - angles_[klo] ) / dThetaVal;

      // value of nuPrime interpolated along z-direction and phi=0°
      // Note: scaling of mu by factor c leads to scaling of nu' by 1/c since:
      // 1) nu = 1/mu; c*mu => 1/c*mu meaning that scaling has to be applied by the reciprocal factor 
      // 2) g(x)=c*f(x) => g'(x)=c*f'(x) meaning that scaling can also be applied to derivative
      Double VALklo, VALkhi;
      nLinFnc_[klo]->GetScalar(VALklo, lpm);
      nLinFnc_[khi]->GetScalar(VALkhi, lpm);
      nuPrimeZ =   ahi * VALklo * (1/zScalings_[klo])
                 + alo * VALkhi * (1/zScalings_[khi]);
//      nuPrimeZ =   ahi * nLinFnc_[klo]->EvaluatePrimeNu(fieldAbs) * (1/zScalings_[klo])
//                 + alo * nLinFnc_[khi]->EvaluatePrimeNu(fieldAbs) * (1/zScalings_[khi]);
    }

    // -----------------------------------------------------------------------------------
    //  Interpolation between calculated nuPrime in xy-plane and z-direction (angle theta)
    // -----------------------------------------------------------------------------------
    
    Double bhi = ( 90 - angleBTheta ) / 90;
    Double blo = ( angleBTheta - 0.0 ) / 90;
    nuPrime = bhi * nuPrimeXY + blo * nuPrimeZ;        
    
    
    // coefMat = B^T [ e_B^T * nu' * |B| * e_B] B 
    Vector<Double> eB(dimDMat_);
    eB = elemB / fieldAbs;
    coefMat.DyadicMult( eB );
    coefMat *= nuPrime * fieldAbs;
    
    
#ifndef NDEBUG
    // some debug output... TODO-avolk: weg damit!
    std::cerr << "GetTensor(): statistics of element " << lpm.ptEl->elemNum;
    std::cerr << ": angleBPhi=" << angleBPhi;
    std::cerr << ", angleBTheta=" << angleBTheta;
    std::cerr << ", zScaling(end)=" << zScalings_[kend];
//    std::cerr << ", startAngle=" << angles_[klo]; 
//    std::cerr << ", stopAngle=" << angles_[khi];
//    std::cerr << ", ahi=" << ahi << ", alo=" << alo << ", ahi+alo=" << ahi+alo;
    std::cerr << ", nuPrimeXY=" << nuPrimeXY;
    std::cerr << ", nuPrimeZ=" << nuPrimeZ;
    std::cerr << ", nuPrime=" << nuPrime;
//    std::cerr << std::setprecision(4); 
    std::cerr << std::endl;
#endif
    
  }
}

std::string CoefFunctionApproxDerivAniso::ToString() const {
  EXCEPTION( "Implement me");
  return "";
}

} // namespace
