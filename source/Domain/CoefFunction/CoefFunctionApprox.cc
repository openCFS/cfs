#include "CoefFunctionApprox.hh"

// classes for function / spline approximation
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"


namespace CoupledField {

// ============================================================================
//  ISOTROPIC VERSIONS
// ============================================================================

CoefFunctionApprox::CoefFunctionApprox() : CoefFunction() {
  // this type of coefficient is nonlinear (i.e. solution dependend)
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  bOperator_ = NULL;
  isComplex_ = false;
}

CoefFunctionApprox::~CoefFunctionApprox(){
  delete bOperator_;
  bOperator_ = NULL;
  ;
}

void CoefFunctionApprox::Init( Double coefScalar, ApproxData * nLinFnc,
                               shared_ptr<FeFunction<Double> > fct,
                               BaseBOperator* bOp ) {

  // set type to scalar
  dimType_ = SCALAR;
  coefScalar_ = coefScalar;
  nLinFnc_ = nLinFnc;
  feFct_ = fct;
  bOperator_ = bOp;
}

//! \see CoefFunction::GetScalar
void CoefFunctionApprox::GetScalar(Double& coefScalar, 
                                   const LocPointMapped& lpm ) {

  // extract element solution from feFunction
  Vector<Double> elemSol, elemOpSol;
  feFct_->GetElemSolution(elemSol, lpm.ptEl);

  // apply b-operator matrix to element solution to obtain field value
  BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
  bOperator_->ApplyOp( elemOpSol, lpm, ptFe, elemSol );

  if ( nLinFnc_->GetMatType() == MAG_PERMEABILITY ) {
    Double fieldAbs = elemOpSol.NormL2();

    if( fieldAbs == 0 ) { 
      coefScalar = coefScalar_;
    } else {
      coefScalar = nLinFnc_->EvaluateFuncNu(fieldAbs);
    }
  }
  else {
    coefScalar = nLinFnc_->EvaluateFunc(elemOpSol[0]);
  }
}

bool IsComplex(){
  return false;
}

std::string CoefFunctionApprox::ToString() const {
  return "";
}


// ========================================================================

CoefFunctionApproxDeriv::CoefFunctionApproxDeriv() : CoefFunction() {
  // this type of coefficient is nonlinear, i.e. spatial and time dependent
  dependType_ = GENERAL;
  isComplex_ = false;
}

CoefFunctionApproxDeriv::~CoefFunctionApproxDeriv(){
  delete bOperator_;
  bOperator_ = NULL;
}

void CoefFunctionApproxDeriv::Init( ApproxData * nLinFnc,
                                   shared_ptr<FeFunction<Double> > fct,
                                   BaseBOperator* bOp ) {

  // set type to TENSOR
  dimType_ = TENSOR;
  nLinFnc_ = nLinFnc;
  feFct_ = fct;
  bOperator_ = bOp;
  dimDMat_ = bOperator_->GetDimDMat();
}

void CoefFunctionApproxDeriv::GetTensor(Matrix<Double>& coefMat, 
                                        const LocPointMapped& lpm ) {


  // extract element solution from feFunction
  Vector<Double> elemSol, elemB;
  feFct_->GetElemSolution(elemSol, lpm.ptEl);

  // apply b-operator matrix to element solution to obtain field value
  BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
  bOperator_->ApplyOp( elemB, lpm, ptFe, elemSol );

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
  bOperator_ = NULL;
}

//! Destructor
CoefFunctionApproxAniso::~CoefFunctionApproxAniso(){
  delete bOperator_;
  bOperator_ = NULL;
  ;
}

//! Initialize with data
void CoefFunctionApproxAniso::Init( Double coefScalar, 
                                    StdVector<ApproxData*>  nLinFnc,
                                    StdVector<Double> angles,
                                    shared_ptr<FeFunction<Double> > fct,
                                    BaseBOperator* bOp ) {

  // set type to scalar
  dimType_ = SCALAR;
  coefScalar_ = coefScalar;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  feFct_ = fct;
  bOperator_ = bOp;
  
  // TODO-avolk: this should not be hardcoded here but specified in the xml file!!
  zScaling_ = ZSCALING;
}

void CoefFunctionApproxAniso::GetScalar(Double& coefScalar, 
                                        const LocPointMapped& lpm ) {

  // extract element solution from feFunction
  Vector<Double> elemSol, elemOpSol;
  feFct_->GetElemSolution(elemSol, lpm.ptEl);

  // apply b-operator matrix to element solution to obtain field value
  BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
  bOperator_->ApplyOp( elemOpSol, lpm, ptFe, elemSol );

  Double fieldAbs = elemOpSol.NormL2();

  if( fieldAbs == 0 ) { 
    coefScalar = coefScalar_;
  } else {

    // x-values have to be positive! TODO-avolk: or is this handled somewhere below??
//    if ( elemOpSol[0] <= 0 ) { 
//      EXCEPTION("CoefFunctionApproxAniso::GetScalar(): x value has to be positive!" );
//    }

    // -------------------------------
    //  compute angle phi of B vector 
    // -------------------------------
    Double angleBPhi;
    if ( abs(elemOpSol[0]) > 1e-5 ) { // why is this done?? TODO-avolk
      angleBPhi = abs( std::atan( elemOpSol[1] / elemOpSol[0] ) );
      angleBPhi *= 180.0/3.141592654; // conversion rad to deg
    }
    else {
      angleBPhi = 90.0;
    }
    
    // ---------------------------------
    //  compute angle theta of B vector
    // ---------------------------------
    Double angleBTheta;
    angleBTheta = std::acos( elemOpSol[2] / sqrt(elemOpSol[0]*elemOpSol[0] + 
        elemOpSol[1]*elemOpSol[1] + elemOpSol[2]*elemOpSol[2] ));
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
    
    
    // ---------------------------------
    //  Nearest neighbour interpolation
    // ---------------------------------
    
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
    
    
    // ---------------------------------------------------
    //  Angular based interpolation in the xy-plane (phi)
    // ---------------------------------------------------
    
    Double coefScalarXY = 0.0;
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleBPhi > angles_[kend] || kend == 0) {
      coefScalarXY = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs);
    }
    else if ( angleBPhi < angles_[0] ) {
      coefScalarXY = nLinFnc_[0]->EvaluateFuncNu(fieldAbs);
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
      coefScalarXY =   ahi * nLinFnc_[klo]->EvaluateFuncNu(fieldAbs) 
                     + alo * nLinFnc_[khi]->EvaluateFuncNu(fieldAbs);      
    }
    
    // --------------------------------------------------------------------
    //  Angular based interpolation according to z direction (angle theta)
    // --------------------------------------------------------------------

    // TODO-avolk: this should be deleted after first tests!!
    // added just to be shure that everything is woring correctly!
    if (zScaling_ != ZSCALING) {
      EXCEPTION("Something went wrong: scaling factor in z-direction was not set correctly!" );
    }
    
    // ASSUMPTION: BH-curve of magn. flux in thickness direction (theta=90°) 
    //             is equivalent to BH-curve in transverse direction (phi=90°)
    // !! NOTE !!  till now the maximum angle provided for phi will also  
    //             be chosen for theta as well as its corresponding BH-curve      
    Double coefScalarZ = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs) * zScaling_;

    Double dThetaVal = angles_[kend] - 0.0;
    Double bhi = ( angles_[kend] - angleBTheta ) / dThetaVal;
    Double blo = ( angleBTheta - 0.0 ) / dThetaVal;

    // interpolate between calculated coefficient in xy-plane and z-direction 
    coefScalar = bhi * coefScalarXY + blo * coefScalarZ;


#ifndef NDEBUG
    // some debug output... TODO-avolk: weg damit!
    std::cerr << "GetScalar(): statistics of element " << lpm.ptEl->elemNum << ":";
//    std::cerr << " B_x=" << elemOpSol[0];
//    std::cerr << " B_y=" << elemOpSol[1];
//    std::cerr << " B_z=" << elemOpSol[2];
    std::cerr << " angleBPhi=" << angleBPhi;
    std::cerr << " angleBTheta=" << angleBTheta;
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
  bOperator_ = NULL;
}

//! Destructor
CoefFunctionApproxDerivAniso::~CoefFunctionApproxDerivAniso(){
  delete bOperator_;
  bOperator_ = NULL;
  ;
}

//! Initialize with data
void CoefFunctionApproxDerivAniso::Init( StdVector<ApproxData*>  nLinFnc,
                                         StdVector<Double> angles,
                                         shared_ptr<FeFunction<Double> > fct,
                                         BaseBOperator* bOp ) {
  // set type to TENSOR
  dimType_ = TENSOR;
  nLinFnc_ = nLinFnc;
  angles_ = angles;
  feFct_ = fct;
  bOperator_ = bOp;
  dimDMat_ = bOperator_->GetDimDMat();
  
  // TODO-avolk: this should not be hardcoded here but specified in the xml file!!
  zScaling_ = ZSCALING;
}

void CoefFunctionApproxDerivAniso::GetTensor(Matrix<Double>& coefMat, 
                                             const LocPointMapped& lpm ) {
  
  // extract element solution from feFunction
  Vector<Double> elemSol, elemB;
  feFct_->GetElemSolution(elemSol, lpm.ptEl);

  // apply b-operator matrix to element solution to obtain field value
  BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
  bOperator_->ApplyOp( elemB, lpm, ptFe, elemSol );

  Double fieldAbs = elemB.NormL2();
  coefMat.Resize( dimDMat_, dimDMat_ );
  if( fieldAbs == 0 ) {
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
    
    
    // ---------------------------------
    //  Nearest neighbour interpolation
    // ---------------------------------
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
    
    
    // ---------------------------------------------------
    //  Angular based interpolation in the xy-plane (phi)
    // ---------------------------------------------------
    Double nuPrime, nuPrimeXY, nuPrimeZ;
    
    // if angle is out of bounds or we have just one entry,
    // return boundary value (i.e.first or last) 
    const UInt kend = angles_.GetSize() - 1;
    if ( angleBPhi > angles_[kend] || kend == 0) {
      nuPrimeXY = nLinFnc_[kend]->EvaluatePrimeNu(fieldAbs);
    }
    else if ( angleBPhi < angles_[0] ) {
      nuPrimeXY = nLinFnc_[0]->EvaluatePrimeNu(fieldAbs);
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
      nuPrimeXY =   ahi * nLinFnc_[klo]->EvaluatePrimeNu(fieldAbs) 
                  + alo * nLinFnc_[khi]->EvaluatePrimeNu(fieldAbs);
    }
    
    // --------------------------------------------------------------------
    //  Angular based interpolation according to z direction (angle theta)
    // --------------------------------------------------------------------

    // TODO-avolk: this should be deleted after first tests!!
    // added just to be shure that everything is woring correctly!
    if (zScaling_ != ZSCALING) {
      EXCEPTION("Something went wrong: scaling factor in z-direction was not set correctly!" );
    }
        
    // ASSUMPTION: BH-curve of magn. flux in thickness direction (theta=90°) 
    //             is equivalent to BH-curve in transverse direction (phi=90°)
    // !! NOTE !!  till now the maximum angle provided for phi will also  
    //             be chosen for theta as well as its corresponding BH-curve
    // Note: scaling of mu by factor c leads to scaling of nu' by 1/c since:
    // 1) nu = 1/mu; c*mu => 1/c*mu meaning that scaling has to be applied by the reciprocal factor 
    // 2) g(x)=c*f(x) => g'(x)=c*f'(x) meaning that scaling can also be applied to derivative
    nuPrimeZ = nLinFnc_[kend]->EvaluateFuncNu(fieldAbs) * (1/zScaling_);

    Double dThetaVal = angles_[kend] - 0.0;
    Double bhi = ( angles_[kend] - angleBTheta ) / dThetaVal;
    Double blo = ( angleBTheta - 0.0 ) / dThetaVal;

    // interpolate between calculated nuPrime in xy-plane and z-direction 
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
