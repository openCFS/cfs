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


    //compute angle 
    Double angleB;
    if ( abs(elemOpSol[0]) > 1e-5 ) {
      angleB = abs( std::atan( elemOpSol[1] / elemOpSol[0] ) );
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
