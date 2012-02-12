#ifndef FILE_COEFFUNCTION_APPROX_HH
#define FILE_COEFFUNCTION_APPROX_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"

// classes for function / spline approximation
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"

namespace CoupledField {

//! Provide a coefficient for approximated sample data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
template<typename B_OP>
class CoefFunctionApprox : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApprox() {
    // this type of coefficient is nonlinear (i.e. solution dependend)
    dependType_ = SOLUTION;
  }

  //! Destructor
  virtual ~CoefFunctionApprox(){
    ;
  }
  
  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc,
             shared_ptr<FeFunction<Double> > fct ) {
    
    // set type to scalar
    dimType_ = SCALAR;
    coefScalar_ = coefScalar;
    nLinFnc_ = nLinFnc;
    feFct_ = fct;
  }

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm ) {
    
    // extract element solution from feFunction
    Vector<Double> elemSol, elemB;
    feFct_->GetElemSolution(elemSol, lpm.ptEl);
    
    // apply b-operator matrix to element solution to obtain field value
    BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    bOperator_.ApplyOp( elemB, lpm, ptFe, elemSol );
    
    Double fieldAbs = elemB.NormL2();
    
     if( fieldAbs == 0 ) { 
        coefScalar = coefScalar_;
     } else {
        coefScalar = nLinFnc_->EvaluateFuncNu(fieldAbs);
     }
  }

  //! \see CoefFunction::IsComplex
  bool IsComplex(){
    return false;
  }
  
  //! \see CoefFunction::ToString
  std::string ToString() const {
    EXCEPTION( "Implement me");
  }

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Differential operator to calculate the field value
  B_OP bOperator_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Depending FeFunction
  shared_ptr<FeFunction<Double> > feFct_; 
};

//! Provide a coefficient for approximated derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
template<typename B_OP>
class CoefFunctionApproxDeriv : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxDeriv() {
    // this type of coefficient is nonlinear, i.e. spatial and time dependent
    dependType_ = GENERAL;
  }

  //! Destructor
  virtual ~CoefFunctionApproxDeriv(){
    ;
  }
  
  //! Initialize with data
  void Init( ApproxData * nLinFnc,
             shared_ptr<FeFunction<Double> > fct ) {
    
    // set type to TENSOR
    dimType_ = TENSOR;
    nLinFnc_ = nLinFnc;
    feFct_ = fct;
  }

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm ) {


    // extract element solution from feFunction
    Vector<Double> elemSol, elemB;
    feFct_->GetElemSolution(elemSol, lpm.ptEl);

    // apply b-operator matrix to element solution to obtain field value
    BaseFE * ptFe = feFct_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    bOperator_.ApplyOp( elemB, lpm, ptFe, elemSol );

    Double fieldAbs = elemB.NormL2();
    coefMat.Resize( B_OP::DIM_D_MAT, B_OP::DIM_D_MAT );
    if( fieldAbs == 0 ) {
      coefMat.Init();
    } else {
      
      // evaluate derivative of reluctivity
      Double nuPrime = nLinFnc_->EvaluatePrimeNu(fieldAbs);
      // coefMat = B^T [ e_B^T * nu' * |B| * e_B] B 
      Vector<Double> eB(B_OP::DIM_D_MAT);
      eB = elemB / fieldAbs;
      coefMat.DyadicMult( eB );
      coefMat *= nuPrime * fieldAbs;
    }
  }

  //! \see CoefFunction::IsComplex
  bool IsComplex(){
    return false;
  }
  
  //! \see CoefFunction::ToString
  std::string ToString() const {
  EXCEPTION( "Implement me");
  }

protected:
  
  //! Differential operator to calculate the field value
  B_OP bOperator_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Depending FeFunction
  shared_ptr<FeFunction<Double> > feFct_; 
};



}
#endif
