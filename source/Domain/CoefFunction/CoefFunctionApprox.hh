#ifndef FILE_COEFFUNCTION_APPROX_HH
#define FILE_COEFFUNCTION_APPROX_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField {

// forward class declaration
class ApproxData;
class BaseBOperator;


// ============================================================================
//  ISOTROPIC VERSIONS
// ============================================================================
//! Provide a coefficient for approximated sample data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApprox : public CoefFunction {
public:

  //! Constructor
  CoefFunctionApprox();

  //! Destructor
  virtual ~CoefFunctionApprox();
  
  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc,
             shared_ptr<FeFunction<Double> > fct,
             BaseBOperator* bOp );

  
  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Differential operator to calculate the field value
  BaseBOperator* bOperator_;
  
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
class CoefFunctionApproxDeriv : public CoefFunction {
public:

  //! Constructor
  CoefFunctionApproxDeriv();

  //! Destructor
  virtual ~CoefFunctionApproxDeriv();
  
  //! Initialize with data
  void Init( ApproxData * nLinFnc,
             shared_ptr<FeFunction<Double> > fct,
             BaseBOperator* bOp );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Differential operator to calculate the field value
  BaseBOperator* bOperator_;
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Depending FeFunction
  shared_ptr<FeFunction<Double> > feFct_; 
};

// ============================================================================
//  ANISOTROPIC VERSIONS
// ============================================================================
//! Provide a coefficient for approximated anisotropic sampled data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! In addition it maps the 
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxAniso : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxAniso();

  //! Destructor
  virtual ~CoefFunctionApproxAniso();

  //! Initialize with data
  void Init( Double coefScalar, 
             StdVector<ApproxData*>  nLinFnc,
             StdVector<Double> angles,
             StdVector<Double> zScalings,
             shared_ptr<FeFunction<Double> > fct,
             BaseBOperator* bOp );

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Differential operator to calculate the field value
  BaseBOperator* bOperator_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Scaling factor of anisotropic behavior in z-direction
  //! -> Since we do not yet have nonlinear curves in z direction we use the same 
  //! curves as given for the xy-plane but scale it with an appropriate factor.
  //! Scaling is meant to be applied to mu (provided via BH-curve) -> nu is scaled by 1/zScaling_
  StdVector<Double> zScalings_;
  
  //! Depending FeFunction
  shared_ptr<FeFunction<Double> > feFct_; 
};

//! Provide a coefficient for approximated anisotropic derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxDerivAniso : public CoefFunction {
public:

  //! Constructor
  CoefFunctionApproxDerivAniso();

  //! Destructor
  virtual ~CoefFunctionApproxDerivAniso();
  
  //! Initialize with data
  void Init( StdVector<ApproxData*>  nLinFnc,
             StdVector<Double> angles,
             StdVector<Double> zScalings,
             shared_ptr<FeFunction<Double> > fct,
             BaseBOperator* bOp );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Differential operator to calculate the field value
  BaseBOperator* bOperator_;
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Scaling factor of anisotropic behavior in z-direction
  //! -> Since we do not yet have nonlinear curves in z direction we use the same 
  //! curves as given for the xy-plane but scale it with an appropriate factor.
  //! Scaling is meant to be applied to mu (provided via BH-curve) -> nu is scaled by 1/zScaling_
  StdVector<Double> zScalings_;
  
  //! Depending FeFunction
  shared_ptr<FeFunction<Double> > feFct_; 
};

}
#endif
