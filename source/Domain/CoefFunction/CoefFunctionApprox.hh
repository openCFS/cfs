#ifndef FILE_COEFFUNCTION_APPROX_HH
#define FILE_COEFFUNCTION_APPROX_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField {

// forward clas declaration
class ApproxData;
class BaseBOperator;
class Grid;
class FeFunctions;


// ============================================================================
//  ISOTROPIC VERSIONS
// ============================================================================
//! Provide a coefficient for approximated sample data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApprox : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApprox();

  //! Destructor
  virtual ~CoefFunctionApprox();
  
  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc,
             PtrCoefFct dependCoef  );

  
  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

// ============================================================================
//  Super Approx Coef Function
// ============================================================================
//! Provide a coefficient for approximated sample data (scalar)

//! This class encapsulates a CoefFunctionApprox and holds data on regions where to evaluate
//! the dependCoef
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxSuper : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxSuper();

  //! Destructor
  virtual ~CoefFunctionApproxSuper();
  
  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc,
             PtrCoefFct dependCoef  );

  
  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    
  void SetLumpedRegions(RegionIdType regA, RegionIdType regB, RegionIdType regC = -1);

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:

  //! the actual approximation coef function
  CoefFunctionApprox * ptAproxCoefFct;

  //! region Ids of the terminals
  RegionIdType terminalA_, terminalB_, terminalC_;
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

//! Provide a coefficient for approximated derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxDeriv : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxDeriv();

  //! Destructor
  virtual ~CoefFunctionApproxDeriv();
  
  //! Initialize with data
  void Init( ApproxData * nLinFnc,
             UInt dimDMat,
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
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
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

//! Provide a coefficient for approximated anisotropic derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxDerivAniso : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxDerivAniso();

  //! Destructor
  virtual ~CoefFunctionApproxDerivAniso();
  
  //! Initialize with data
  void Init( StdVector<ApproxData*>  nLinFnc,
             StdVector<Double> angles,
             UInt dimDMat,
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

}
#endif
